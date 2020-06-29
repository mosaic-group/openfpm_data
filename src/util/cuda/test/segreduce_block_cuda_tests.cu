//
// Created by tommaso on 20/05/19.
//

#include "config.h"

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "Vector/map_vector.hpp"
#include "SparseGridGpu/BlockMapGpu.hpp"
//todo: here include SparseGridGpu_kernels and remove local kernel definitions

template<typename op, typename ScalarT>
__device__ inline ScalarT applyOp(ScalarT a, ScalarT b, bool aExist, bool bExist)
{
    op op_;
    if (aExist && bExist)
    {
        a = op_(a, b);
    }
    else if (bExist)
    {
        a = b;
    }
    return a;
}

// GridSize = number of segments
// BlockSize = chunksPerBlock * chunkSize
//
template<unsigned int chunksPerBlock, typename op, typename SegType, typename DataType, typename MaskType>
__global__ void
segreduce_block(
        DataType *data,
        SegType *segments,
        MaskType *masks,
        DataType *output
        )
{
    unsigned int segmentId = blockIdx.x;
    int segmentSize = segments[segmentId + 1] - segments[segmentId];

    unsigned int start = segments[segmentId];

    unsigned int chunkId = threadIdx.x / DataType::size;
    unsigned int offset = threadIdx.x % DataType::size;

    __shared__ DataType A[chunksPerBlock];
    __shared__ MaskType AMask[chunksPerBlock];
    typename DataType::scalarType bReg;
    typename MaskType::scalarType aMask, bMask;

    // Phase 0: Load chunks as first operand of the reduction
    if (chunkId < segmentSize)
    {
        A[chunkId][offset] = data[start + chunkId][offset];
        aMask = masks[start + chunkId][offset];
    }

    int i = chunksPerBlock;
    for ( ; i < segmentSize - (int) (chunksPerBlock); i += chunksPerBlock)
    {
        bReg = data[start + i + chunkId][offset];
        bMask = masks[start + i + chunkId][offset];

        A[chunkId][offset] = applyOp<op>(A[chunkId][offset],
                                         bReg,
                                         BlockMapGpu_ker<>::exist(aMask),
                                         BlockMapGpu_ker<>::exist(bMask));
        aMask = aMask | bMask;
    }

    if (i + chunkId < segmentSize)
    {
        bReg = data[start + i + chunkId][offset];
        bMask = masks[start + i + chunkId][offset];

        A[chunkId][offset] = applyOp<op>(A[chunkId][offset],
                                         bReg,
                                         BlockMapGpu_ker<>::exist(aMask),
                                         BlockMapGpu_ker<>::exist(bMask));
        aMask = aMask | bMask;
    }

    AMask[chunkId][offset] = aMask;

    __syncthreads();

    // Horizontal reduction finished
    // Now vertical reduction
    for (int j = 2; j <= chunksPerBlock && j <= segmentSize; j *= 2)
    {
        if (chunkId % j == 0 && chunkId < segmentSize)
        {
            unsigned int otherChunkId = chunkId + (j / 2);
            if (otherChunkId < segmentSize)
            {
                aMask = AMask[chunkId][offset];
                bMask = AMask[otherChunkId][offset];
                A[chunkId][offset] = applyOp<op>(A[chunkId][offset],
                                                 A[otherChunkId][offset],
                                                 BlockMapGpu_ker<>::exist(aMask),
                                                 BlockMapGpu_ker<>::exist(bMask));
                AMask[chunkId][offset] = aMask | bMask;
            }
        }
        __syncthreads();
    }

    // Write output
    if (chunkId == 0)
    {
        output[segmentId][offset] = A[chunkId][offset];
    }
}

BOOST_AUTO_TEST_SUITE(segreduce_block_cuda_tests)

    BOOST_AUTO_TEST_CASE (segreduce_block_test)
    {
        typedef float ScalarT;
        typedef DataBlock<unsigned char, 64> MaskBlockT;
        typedef DataBlock<ScalarT, 64> BlockT;
        openfpm::vector_gpu<aggregate<int>> segments;
        segments.resize(8);
        segments.template get<0>(0) = 0;
        segments.template get<0>(1) = 4;
        segments.template get<0>(2) = 5;
        segments.template get<0>(3) = 7;
        segments.template get<0>(4) = 8;
        segments.template get<0>(5) = 11;
        segments.template get<0>(6) = 17;
        segments.template get<0>(7) = 18; // Id of first non-existent data

        segments.template hostToDevice<0>();

        const unsigned int BITMASK = 0, BLOCK = 1;
        BlockT block;
        MaskBlockT mask;
        for (int i = 0; i < 32; ++i)
        {
            block[i] = i + 1;
            mask[i] = 1;
        }
        for (int i = 32; i < 64; ++i)
        {
            block[i] = 666;
            mask[i] = 0;
        }

        openfpm::vector_gpu<aggregate<MaskBlockT, BlockT>> data;
        data.resize(18);
        for (int i = 0; i < 18; ++i)
        {
            data.template get<BITMASK>(i) = mask;
            data.template get<BLOCK>(i) = block;
        }

        data.template hostToDevice<BITMASK, BLOCK>();

        // Allocate output buffer
        openfpm::vector_gpu<aggregate<MaskBlockT, BlockT>> outputData;
        outputData.resize(segments.size()-1);

        // template<unsigned int chunksPerBlock, typename op, typename SegType, typename DataType, typename MaskType>
        // segreduce(DataType *data, SegType *segments, MaskType *masks, DataType *output, MaskType *outputMasks)
//        segreduce<2, mgpu::maximum_t<ScalarT>> <<< outputData.size(), 2*BlockT::size >>> (
        segreduce_block<2, mgpu::plus_t<ScalarT>> <<< outputData.size(), 2*BlockT::size >>> (
                        (BlockT *) data.template getDeviceBuffer<BLOCK>(),
                        (int *) segments.template getDeviceBuffer<0>(),
                        (MaskBlockT *) data.template getDeviceBuffer<BITMASK>(),
                        (BlockT *) outputData.template getDeviceBuffer<BLOCK>()
                                );

        // Segreduce on mask
        segreduce_block<2, mgpu::maximum_t<unsigned char>> <<< outputData.size(), 2*BlockT::size >>> (
                        (MaskBlockT *) data.template getDeviceBuffer<BITMASK>(),
                        (int *) segments.template getDeviceBuffer<0>(),
                        (MaskBlockT *) data.template getDeviceBuffer<BITMASK>(),
                        (MaskBlockT *) outputData.template getDeviceBuffer<BITMASK>()
        );

        outputData.template deviceToHost<BITMASK, BLOCK>();

        for (int j = 0; j < outputData.size(); ++j)
        {
            BlockT outBlock = outputData.template get<BLOCK>(j);
            MaskBlockT outMask = outputData.template get<BITMASK>(j);

            int seg_sz = segments.template get<0>(j+1) - segments.template get<0>(j);

            for (int i = 0; i < BlockT::size; ++i)
            {
                if (i < 32)
                {
                	BOOST_REQUIRE_EQUAL(outBlock[i],seg_sz*(i+1));
                	BOOST_REQUIRE_EQUAL(outMask[i],1);
                }
                else
                {
                	BOOST_REQUIRE_EQUAL(outMask[i],0);
                }
            }
        }
    }

BOOST_AUTO_TEST_SUITE_END()
