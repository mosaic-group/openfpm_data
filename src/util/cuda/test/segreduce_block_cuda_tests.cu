//
// Created by tommaso on 20/05/19.
//

#include "config.h"

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "Vector/map_vector.hpp"
#include "SparseGridGpu/SparseGridGpu.hpp"

template<unsigned int DataBlockSize, typename ScalarT>
struct DataBlockLoc
{
    typedef ScalarT scalarType;

    static const unsigned int size = DataBlockSize;
    ScalarT block[size];

    __device__ __host__ DataBlockLoc() = default;

    __device__ __host__ DataBlockLoc(const DataBlockLoc &other)
    {
        memcpy(block, other.block, size * sizeof(ScalarT));
    }

    __device__ __host__ DataBlockLoc operator=(const DataBlockLoc &other)
    {
        memcpy(block, other.block, size * sizeof(ScalarT));
        return *this;
    }

    __device__ __host__ DataBlockLoc operator=(float v) // Hack to make things compile
    {
        return *this;
    }

    __device__ __host__ inline ScalarT &operator[](unsigned int i)
    {
        return block[i];
    }

    __device__ __host__ inline const ScalarT &operator[](unsigned int i) const
    {
        return block[i];
    }

    __device__ __host__ inline static bool exist(unsigned int i, size_t &bitMask)
    {
        return (bitMask >> i) & ((size_t) 1);
    }

    __device__ __host__ inline static void setElement(unsigned int i, size_t &bitMask)
    {
        bitMask = bitMask | ((size_t) 1) << i;
    }
};

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
        DataType *output,
        MaskType *outputMasks
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
    MaskType aMask, bMask;

    // Phase 0: Load chunks as first operand of the reduction
    if (chunkId < segmentSize)
    {
        A[chunkId][offset] = data[start + chunkId][offset];
        aMask = masks[start + chunkId];
    }

    int i = chunksPerBlock;
    for ( ; i < segmentSize - (int) (chunksPerBlock); i += chunksPerBlock)
    {
        bReg = data[start + i + chunkId][offset];
        bMask = masks[start + i + chunkId];

        A[chunkId][offset] = applyOp<op>(A[chunkId][offset],
                                         bReg,
                                         DataType::exist(offset, aMask),
                                         DataType::exist(offset, bMask));
        aMask = aMask | bMask;
    }

    if (i + chunkId < segmentSize)
    {
        bReg = data[start + i + chunkId][offset];
        bMask = masks[start + i + chunkId];

        A[chunkId][offset] = applyOp<op>(A[chunkId][offset],
                                         bReg,
                                         DataType::exist(offset, aMask),
                                         DataType::exist(offset, bMask));
        aMask = aMask | bMask;
    }

    if (offset == 0) // && chunkId < segmentSize)  //todo: check if second condition is necessary
    {
        AMask[chunkId] = aMask;
    }

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
                aMask = AMask[chunkId];
                bMask = AMask[otherChunkId];
                A[chunkId][offset] = applyOp<op>(A[chunkId][offset],
                                                 A[otherChunkId][offset],
                                                 DataType::exist(offset, aMask),
                                                 DataType::exist(offset, bMask));
                AMask[chunkId] = aMask | bMask;
            }
        }
        __syncthreads();
    }

    // Write output
    if (chunkId == 0)
    {
        output[segmentId][offset] = A[chunkId][offset];
        if (offset == 0)
        {
            outputMasks[segmentId] = AMask[chunkId];
        }
    }
}

BOOST_AUTO_TEST_SUITE(segreduce_block_cuda_tests)

    BOOST_AUTO_TEST_CASE (segreduce_block_test)
    {
        typedef float ScalarT;
        typedef DataBlockLoc<64, ScalarT> BlockT;
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
        const size_t mask = 0xFFFFFFFF;
        BlockT block;
        for (int i = 0; i < 32; ++i)
        {
            block[i] = i + 1;
        }
        for (int i = 32; i < 64; ++i)
        {
            block[i] = 666;
        }

        openfpm::vector_gpu<aggregate<size_t, BlockT>> data;
        data.resize(18);
        for (int i = 0; i < 18; ++i)
        {
            data.template get<BITMASK>(i) = mask;
            data.template get<BLOCK>(i) = block;
        }

        data.template hostToDevice<BITMASK, BLOCK>();

        // Allocate output buffer
        openfpm::vector_gpu<aggregate<size_t, BlockT>> outputData;
        outputData.resize(segments.size()-1);

        // template<unsigned int chunksPerBlock, typename op, typename SegType, typename DataType, typename MaskType>
        // segreduce_block(DataType *data, SegType *segments, MaskType *masks, DataType *output, MaskType *outputMasks)
//        segreduce_block<2, mgpu::maximum_t<ScalarT>> <<< outputData.size(), 2*BlockT::size >>> (
        segreduce_block<2, mgpu::plus_t<ScalarT>> <<< outputData.size(), 2*BlockT::size >>> (
                        (BlockT *) data.template getDeviceBuffer<BLOCK>(),
                        (int *) segments.template getDeviceBuffer<0>(),
                        (size_t *) data.template getDeviceBuffer<BITMASK>(),
                        (BlockT *) outputData.template getDeviceBuffer<BLOCK>(),
                        (size_t *) outputData.template getDeviceBuffer<BITMASK>()
                                );

        outputData.template deviceToHost<BITMASK, BLOCK>();

        for (int j = 0; j < outputData.size(); ++j)
        {
            BlockT outBlock = outputData.template get<BLOCK>(j);
            size_t outMask = outputData.template get<BITMASK>(j);

            std::cout << std::bitset<64>(outMask) << std::endl;
            for (int i = 0; i < BlockT::size; ++i)
            {
                std::cout << outBlock[i] << " ";
            }
            std::cout << std::endl;
        }
    }

BOOST_AUTO_TEST_SUITE_END()