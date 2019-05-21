//
// Created by tommaso on 16/05/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH

//#ifdef __NVCC__

#include <cstdlib>
#include <crt/host_defines.h>
#include <device_launch_parameters.h>
#include "SparseGridGpu.hpp"

namespace SparseGridGpuKernels
{
    template<unsigned int p, unsigned int maskProp, typename InsertBufferT, typename ScalarT>
    __global__ void initializeInsertBuffer(InsertBufferT insertBuffer, ScalarT backgroundValue)
    {
        typedef typename InsertBufferT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, p> BlockT;

        int pos = blockIdx.x * blockDim.x + threadIdx.x;
        unsigned int dataBlockId = pos / BlockT::size;
        unsigned int offset = pos % BlockT::size;

        insertBuffer.template get<maskProp>(dataBlockId)[offset] = 0;
        insertBuffer.template get<p>(dataBlockId)[offset] = backgroundValue;
    }

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
    
    /**
     * Compact various memory pools (with empty segments) into a single contiguous memory region.
     * NOTE: Each thread block is in charge of one pool
     * 
     * @tparam p The index of the property of the aggregate type to be considered.
     * @tparam AggregateT The aggregate type.
     * @tparam IndexT The type of the indices.
     * @param src The pointer to the start of the first pool (pools are contiguous).
     * @param blockPoolSize The size of each pool (number of p-th property objects in each pool).
     * @param starts The exclusive scan of actually used number of elements in each pool.
     * @param dst The pointer to the start of the destination memory region.
     */
    template <typename DataVectorT, typename IndexVectorT>
    __global__ void compact(DataVectorT & src, size_t blockPoolSize, IndexVectorT & starts, DataVectorT & dst)
    {
        typedef typename DataVectorT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, 0> BlockT0; // The type of the 0-th property
        unsigned int chunkSize = BlockT0::size;

        unsigned int poolId = blockIdx.x; // Each block is in charge of one pool
        unsigned int poolStartPos = poolId * blockPoolSize; // Pos at which pool starts in src
        unsigned int chunksPerBlock = blockDim.x / chunkSize;
        unsigned int numChunksToProcess = starts[poolId+1] + starts[poolId];
        unsigned int chunkOffset = threadIdx.x / chunkSize; // The thread block can work on several chunks in parallel
//        unsigned int elementOffset = threadIdx.x % chunkSize; // Each thread gets one element of a chunk to work on

        unsigned int curChunk = 0;
        for (curChunk; curChunk < numChunksToProcess - chunksPerBlock; curChunk += chunksPerBlock)
        {
            unsigned int dstId = starts[poolId] + curChunk + chunkOffset;
            unsigned int srcId = poolStartPos + curChunk + chunkOffset;
            dst.get(dstId) = src.get(srcId);
        }
        if (curChunk + chunkOffset < numChunksToProcess)
        {
            unsigned int dstId = starts[poolId] + curChunk + chunkOffset;
            unsigned int srcId = poolStartPos + curChunk + chunkOffset;
            dst.get(dstId) = src.get(srcId);
        }
    }
    
    template <typename IndexT>
    __global__ void reorder(IndexT* dstIndices)
    {
        //todo
    }
    
    // Below the kernels to be used inside the "compute segments" part of the solveConflicts
    template <typename IndexT>
    __global__ void computePredicates(IndexT* keys, IndexT* predicates)
    {
        //todo
    }
    
    template <typename IndexT>
    __global__ void copyPositionToDestIfPredicate(IndexT* predicates, IndexT* exclusiveScanOfPredicates)
    {
        //todo
    }
    
    //todo: Build the functors for the vector sparse flush() in some place
}

//#endif //__NVCC__

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
