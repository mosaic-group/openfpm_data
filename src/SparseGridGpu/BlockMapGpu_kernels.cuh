//
// Created by tommaso on 16/05/19.
//

#ifndef OPENFPM_PDATA_BLOCKMAPGPU_KERNELS_CUH
#define OPENFPM_PDATA_BLOCKMAPGPU_KERNELS_CUH

//#ifdef __NVCC__

#include <cstdlib>
#include "BlockMapGpu.hpp"
#include "util/cuda_util.hpp"
#include "BlockMapGpu_dimensionalityWrappers.cuh"
#include "Vector/map_vector_sparse.hpp"

namespace BlockMapGpuKernels
{
    template<unsigned int maskProp, unsigned int chunksPerBlock, typename InsertBufferT>
    __global__ void initializeInsertBuffer(InsertBufferT insertBuffer)
    {
#ifdef __NVCC__
        typedef typename InsertBufferT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, maskProp> MaskT;

        int pos = blockIdx.x * blockDim.x + threadIdx.x;
        const unsigned int dataBlockId = pos / MaskT::size;
        const unsigned int offset = pos % MaskT::size;
        const unsigned int chunkOffset = dataBlockId % chunksPerBlock;

        __shared__ MaskT mask[chunksPerBlock];

        if (dataBlockId < insertBuffer.size())
        {
            mask[chunkOffset][offset] = 0;

            __syncthreads();

            // Here the operator= spreads nicely across threads...
            insertBuffer.template get<maskProp>(dataBlockId) = mask[chunkOffset];
        }
        else
        {
            __syncthreads();
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }

    /**
     * Apply a binary operator on 2 operands and put the result on the first one.
     *
     * @tparam op The binary operator functor.
     * @tparam ScalarT The type accepted by the operator.
     * @param a First operand. It will be filled with the result.
     * @param b Second operand.
     * @param aExist If the first operand's value is valid.
     * @param bExist If the second operand's value is valid.
     */
    template<typename op, typename ScalarT>
    __device__ inline void applyOp(ScalarT &a, ScalarT b, bool aExist, bool bExist)
    {
#ifdef __NVCC__
        op op_;
        if (aExist && bExist)
        {
            a = op_(a, b);
        }
        else if (bExist)
        {
            a = b;
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }

    // GridSize = number of segments
    // BlockSize = chunksPerBlock * chunkSize
    //
    template<unsigned int p,
            unsigned int pSegment,
            unsigned int pMask,
            unsigned int chunksPerBlock,
            typename op,
            typename IndexVectorT, typename DataVectorT, typename MaskVectorT>
    __global__ void
    segreduce(
            DataVectorT data,
            IndexVectorT segments,
            MaskVectorT masks,
            DataVectorT output
    )
    {
#ifdef __NVCC__
        typedef typename DataVectorT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, p> DataType;
        typedef BlockTypeOf<AggregateT, pMask> MaskType;
        typedef typename std::remove_all_extents<DataType>::type BaseBlockType;
        constexpr unsigned int chunkSize = BaseBlockType::size;

        unsigned int segmentId = blockIdx.x;
        int segmentSize = segments.template get<pSegment>(segmentId + 1)
                          - segments.template get<pSegment>(segmentId);

        unsigned int start = segments.template get<pSegment>(segmentId);

        unsigned int chunkId = threadIdx.x / chunkSize;
        unsigned int offset = threadIdx.x % chunkSize;

        __shared__ ArrayWrapper<DataType> A[chunksPerBlock];
        __shared__ MaskType AMask[chunksPerBlock];
        typename ComposeArrayType<DataType>::type bReg;
        typename MaskType::scalarType aMask, bMask;

        // Phase 0: Load chunks as first operand of the reduction
        if (chunkId < segmentSize)
        {
            A[chunkId][offset] = RhsBlockWrapper<DataType>(data.template get<p>(start + chunkId), offset).value;
            aMask = masks.template get<pMask>(start + chunkId)[offset];
        }

        int i = chunksPerBlock;
        for (; i < segmentSize - (int) (chunksPerBlock); i += chunksPerBlock)
        {
            generalDimensionFunctor<decltype(bReg)>::assignWithOffsetRHS(bReg,
                                                                         data.template get<p>(start + i + chunkId),
                                                                         offset);
            bMask = masks.template get<pMask>(start + i + chunkId)[offset];

            generalDimensionFunctor<DataType>::template applyOp<op>(A[chunkId][offset],
                                                                    bReg,
                                                                    BlockMapGpu_ker<>::exist(aMask),
                                                                    BlockMapGpu_ker<>::exist(bMask));
            aMask = aMask | bMask;
        }

        if (i + chunkId < segmentSize)
        {
            generalDimensionFunctor<decltype(bReg)>::assignWithOffsetRHS(bReg,
                                                                         data.template get<p>(start + i + chunkId),
                                                                         offset);
            bMask = masks.template get<pMask>(start + i + chunkId)[offset];

            generalDimensionFunctor<DataType>::template applyOp<op>(A[chunkId][offset],
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
                    generalDimensionFunctor<DataType>::template applyOp<op>(A[chunkId][offset],
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
            generalDimensionFunctor<DataType>::assignWithOffset(output.template get<p>(segmentId), A[chunkId].data,
                                                                offset);
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }

    // GridSize = number of segments
    // BlockSize = chunksPerBlock * chunkSize
    //
    template<unsigned int p,
            unsigned int pSegment,
            unsigned int pMask,
            unsigned int chunksPerBlock,
            typename op,
            typename IndexVectorT, typename IndexVectorT2, typename DataVectorT, typename MaskVectorT>
    __global__ void
    segreduce_beta(
            DataVectorT data_new,
            IndexVectorT2 segments_data,
            IndexVectorT segments_dataMap,
            MaskVectorT masks,
            DataVectorT output
    )
    {
#ifdef __NVCC__
        typedef typename DataVectorT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, p> DataType;
        typedef BlockTypeOf<AggregateT, pMask> MaskType;
        typedef typename std::remove_all_extents<DataType>::type BaseBlockType;
        constexpr unsigned int chunkSize = BaseBlockType::size;

        unsigned int segmentId = blockIdx.x;
        int segmentSize = segments_data.template get<pSegment>(segmentId + 1)
                          - segments_data.template get<pSegment>(segmentId);

        unsigned int start = segments_data.template get<pSegment>(segmentId);

        unsigned int chunkId = threadIdx.x / chunkSize;
        unsigned int offset = threadIdx.x % chunkSize;

        __shared__ ArrayWrapper<DataType> A[chunksPerBlock];
        __shared__ MaskType AMask[chunksPerBlock];
        typename ComposeArrayType<DataType>::type bReg;
        typename MaskType::scalarType aMask, bMask;

        // Phase 0: Load chunks as first operand of the reduction
        if (chunkId < segmentSize)
        {
        	unsigned int m_chunkId = segments_dataMap.template get<0>(start + chunkId);
            A[chunkId][offset] = RhsBlockWrapper<DataType>(data_new.template get<p>(m_chunkId), offset).value;
            aMask = data_new.template get<pMask>(m_chunkId)[offset];
        }

        //////////////////////////// Horizontal reduction (work on data) //////////////////////////////////////////////////
        //////////////////////////// We reduce

        int i = chunksPerBlock;
        for (; i < segmentSize - (int) (chunksPerBlock); i += chunksPerBlock)
        {
        	unsigned int m_chunkId = segments_dataMap.template get<0>(start + i + chunkId);

        	// it breg = data_new.template get<p>(m_chunkId)
            generalDimensionFunctor<decltype(bReg)>::assignWithOffsetRHS(bReg,
                                                                         data_new.template get<p>(m_chunkId),
                                                                         offset);
            bMask = data_new.template get<pMask>(m_chunkId)[offset];

            // it reduce A[chunkId][offset] with breg
            generalDimensionFunctor<DataType>::template applyOp<op>(A[chunkId][offset],
                                                                    bReg,
                                                                    BlockMapGpu_ker<>::exist(aMask),
                                                                    BlockMapGpu_ker<>::exist(bMask));
            // reduce aMask with bMask
            aMask = aMask | bMask;
        }


        if (i + chunkId < segmentSize)
        {
        	unsigned int m_chunkId = segments_dataMap.template get<0>(start + i + chunkId);

        	// it breg = data_new.template get<p>(m_chunkId)
            generalDimensionFunctor<decltype(bReg)>::assignWithOffsetRHS(bReg,
                                                                         data_new.template get<p>(m_chunkId),
                                                                         offset);
            bMask = data_new.template get<pMask>(m_chunkId)[offset];

            // it reduce A[chunkId][offset] with breg
            generalDimensionFunctor<DataType>::template applyOp<op>(A[chunkId][offset],
                                                                    bReg,
                                                                    BlockMapGpu_ker<>::exist(aMask),
                                                                    BlockMapGpu_ker<>::exist(bMask));

            // reduce aMask with bMask
            aMask = aMask | bMask;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
                    generalDimensionFunctor<DataType>::template applyOp<op>(A[chunkId][offset],
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
            generalDimensionFunctor<DataType>::assignWithOffset(output.template get<p>(segmentId), A[chunkId].data,
                                                                offset);
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }


    // GridSize = number of segments
    // BlockSize = chunksPerBlock * chunkSize
    //
    template<unsigned int p,
            unsigned int pSegment,
            unsigned int pMask,
            unsigned int chunksPerBlock,
            typename op,
            typename IndexVectorT, typename IndexVectorT2, typename DataVectorT>
    __global__ void
    segreduce_total(
            DataVectorT data_new,
            DataVectorT data_old,
            IndexVectorT2 segments_data,
            IndexVectorT segments_dataMap,
            IndexVectorT segments_dataOld,
            IndexVectorT outputMap,
            DataVectorT output
    )
    {
#ifdef __NVCC__
        typedef typename DataVectorT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, p> DataType;
        typedef BlockTypeOf<AggregateT, pMask> MaskType;
        typedef typename std::remove_all_extents<DataType>::type BaseBlockType;
        constexpr unsigned int chunkSize = BaseBlockType::size;

        unsigned int segmentId = blockIdx.x;
        int segmentSize = segments_data.template get<pSegment>(segmentId + 1)
                          - segments_data.template get<pSegment>(segmentId);

        unsigned int start = segments_data.template get<pSegment>(segmentId);

        unsigned int chunkId = threadIdx.x / chunkSize;
        unsigned int offset = threadIdx.x % chunkSize;

        __shared__ ArrayWrapper<DataType> A[chunksPerBlock];
        __shared__ MaskType AMask[chunksPerBlock];
        typename ComposeArrayType<DataType>::type bReg;
        typename MaskType::scalarType aMask, bMask;

        // Phase 0: Load chunks as first operand of the reduction
        if (chunkId < segmentSize)
        {
        	unsigned int m_chunkId = segments_dataMap.template get<0>(start + chunkId);
            A[chunkId][offset] = RhsBlockWrapper<DataType>(data_new.template get<p>(m_chunkId), offset).value;
            aMask = data_new.template get<pMask>(m_chunkId)[offset];
        }

        //////////////////////////// Horizontal reduction (work on data) //////////////////////////////////////////////////
        //////////////////////////// We reduce

        int i = chunksPerBlock;
        for (; i < segmentSize - (int) (chunksPerBlock); i += chunksPerBlock)
        {
        	unsigned int m_chunkId = segments_dataMap.template get<0>(start + i + chunkId);

        	// it breg = data_new.template get<p>(m_chunkId)
            generalDimensionFunctor<decltype(bReg)>::assignWithOffsetRHS(bReg,
                                                                         data_new.template get<p>(m_chunkId),
                                                                         offset);
            bMask = data_new.template get<pMask>(m_chunkId)[offset];

            // it reduce A[chunkId][offset] with breg
            generalDimensionFunctor<DataType>::template applyOp<op>(A[chunkId][offset],
                                                                    bReg,
                                                                    BlockMapGpu_ker<>::exist(aMask),
                                                                    BlockMapGpu_ker<>::exist(bMask));
            // reduce aMask with bMask
            aMask = aMask | bMask;
        }


        if (i + chunkId < segmentSize)
        {
        	unsigned int m_chunkId = segments_dataMap.template get<0>(start + i + chunkId);

        	// it breg = data_new.template get<p>(m_chunkId)
            generalDimensionFunctor<decltype(bReg)>::assignWithOffsetRHS(bReg,
                                                                         data_new.template get<p>(m_chunkId),
                                                                         offset);
            bMask = data_new.template get<pMask>(m_chunkId)[offset];

            // it reduce A[chunkId][offset] with breg
            generalDimensionFunctor<DataType>::template applyOp<op>(A[chunkId][offset],
                                                                    bReg,
                                                                    BlockMapGpu_ker<>::exist(aMask),
                                                                    BlockMapGpu_ker<>::exist(bMask));

            // reduce aMask with bMask
            aMask = aMask | bMask;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
                    generalDimensionFunctor<DataType>::template applyOp<op>(A[chunkId][offset],
                                                                            A[otherChunkId][offset],
                                                                            BlockMapGpu_ker<>::exist(aMask),
                                                                            BlockMapGpu_ker<>::exist(bMask));
                    AMask[chunkId][offset] = aMask | bMask;
                }
            }
            __syncthreads();
        }

        //////////////////////////////////////// Reduce now with old data if present link

        int seg_old = segments_dataOld.template get<0>(segmentId);

        if (seg_old != -1)
        {
        	aMask = AMask[0][offset];
        	bMask = data_old.template get<pMask>(seg_old)[offset];
        	generalDimensionFunctor<DataType>::template applyOp<op>(A[0][offset],
        														data_old.template get<p>(seg_old)[offset],
                                                                BlockMapGpu_ker<>::exist(aMask),
                                                                BlockMapGpu_ker<>::exist(bMask));
        	AMask[0][offset] = aMask | bMask;
        }

        __syncthreads();

        ///////////////////////////////////////////////////////////////////////////////////

        // Write output
        if (chunkId == 0)
        {
        	unsigned int out_id = outputMap.template get<0>(segmentId);
            generalDimensionFunctor<DataType>::assignWithOffset(output.template get<p>(out_id), A[chunkId].data,
                                                                offset);
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }

    /**
     * Reorder blocks of data according to the permutation given by the input srcIndices vector.
     * NOTE: Each thread block is in charge of a fixed amount (automatically determined) of data blocks.
     *
     * @tparam DataVectorT The type of the (OpenFPM) vector of values.
     * @tparam IndexVectorT The type of the (OpenFPM) vector of indices.
     * @param src The input vector of data, which needs to be reordered.
     * @param srcIndices Vector which specifies the reordering permutation, i.e. dst[i] = src[srcIndices[i]].
     * @param dst The output vector of reordered data.
     */
    template<typename DataVectorT, typename IndexVectorT>
    __global__ void reorder(DataVectorT src, IndexVectorT srcIndices, DataVectorT dst)
    {
#ifdef __NVCC__
        typedef typename DataVectorT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, 0> BlockT0; // The type of the 0-th property
        unsigned int chunkSize = BlockT0::size;

        unsigned int chunksPerBlock = blockDim.x / chunkSize;
        unsigned int chunkOffset = threadIdx.x / chunkSize; // The thread block can work on several chunks in parallel

        unsigned int chunkBasePos = blockIdx.x * chunksPerBlock;

        unsigned int dstId = chunkBasePos + chunkOffset;
        if (dstId < srcIndices.size())
        {
            unsigned int srcId = srcIndices.template get<0>(dstId);
            dst.get(dstId) = src.get(srcId); //todo How does this = spread across threads...?
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }

    template<typename DataVectorT, typename IndexVectorT>
    __global__ void mergeData(DataVectorT data1, DataVectorT data2, IndexVectorT mergeIndices, DataVectorT dataOut)
    {
#ifdef __NVCC__
        typedef typename DataVectorT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, 0> BlockT0; // The type of the 0-th property
        unsigned int chunkSize = BlockT0::size;

        unsigned int chunksPerBlock = blockDim.x / chunkSize;
        unsigned int chunkOffset = threadIdx.x / chunkSize; // The thread block can work on several chunks in parallel

        unsigned int chunkBasePos = blockIdx.x * chunksPerBlock;

        unsigned int dstId = chunkBasePos + chunkOffset;
        if (dstId < mergeIndices.size())
        {
            auto size1 = data1.size();
            auto size2 = data2.size();
            unsigned int mrgId = mergeIndices.template get<0>(dstId);
            unsigned int srcId = mrgId - size1; // This is actually only used on data2 branch!
            if (mrgId < size1) // We need to read from data1
            {
                dataOut.get(dstId) = data1.get(mrgId); //todo How does this = spread across threads...?
            }
            else if (srcId < size2) // We need to read from data2
            {
                dataOut.get(dstId) = data2.get(srcId); //todo How does this = spread across threads...?
            }
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }

    // Below the kernels to be used inside the "compute segments" part of the solveConflicts
    /**
     * Fill a vector of predicates specifying if current key is different from the previous one.
     *
     * @tparam IndexVectorT The type of the (OpenFPM) vector of indices.
     * @param keys The keys to compute the predicate on.
     * @param predicates The output vector of {0,1} predicates.
     */
    template<typename IndexVectorT>
    __global__ void computePredicates(IndexVectorT keys, IndexVectorT predicates)
    {
#ifdef __NVCC__
       // predicates[i] must be 1 if keys[i] != keys[i-1] or if i == 0
        // else it must be 0
        unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
        if (pos < keys.size())
        {
            if (pos > 0)
            {
                predicates.template get<0>(pos) = (keys.template get<0>(pos - 1) != keys.template get<0>(pos)) ? 1 : 0;
            }
            else
            {
                predicates.template get<0>(pos) = 1;
            }
        }
        else if (pos == keys.size())
        {
            // This is just to get one extra trailing element in the exScan
            // (value here doesn't matter, as the scan is exclusive!)
            predicates.template get<0>(pos) = 1;
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }

    /**
     * Copy the thread id to the destination position in the out array specified by dstIndices at the corresponding position.
     * 
     * @tparam IndexVectorT 
     * @param keysSize 
     * @param dstIndices 
     * @param out 
     */
    template<typename IndexVectorT>
    __global__ void copyIdToDstIndexIfPredicate(size_t keysSize, IndexVectorT dstIndices, IndexVectorT out)
    {
#ifdef __NVCC__
        // dstIndices is exclusive scan of predicates
        unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
        if (pos < keysSize) // dstIndices.size() must be keysSize+1, so a trailing element is required
        {
            bool predicate = dstIndices.template get<0>(pos) != dstIndices.template get<0>(pos+1);
            if (predicate)
            {
                auto dstPos = dstIndices.template get<0>(pos);
                out.template get<0>(dstPos) = pos;
            }
        }
        else if (pos == keysSize)
        {
            // Append a trailing index as end of last segment
            auto lastDst = dstIndices.template get<0>(pos - 1);
            auto lastPredicate = dstIndices.template get<0>(pos - 1) != dstIndices.template get<0>(pos);
            out.template get<0>(lastDst +
                                lastPredicate) = pos; // We need to increment position if the last element was also written
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }


    /**
     * Copy the key to the destination position in the out array specified by dstIndices.
     * NOTE: it is not possible to reorder in-place!
     * 
     * @tparam IndexVectorT 
     * @param keys
     * @param dstIndices 
     * @param out 
     */
    template<typename IndexVectorT>
    __global__ void copyKeyToDstIndexIfPredicate(IndexVectorT keys, IndexVectorT dstIndices, IndexVectorT out)
    {
#ifdef __NVCC__
       // dstIndices is exclusive scan of predicates
        unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
        if (pos < keys.size()) // dstIndices.size() must be keysSize+1, so a trailing element is required
        {
            bool predicate = dstIndices.template get<0>(pos) != dstIndices.template get<0>(pos+1);
            if (predicate)
            {
                auto dstPos = dstIndices.template get<0>(pos);
                out.template get<0>(dstPos) = keys.template get<0>(pos);
            }
        }
    }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
}

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<unsigned int blockSize,
		typename vector_data_type,
		typename vector_index_type,
        typename vector_index_type2,
        typename vector_reduction,
        typename block_functor,
        unsigned int impl2, unsigned int pSegment=1>
struct sparse_vector_reduction_solve_conflict
{
	//! Vector in which to the reduction
	vector_data_type & vector_data_red;

	//! new datas
	vector_data_type & vector_data;

	//! new data in an unsorted way
	vector_data_type & vector_data_unsorted;

	//! segment of offsets
	vector_index_type2 & segment_offset;

	//! map of the data
	vector_index_type & vector_data_map;

	//! gpu context
	mgpu::ofp_context_t & context;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline sparse_vector_reduction_solve_conflict(vector_data_type & vector_data_red,
								   vector_data_type & vector_data,
								   vector_data_type & vector_data_unsorted,
								   vector_index_type & vector_data_map,
								   vector_index_type2 & segment_offset,
								   mgpu::ofp_context_t & context)
	:vector_data_red(vector_data_red),vector_data(vector_data),vector_data_unsorted(vector_data_unsorted),segment_offset(segment_offset),vector_data_map(vector_data_map),context(context)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
#ifdef __NVCC__

        typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
        typedef typename boost::mpl::at<typename ValueTypeOf<vector_data_type>::type,typename reduction_type::prop>::type red_type;
        if (reduction_type::is_special() == false)
		{
            typedef typename std::remove_reference<vector_data_type>::type::value_type AggregateT;

            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<
                    typename std::remove_reference<vector_data_type>::type::value_type::type,
                    typename reduction_type::prop
                    >::type::scalarType red_type;
            typedef typename reduction_type::op_red<red_type> red_op;

            constexpr unsigned int p = reduction_type::prop::value;
            constexpr unsigned int pMask = AggregateT::max_prop_real - 1;

            typedef BlockTypeOf<AggregateT, p> BlockT; // The type of the 0-th property

            constexpr unsigned int chunksPerBlock = blockSize / BlockT::size;
            const unsigned int gridSize =
                    segment_offset.size() - 1; // This "-1" is because segments has a trailing extra element

            CUDA_LAUNCH_DIM3((BlockMapGpuKernels::segreduce<p, pSegment, pMask, chunksPerBlock, red_op>),gridSize, blockSize,
                    vector_data.toKernel(),
                    segment_offset.toKernel()/*,
                    src_map.toKernel()*/,
                    vector_data.toKernel(),
                    vector_data_red.toKernel());
		}
#else
		std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif
	}
};

namespace BlockMapGpuFunctors
{
    /**
     * This functor is used in the sparse vector flush method to achieve the right blocked behaviour
     */
    template<unsigned int blockSize>
    struct BlockFunctor
    {
        template<unsigned int pSegment, typename vector_reduction, typename T, typename vector_index_type, typename vector_index_type2, typename vector_data_type>
        static bool seg_reduce(vector_index_type2 &segments, vector_data_type &src, vector_data_type &src_unsorted, vector_index_type & src_map, vector_data_type &dst)
        {
#ifdef __NVCC__
            typedef typename vector_data_type::value_type AggregateT;

            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<
                    typename vector_data_type::value_type::type,
                    typename reduction_type::prop
                    >::type::scalarType red_type;
            typedef typename reduction_type::op_red<red_type> red_op;

            constexpr unsigned int p = reduction_type::prop::value;
            constexpr unsigned int pMask = AggregateT::max_prop_real - 1;

            typedef BlockTypeOf<AggregateT, p> BlockT; // The type of the 0-th property

            constexpr unsigned int chunksPerBlock = blockSize / BlockT::size;
            const unsigned int gridSize =
                    segments.size() - 1; // This "-1" is because segments has a trailing extra element

            CUDA_LAUNCH_DIM3((BlockMapGpuKernels::segreduce_beta<p, pSegment, pMask, chunksPerBlock, red_op>),gridSize, blockSize,
                    src_unsorted.toKernel(),
                    segments.toKernel(),
                    src_map.toKernel(),
                    src.toKernel(),
                    dst.toKernel());

            return true;
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
            return true;
#endif // __NVCC__
        }

        template<typename vector_index_type, typename vector_data_type, typename ... v_reduce>
        static bool solve_conflicts(vector_index_type &keys, vector_index_type &mergeIndices,
                                    vector_data_type &dataOld, vector_data_type &dataNew,
                                    vector_index_type &tmpIndices, vector_data_type &tmpData,
                                    vector_index_type &keysOut, vector_data_type &dataOut,
                                    mgpu::ofp_context_t & context)
        {
#ifdef __NVCC__
            typedef ValueTypeOf<vector_data_type> AggregateT;
            typedef ValueTypeOf<vector_index_type> AggregateIndexT;

            typedef BlockTypeOf<AggregateIndexT, 0> IndexT;

            typedef BlockTypeOf<AggregateT, 0> BlockT0; // The type of the 0-th property

            constexpr unsigned int chunksPerBlock = blockSize / BlockT0::size;
            const unsigned int gridSize = mergeIndices.size() / chunksPerBlock + ((mergeIndices.size() % chunksPerBlock) != 0);

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Calculate maps

/*        	CUDA_LAUNCH(compute_predicate,ite,mergeIndices.toKernel(),dataOld.size(),p_ids.toKernel());

        	mgpu::standard_context_t context(false);
        	mgpu::scan((int *)p_ids.template getDeviceBuffer<0>(),
        				s_ids.size(),
        	            (int *)s_ids.template getDeviceBuffer<0>(),
                        context);

        	mgpu::scan((int *)p_ids.template getDeviceBuffer<1>(),
        				s_ids.size(),
        	            (int *)s_ids.template getDeviceBuffer<1>(),
                        context);

        	mgpu::scan((int *)p_ids.template getDeviceBuffer<2>(),
        				s_ids.size(),
        	            (int *)s_ids.template getDeviceBuffer<2>(),
                        context);

        	mgpu::scan((int *)p_ids.template getDeviceBuffer<3>(),
        				s_ids.size(),
        	            (int *)s_ids.template getDeviceBuffer<3>(),
                        context);

        	s_ids.template deviceToHost<0,1,2,3>();
        	p_ids.template deviceToHost<0,1,2,3>();

        	size_t copy_old_size = s_ids.template get<3>(s_ids.size()-1) + p_ids.template get<3>(p_ids.size()-1);
        	size_t seg_old_size = s_ids.template get<1>(s_ids.size()-1) + p_ids.template get<1>(p_ids.size()-1);
        	size_t out_map_size = s_ids.template get<1>(s_ids.size()-1) + p_ids.template get<1>(p_ids.size()-1);

        	segments_oldData.resize(seg_old_size);
        	outputMap.resize(out_map_size);
        	copy_old.resize(copy_old_size);

        	CUDA_LAUNCH(maps_create,ite,s_ids.toKernel(),p_ids.toKernel(),segments_oldData.toKernel(),outputMap.toKernel(),copy_old.toKernel());*/



            /////////////////////////////////////////////////////////////////////////////////////////////////////////////

            //todo Work plan
            // Phases of solve_conflicts:
            // 0) merge data
            // 1) compute segments
            //      a) compute predicates
            //      b) run exclusive scan on predicates
            //      c) copyIdToDst...
            //      d) update the keys to the soon-to-be segreduced ones
            // 2) perform segreduce

            // First ensure we have the right sizes on the buffers

            tmpData.resize(mergeIndices.size(),EXACT_RESIZE); //todo: check if we need some other action to actually have the right space on gpu

            // Phase 0 - merge data
            BlockMapGpuKernels::mergeData<<< gridSize, blockSize >>>(dataOld.toKernel(),
                    dataNew.toKernel(), mergeIndices.toKernel(), tmpData.toKernel());
            // Now we can recycle the mergeIndices buffer for other things.

            // Phase 1 - compute segments
            mergeIndices.resize(mergeIndices.size() + 1); // This is to get space for the extra trailing predicate
            tmpIndices.resize(mergeIndices.size());

            unsigned int gridSize2 = static_cast<unsigned int>(std::ceil(static_cast<float>(mergeIndices.size()) / blockSize));
            BlockMapGpuKernels::computePredicates << <
                        gridSize2,
                        blockSize>>>
                        (keys.toKernel(), mergeIndices.toKernel());

            mgpu::scan( (IndexT*) mergeIndices.template getDeviceBuffer<0>(),
                        mergeIndices.size(), // Here note that mergeIndices is the predicate vector
                        (IndexT*) tmpIndices.template getDeviceBuffer<0>(),
                        context); // mgpu scan is exclusive by default
            // Now it's important to resize mergeIndices in the right way, otherwise dragons ahead!
            tmpIndices.template deviceToHost<0>(tmpIndices.size()-1, tmpIndices.size()-1); //getting only the last element from device
            mergeIndices.resize(tmpIndices.template get<0>(tmpIndices.size()-1) + 1);
            BlockMapGpuKernels::copyIdToDstIndexIfPredicate<< < gridSize2, blockSize >> >
                                                                            (keys.size(), tmpIndices.toKernel(), mergeIndices.toKernel());
            // mergeIndices now contains the segments

            // Now update the keys
            keysOut.resize(mergeIndices.size()-1); // The final number of keys is one less than the segments values
            BlockMapGpuKernels::copyKeyToDstIndexIfPredicate<< < gridSize2, blockSize >> >
                                                                             (keys.toKernel(), tmpIndices.toKernel(), keysOut.toKernel());
            // the new keys are now in keysOut

            // Phase 2 - segreduce on all properties
            dataOut.resize(mergeIndices.size()-1); // Right size for output, i.e. the number of segments
            typedef boost::mpl::vector<v_reduce...> vv_reduce;
            constexpr unsigned int pSegment = 0;
            sparse_vector_reduction_solve_conflict<blockSize,decltype(dataOut),decltype(mergeIndices),decltype(mergeIndices),vv_reduce,BlockFunctor,2, pSegment>
                    svr(dataOut,tmpData,tmpData,mergeIndices,mergeIndices,context);
            boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(svr);


            return true; //todo: check if error in kernel
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
            return true;
#endif // __NVCC__
        }
    };
}



namespace BlockMapGpuFunctorsNew
{
    /**
     * This functor is used in the sparse vector flush method to achieve the right blocked behaviour
     */
    template<unsigned int blockSize>
    struct BlockFunctor
    {
        template<typename vector_index_type, typename vector_data_type>
        static bool reorder(vector_index_type &src_id, vector_data_type &data, vector_data_type &data_reord)
        {
            return true;
        }

        template<unsigned int pSegment, typename vector_reduction, typename T, typename vector_index_type, typename vector_data_type>
        static bool seg_reduce(vector_index_type &segments, vector_data_type &data_old, vector_data_type &data_output,
        					   vector_index_type &vct_old_copy,
        					   vector_data_type &data_new, vector_index_type data_new_map,
        					   vector_index_type & data_out_map,
        					   vector_index_type & data_seg_old_map)
        {
#ifdef __NVCC__
            typedef typename vector_data_type::value_type AggregateT;

            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<
                    typename vector_data_type::value_type::type,
                    typename reduction_type::prop
                    >::type::scalarType red_type;
            typedef typename reduction_type::op_red<red_type> red_op;

            constexpr unsigned int p = reduction_type::prop::value;
            constexpr unsigned int pMask = AggregateT::max_prop_real - 1;

            typedef BlockTypeOf<AggregateT, p> BlockT; // The type of the 0-th property

            constexpr unsigned int chunksPerBlock = blockSize / BlockT::size;
            const unsigned int gridSize =
                    segments.size() - 1; // This "-1" is because segments has a trailing extra element

			// Copy old elements

			CUDA_LAUNCH_DIM3(copy_old,gridSize,data_output.toKernel(),vct_old_copy.toKernel(),data_old.toKernel());


    		CUDA_LAUNCH_DIM3((BlockMapGpuKernels::segreduce_total<p, pSegment, pMask, 2, chunksPerBlock, red_op>),gridSize, blockSize,
    						data_new.toKernel(),
    						data_old.toKernel(),
    						segments.toKernel(),
    						data_new_map.toKernel(),
    						data_seg_old_map.toKernel(),
    						data_out_map.toKernel(),
    						data_output.toKernel());

            return true; //todo: check if error in kernel
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
            return true;
#endif // __NVCC__
        }

        template<unsigned int pSegment, typename vector_reduction, typename T,
                 typename vector_index_type,
                 typename vector_index_type2,
                 typename vector_data_type>
        static bool seg_reduce_total(vector_index_type2 &segments, vector_data_type &data_old, vector_data_type &data_output,
        					   vector_index_type &vct_old_copy,
        					   vector_data_type &data_new, vector_index_type data_new_map,
        					   vector_index_type & data_out_map,
        					   vector_index_type & data_seg_old_map)
        {
#ifdef __NVCC__
            typedef typename vector_data_type::value_type AggregateT;

            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<
                    typename vector_data_type::value_type::type,
                    typename reduction_type::prop
                    >::type::scalarType red_type;
            typedef typename reduction_type::op_red<red_type> red_op;

            constexpr unsigned int p = reduction_type::prop::value;
            constexpr unsigned int pMask = AggregateT::max_prop_real - 1;

            typedef BlockTypeOf<AggregateT, p> BlockT; // The type of the 0-th property

            constexpr unsigned int chunksPerBlock = blockSize / BlockT::size;
            const unsigned int gridSize =
                    segments.size() - 1; // This "-1" is because segments has a trailing extra element

			// Copy old elements

			CUDA_LAUNCH_DIM3(copy_old,gridSize,blockSize,data_output.toKernel(),vct_old_copy.toKernel(),data_old.toKernel());


    		CUDA_LAUNCH_DIM3((BlockMapGpuKernels::segreduce_total<p, pSegment, pMask, chunksPerBlock, red_op>),gridSize, blockSize,
    						data_new.toKernel(),
    						data_old.toKernel(),
    						segments.toKernel(),
    						data_new_map.toKernel(),
    						data_seg_old_map.toKernel(),
    						data_out_map.toKernel(),
    						data_output.toKernel());

            return true; //todo: check if error in kernel
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
            return true;
#endif // __NVCC__
        }

        template<unsigned int pSegment, typename vector_reduction, typename T, typename vector_index_type, typename vector_data_type>
        static bool seg_reduce(vector_index_type &segments, vector_data_type &src, vector_data_type &dst)
        {
#ifdef __NVCC__
            typedef typename vector_data_type::value_type AggregateT;

            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<
                    typename vector_data_type::value_type::type,
                    typename reduction_type::prop
                    >::type::scalarType red_type;
            typedef typename reduction_type::op_red<red_type> red_op;

            constexpr unsigned int p = reduction_type::prop::value;
            constexpr unsigned int pMask = AggregateT::max_prop_real - 1;

            typedef BlockTypeOf<AggregateT, p> BlockT; // The type of the 0-th property

            constexpr unsigned int chunksPerBlock = blockSize / BlockT::size;
            const unsigned int gridSize =
                    segments.size() - 1; // This "-1" is because segments has a trailing extra element

            BlockMapGpuKernels::segreduce<p, pSegment, pMask, chunksPerBlock, red_op> <<< gridSize, blockSize >>> (
                    src.toKernel(),
                    segments.toKernel(),
                    src.toKernel(),
                    dst.toKernel()
            );
            return true; //todo: check if error in kernel
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
            return true;
#endif // __NVCC__
        }

        template<typename vector_index_type, typename vector_data_type, typename ... v_reduce>
        static bool solve_conflicts(vector_index_type &keys, vector_index_type &mergeIndices,
                                    vector_data_type &dataOld, vector_data_type &dataNew,
                                    vector_index_type &tmpIndices, vector_data_type &tmpData,
                                    vector_index_type &keysOut, vector_data_type &dataOut,
                                    mgpu::ofp_context_t & context)
        {
            return true;
        }
    };
}

//#endif //__NVCC__

#endif //OPENFPM_PDATA_BLOCKMAPGPU_KERNELS_CUH
