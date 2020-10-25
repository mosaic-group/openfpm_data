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
#include "util/cuda/scan_ofp.cuh"

namespace BlockMapGpuKernels
{
	// for each segments
	template<typename vector_index_type, typename vector_index_type2>
	__global__ void compute_predicate(vector_index_type vct_keys_merge,
									vector_index_type vct_index_merge,
									unsigned int m,
									vector_index_type2 pred_ids)
	{
		int p = blockIdx.x * blockDim.x + threadIdx.x;

		if (p >= vct_index_merge.size())	return;

		unsigned int pp1 = (p+1 == vct_index_merge.size())?p:p+1;
		unsigned int pm1 = (p == 0)?p:p-1;

		auto id0 = vct_index_merge.template get<0>(pm1);
		auto id1 = vct_index_merge.template get<0>(p);

		auto k0 = vct_keys_merge.template get<0>(pm1);
		auto k1 = vct_keys_merge.template get<0>(p);
		auto k2 = vct_keys_merge.template get<0>(pp1);

		// predicate 0 count old chunks, but when is merged to new data the 1 must be shifted to the new element
		//
		pred_ids.template get<0>(p) = ((k0 == k1) && (p != 0) && id0 < m) || (id1 < m && (k1 != k2));

		//predicate 1 is used to count the new index segments
		pred_ids.template get<1>(p) = id1 >= m;

		// predicate 2 is used is used to count everything does not merge
		pred_ids.template get<2>(p) = (k1 != k2) | (p == vct_index_merge.size()-1);

		// predicate 3 is used to count old keys that does not reduce with new data
		pred_ids.template get<3>(p) = id1 < m && ((k1 != k2) | (p == vct_index_merge.size()-1));

		//predicate 1 is used to count the new index segments
		pred_ids.template get<4>(p) = id1 < m;
	}

	// for each segments
	template<typename vector_index_type, typename vector_index_type2, typename vector_index_map_type>
	__global__ void maps_create(vector_index_type2 scan_ids,
									vector_index_type2 pred_ids,
									vector_index_type vct_seg_old,
									vector_index_map_type vct_out_map,
									vector_index_type vct_copy_old_dst,
									vector_index_type vct_copy_old_src)
	{
		int p = blockIdx.x * blockDim.x + threadIdx.x;

		if (p >= scan_ids.size())	return;

		auto id1 = scan_ids.template get<0>(p);
		bool pred_id1 = pred_ids.template get<0>(p);
		auto id2 = scan_ids.template get<1>(p);
		bool pred_id2 = pred_ids.template get<1>(p);
		auto id3 = scan_ids.template get<2>(p);
		auto id4 = scan_ids.template get<3>(p);
		bool pred_id4 = pred_ids.template get<3>(p);
		auto id5 = scan_ids.template get<4>(p);

		if (pred_id2 == true)
		{vct_seg_old.template get<0>(id2) = (pred_id1 == true)?id1:-1;}

		if (pred_id2 == true)
		{vct_out_map.template get<0>(id2) = id3;}

		if (pred_id4 == true)
		{
			vct_copy_old_dst.template get<0>(id4) = id3;
			vct_copy_old_src.template get<0>(id4) = id5;
		}
	}



	// for each segments
	template<typename vector_index_type, typename vector_data_type>
	__global__ void copy_old(vector_data_type vct_data_new,
									vector_index_type index_to_copy,
									vector_data_type vct_data_old)
	{
		int p = blockIdx.x * blockDim.x + threadIdx.x;

		if (p >= index_to_copy.size())	return;

		auto id = index_to_copy.template get<0>(p);

		vct_data_new.get(id) = vct_data_old.get(p);
	}

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
            typename IndexVector_segdataT, typename IndexVector_datamapT,
            typename IndexVector_segoldT, typename IndexVector_outmapT, typename DataVectorT>
    __global__ void
    segreduce_total(
            DataVectorT data_new,
            DataVectorT data_old,
            IndexVector_segdataT segments_data,
            IndexVector_datamapT segments_dataMap,
            IndexVector_segoldT segments_dataOld,
            IndexVector_outmapT outputMap,
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

        if (seg_old != -1 && chunkId == 0)
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

    // GridSize = number of segments
    // BlockSize = chunksPerBlock * chunkSize
    //
    template<unsigned int p,
            unsigned int pSegment,
            unsigned int pMask,
            unsigned int chunksPerBlock,
            typename op,
            typename IndexVector_segdataT, typename IndexVector_datamapT,
            typename IndexVector_segoldT, typename IndexVector_outmapT, typename DataVectorT>
    __global__ void
    segreduce_total_with_mask(
            DataVectorT data_new,
            DataVectorT data_old,
            IndexVector_segdataT segments_data,
            IndexVector_datamapT segments_dataMap,
            IndexVector_segoldT segments_dataOld,
            IndexVector_outmapT outputMap,
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

        if (seg_old != -1 && chunkId == 0)
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

            generalDimensionFunctor<MaskType>::assignWithOffset(output.template get<pMask>(out_id), AMask[chunkId],
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
    __global__ void copy_old_ker(IndexVectorT srcIndices, DataVectorT src, IndexVectorT dstIndices, DataVectorT dst)
    {
#ifdef __NVCC__
        typedef typename DataVectorT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, 0> BlockT0; // The type of the 0-th property
        unsigned int chunkSize = BlockT0::size;

        unsigned int chunksPerBlock = blockDim.x / chunkSize;
        unsigned int chunkOffset = threadIdx.x / chunkSize; // The thread block can work on several chunks in parallel

        unsigned int chunkBasePos = blockIdx.x * chunksPerBlock;

        unsigned int p = chunkBasePos + chunkOffset;
        if (p < srcIndices.size())
        {
            auto dstId = dstIndices.template get<0>(p);
            auto srcId = srcIndices.template get<0>(p);

            dst.get(dstId) = src.get(srcId);
        }
#else // __NVCC__
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif // __NVCC__
    }


    template<typename IndexVectorT, typename IndexVectorT2>
    __global__ void copyKeyToDstIndexIfPredicate(IndexVectorT keys, IndexVectorT2 dstIndices, IndexVectorT out)
    {
#ifdef __NVCC__
       // dstIndices is exclusive scan of predicates
        unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

        if (pos >= dstIndices.size())	{return;}

        unsigned int pm1 = (pos == 0)?0:pos-1;

        bool predicate = dstIndices.template get<2>(pos) != dstIndices.template get<2>(pm1) || pos == 0;

		if (predicate)
		{
			auto dstPos = dstIndices.template get<2>(pos);
			out.template get<0>(dstPos) = keys.template get<0>(pos);
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
		typename vector_datamap_type,
        typename vector_segoffset_type,
        typename vector_outmap_type,
        typename vector_segolddata_type,
        typename vector_reduction,
        typename block_functor,
        unsigned int impl2, unsigned int pSegment=1>
struct sparse_vector_reduction_solve_conflict
{
	//! Vector in which to the reduction
	vector_data_type & vector_data_red;

	//! new datas
	vector_data_type & vector_data;

	//! new datas
	vector_data_type & vector_data_old;

	//! new data in an unsorted way
	vector_data_type & vector_data_unsorted;

	//! segment of offsets
	vector_segoffset_type & segment_offset;

	//! map of the data
	vector_datamap_type & vector_data_map;

	//! output map
	vector_outmap_type & out_map;

	//! old data segments
	vector_segolddata_type & segments_oldData;

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
								   vector_data_type & vector_data_old,
								   vector_datamap_type & vector_data_map,
								   vector_segoffset_type & segment_offset,
								   vector_outmap_type & out_map,
								   vector_segolddata_type & segments_oldData,
								   mgpu::ofp_context_t & context)
	:vector_data_red(vector_data_red),
	 vector_data(vector_data),
	 vector_data_unsorted(vector_data_unsorted),
	 vector_data_old(vector_data_old),
	 segment_offset(segment_offset),
	 vector_data_map(vector_data_map),
	 out_map(out_map),
	 segments_oldData(segments_oldData),
	 context(context)
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
            typedef typename std::remove_all_extents<typename boost::mpl::at<
                    typename std::remove_reference<vector_data_type>::type::value_type::type,
                    typename reduction_type::prop
                    >::type>::type::scalarType red_type;
            typedef typename reduction_type::template op_red<red_type> red_op;

            constexpr unsigned int p = reduction_type::prop::value;
            constexpr unsigned int pMask = AggregateT::max_prop_real - 1;

            typedef typename std::remove_all_extents<BlockTypeOf<AggregateT, p>>::type BlockT; // The type of the 0-th property

            constexpr unsigned int chunksPerBlock = blockSize / BlockT::size;
            const unsigned int gridSize =
                    segment_offset.size()  - 1; // This "-1" is because segments has a trailing extra element

            if (T::value == 0)
            {
            	CUDA_LAUNCH_DIM3((BlockMapGpuKernels::segreduce_total_with_mask<p, pSegment, pMask, chunksPerBlock, red_op>),gridSize, blockSize,
                    vector_data.toKernel(),
            		  vector_data_old.toKernel(),
                    segment_offset.toKernel(),
                    vector_data_map.toKernel(),
                    segments_oldData.toKernel(),
                    out_map.toKernel(),
                    vector_data_red.toKernel());
            }
            else
            {
            	CUDA_LAUNCH_DIM3((BlockMapGpuKernels::segreduce_total<p, pSegment, pMask, chunksPerBlock, red_op>),gridSize, blockSize,
                    vector_data.toKernel(),
            		  vector_data_old.toKernel(),
                    segment_offset.toKernel(),
                    vector_data_map.toKernel(),
                    segments_oldData.toKernel(),
                    out_map.toKernel(),
                    vector_data_red.toKernel());
            }
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
    	openfpm::vector_gpu<aggregate<int,int,int,int,int>> p_ids;
    	openfpm::vector_gpu<aggregate<int,int,int,int,int>> s_ids;

    	openfpm::vector_gpu<aggregate<int>> copy_old_dst;
    	openfpm::vector_gpu<aggregate<int>> copy_old_src;
    	openfpm::vector_gpu<aggregate<unsigned int>> outputMap;
    	openfpm::vector_gpu<aggregate<int>> segments_oldData;

    	/*! \brief Create the array of the merged datas with no more conflicts (repeated chunk ids)
    	 *
    	 * \param keys Old chunk id already present in the data structure
    	 * \param mergeIndeces array of the merged indexes (with conflicts repeated chunks ids)
    	 * \param segments_new
    	 * \param dataMap it store the index where the data is located in the original add buffer
    	 * \param dataOld chunk data for the old chunks
    	 * \param dataNew chunk data for the new chunk
    	 * \param keysOut array with the chunks ids (without conflicts)
    	 * \param dataOut array of the merged data (without conflicts)
    	 *
    	 */
        template<unsigned int pSegment, typename vector_index_type, typename vector_index_type2, typename vector_data_type, typename ... v_reduce>
        bool solve_conflicts(vector_index_type &keys, vector_index_type &mergeIndices, vector_index_type2 &segments_new, vector_index_type &data_map,
                                    vector_data_type &dataOld, vector_data_type &dataNew,
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

            auto ite = mergeIndices.getGPUIterator();

            p_ids.resize(mergeIndices.size());
            s_ids.resize(mergeIndices.size());

            // shut-up valgrind uninitialized

            p_ids.template get<1>(p_ids.size()-1) = 0;
        	CUDA_LAUNCH(BlockMapGpuKernels::compute_predicate,ite,keys.toKernel(),mergeIndices.toKernel(),dataOld.size(),p_ids.toKernel());

        	openfpm::scan((int *)p_ids.template getDeviceBuffer<0>(),
        				s_ids.size(),
        	            (int *)s_ids.template getDeviceBuffer<0>(),
                        context);

        	openfpm::scan((int *)p_ids.template getDeviceBuffer<1>(),
        				s_ids.size(),
        	            (int *)s_ids.template getDeviceBuffer<1>(),
                        context);

        	openfpm::scan((int *)p_ids.template getDeviceBuffer<2>(),
        				s_ids.size(),
        	            (int *)s_ids.template getDeviceBuffer<2>(),
                        context);

        	openfpm::scan((int *)p_ids.template getDeviceBuffer<3>(),
        				s_ids.size(),
        	            (int *)s_ids.template getDeviceBuffer<3>(),
                        context);

        	openfpm::scan((int *)p_ids.template getDeviceBuffer<4>(),
        				s_ids.size(),
        	            (int *)s_ids.template getDeviceBuffer<4>(),
                        context);

        	s_ids.template deviceToHost<0,1,2,3>(s_ids.size()-1,s_ids.size()-1);
        	p_ids.template deviceToHost<0,1,2,3>(p_ids.size()-1,p_ids.size()-1);

        	size_t copy_old_size = s_ids.template get<3>(s_ids.size()-1) + p_ids.template get<3>(p_ids.size()-1);
        	size_t seg_old_size = s_ids.template get<1>(s_ids.size()-1) + p_ids.template get<1>(p_ids.size()-1);
        	size_t out_map_size = s_ids.template get<1>(s_ids.size()-1) + p_ids.template get<1>(p_ids.size()-1);
        	size_t data_out_size = s_ids.template get<2>(s_ids.size()-1) + p_ids.template get<2>(p_ids.size()-1);

        	segments_oldData.resize(seg_old_size);
        	outputMap.resize(out_map_size);
        	copy_old_dst.resize(copy_old_size);
        	copy_old_src.resize(copy_old_size);

        	CUDA_LAUNCH(BlockMapGpuKernels::maps_create,ite,s_ids.toKernel(),p_ids.toKernel(),segments_oldData.toKernel(),outputMap.toKernel(),copy_old_dst.toKernel(),copy_old_src.toKernel());

            // Create the output for the keys
            keysOut.resize(data_out_size); // The final number of keys is one less than the segments values

            ite = keys.getGPUIterator();
            CUDA_LAUNCH(BlockMapGpuKernels::copyKeyToDstIndexIfPredicate,ite,keys.toKernel(), s_ids.toKernel(), keysOut.toKernel());



            // the new keys are now in keysOut

            // Phase 2 - segreduce on all properties
            dataOut.reserve(data_out_size+1);
            dataOut.resize(data_out_size); // Right size for output, i.e. the number of segments
            typedef boost::mpl::vector<v_reduce...> vv_reduce;

			sparse_vector_reduction_solve_conflict<blockSize,decltype(dataOut),
															 decltype(data_map),
															 decltype(segments_new),
															 decltype(outputMap),
															 decltype(segments_oldData),
															 vv_reduce,BlockFunctor,2, pSegment>
					svr(dataOut,dataNew,dataNew,dataOld,data_map,segments_new,outputMap,segments_oldData,context);

            boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(svr);

            //copy the old chunks
            if (copy_old_dst.size() != 0)
            {
            	CUDA_LAUNCH_DIM3(BlockMapGpuKernels::copy_old_ker,copy_old_dst.size(),blockSize,copy_old_src.toKernel(),dataOld.toKernel(),copy_old_dst.toKernel(),dataOut.toKernel());
            }

            return true; //todo: check if error in kernel
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
            return true;
#endif // __NVCC__
        }

        openfpm::vector_gpu<aggregate<unsigned int>> & get_outputMap()
		{
        	return outputMap;
		}

        const openfpm::vector_gpu<aggregate<unsigned int>> & get_outputMap() const
		{
        	return outputMap;
		}
    };
}


//#endif //__NVCC__

#endif //OPENFPM_PDATA_BLOCKMAPGPU_KERNELS_CUH
