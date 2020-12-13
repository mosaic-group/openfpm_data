/*
 * SparseGridGpu_ker_util.hpp
 *
 *  Created on: Aug 7, 2019
 *      Author: i-bird
 */

#ifndef SPARSEGRIDGPU_KER_UTIL_HPP_
#define SPARSEGRIDGPU_KER_UTIL_HPP_

#include "util/variadic_to_vmpl.hpp"

template<bool to_set>
struct set_compile_condition
{
	template<unsigned int p, typename SrcType, typename AggrType>
	__device__ __host__ static inline void set(SrcType & src,AggrType & aggr)
	{
		src = aggr.template get<p>();
	}
};

template<>
struct set_compile_condition<false>
{
	template<unsigned int p, typename SrcType, typename AggrType>
	__device__ __host__ static inline void set(SrcType & src,AggrType & aggr)
	{}
};

template<unsigned int dim,typename T>
struct cross_stencil
{
	T xm[dim];
	T xp[dim];
};

/*template<unsigned int dim, unsigned int block_edge_size>
struct shift_position
{
    __device__ static inline int shift(int pos, int stencilRadius)
    {
        int accu = 1;
        int pos_s = 0;
        for (int i = 0 ; i < dim ; i++)
        {
            pos_s += (pos % block_edge_size + stencilRadius)*accu;
            accu *= (block_edge_size + 2*stencilRadius);
            pos /= block_edge_size;
        }

        return pos_s;
    }
};

template<unsigned int block_edge_size>
struct shift_position<2,block_edge_size>
{
    __device__ static inline int shift(int pos, int stencilRadius)
    {
        unsigned int x = pos % block_edge_size;
        unsigned int y = (pos / block_edge_size);

        unsigned int g_sz = block_edge_size + 2*stencilRadius;

        return (x+stencilRadius) + (y+stencilRadius)*g_sz;
    }
};


template<unsigned int block_edge_size>
struct shift_position<3,block_edge_size>
{
    __device__ static inline int shift(int pos, int stencilRadius)
    {
        unsigned int x = pos % block_edge_size;
        unsigned int y = (pos / block_edge_size) % block_edge_size;
        unsigned int z = (pos / (block_edge_size*block_edge_size));

        unsigned int g_sz = block_edge_size + 2*stencilRadius;

        return (x+stencilRadius) + (y+stencilRadius)*g_sz + (z+stencilRadius)*g_sz*g_sz;
    }
};*/

template<unsigned int dim>
struct NNStar
{
	static const int nNN = IntPow<2, dim>::value;

	template<typename indexT, typename blockCoord_type, typename blockMap_type, typename SparseGrid_type>
	__device__ static inline indexT getNNpos(blockCoord_type & blockCoord,
								  blockMap_type & blockMap,
								  SparseGrid_type & sparseGrid,
								  const unsigned int offset)
	{
        //todo: also do the full neighbourhood version, this is just cross
        int neighbourPos = -1;
        if (offset < 2*dim)
        {
            unsigned int d = offset/2;
            int dPos = blockCoord.get(d) + (offset%2)*2 - 1;
            blockCoord.set_d(d, dPos);

            int bl = sparseGrid.getBlockLinId(blockCoord);

            bl = (dPos < 0)?-1:bl;

            neighbourPos = blockMap.get_sparse(bl).id;
        }
        return neighbourPos;
	}

	template<typename indexT, unsigned int blockEdgeSize, typename coordType>
	__host__ static inline indexT getNNskin(coordType & coord, int stencilSupportRadius)
	{
        int neighbourNum = -1;
        int ctr = 0;
        for (int j = 0; j < dim; ++j)
        {
            int c = static_cast<int>(coord.get(j)) - static_cast<int>(stencilSupportRadius);
            if (c < 0)
            {
                neighbourNum = 2*j;
                ++ctr;
            }
            else if (c >= blockEdgeSize)
            {
                neighbourNum = 2*j + 1;
                ++ctr;
            }
        }
        if (ctr > 1) // If we are on a "corner"
        {
            neighbourNum = 0;
        }

        return neighbourNum;
	}

	template<typename sparseGrid_type, typename coord_type, typename Mask_type,unsigned int eb_size>
	__device__ static inline bool isPadding(sparseGrid_type & sparseGrid, coord_type & coord, Mask_type (& enlargedBlock)[eb_size])
	{
		bool isPadding_ = false;
		for (int d=0; d<dim; ++d)
		{
			auto nPlusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, 1);
			auto nMinusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, -1);
			typename std::remove_all_extents<Mask_type>::type neighbourPlus = enlargedBlock[nPlusId];
			typename std::remove_all_extents<Mask_type>::type neighbourMinus = enlargedBlock[nMinusId];
			isPadding_ = isPadding_ || (!sparseGrid.exist(neighbourPlus));
			isPadding_ = isPadding_ || (!sparseGrid.exist(neighbourMinus));
			if (isPadding_) break;
		}

		return isPadding_;
	}

	/*! \brief given a coordinate give the neighborhood chunk position and the offset in the neighborhood chunk
	 *
	 *
	 */
	__device__ static inline bool getNNindex_offset()
	{
		return false;
	}
};

template<unsigned int n_it>
struct arr_ptr
{
	void * ptr[n_it];
};

template<unsigned int n_it,unsigned int n_prp>
struct arr_arr_ptr
{
	void * ptr[n_it][n_prp+1];
};

template<typename copy_type, unsigned int nprp, unsigned int prp_val, unsigned int prp_id>
struct meta_copy_block
{
	template<typename dataBuffer_type>
	__device__ __host__ static void copy(void * (& data_ptr)[nprp], dataBuffer_type & dataBuff, unsigned int ppos, unsigned int dataBlockPos, unsigned int offset, unsigned int n_pnt)
	{
		((copy_type *)data_ptr[prp_id])[ppos] = dataBuff.template get<prp_val>(dataBlockPos)[offset];
	}

	template<typename dataBuffer_type>
	__device__ __host__ static void copy_inv(arr_arr_ptr<1,nprp> & data_ptr, dataBuffer_type & dataBuff, unsigned int ppos, unsigned int dataBlockPos, unsigned int offset, unsigned int n_pnt)
	{
		dataBuff.template get<prp_val>(dataBlockPos)[offset] = ((copy_type *)data_ptr.ptr[0][prp_id])[ppos];
	}
};

template<typename copy_type, unsigned int nprp, unsigned int prp_val, unsigned int prp_id, unsigned int N1>
struct meta_copy_block<copy_type[N1],nprp,prp_val,prp_id>
{
	template<typename dataBuffer_type>
	__device__ __host__ static void copy(void * (& data_ptr)[nprp], dataBuffer_type & dataBuff, unsigned int ppos, unsigned int dataBlockPos, unsigned int offset, unsigned int n_pnt)
	{
		for (int i = 0 ; i < N1 ; i++)
		{
			((copy_type *)data_ptr[prp_id])[ppos+i*n_pnt] = dataBuff.template get<prp_val>(dataBlockPos)[i][offset];
		}
	}

	template<typename dataBuffer_type>
	__device__ __host__ static void copy_inv(arr_arr_ptr<1,nprp> & data_ptr, dataBuffer_type & dataBuff, unsigned int ppos, unsigned int dataBlockPos, unsigned int offset, unsigned int n_pnt)
	{
		for (int i = 0 ; i < N1 ; i++)
		{
			dataBuff.template get<prp_val>(dataBlockPos)[i][offset] = ((copy_type *)data_ptr.ptr[0][prp_id])[ppos+i*n_pnt];
		}
	}
};

template<typename copy_type, unsigned int nprp, unsigned int prp_val, unsigned int prp_id, unsigned int N1, unsigned int N2>
struct meta_copy_block<copy_type[N1][N2],nprp,prp_val,prp_id>
{
	template<typename dataBuffer_type>
	__device__ __host__ static void copy(void * (& data_ptr)[nprp], dataBuffer_type & dataBuff, unsigned int ppos, unsigned int dataBlockPos, unsigned int offset, unsigned int n_pnt)
	{
		for (int i = 0 ; i < N1 ; i++)
		{
			for (int j = 0 ; j < N2 ; j++)
			{
				((copy_type *)data_ptr[prp_id])[ppos + (i*N2 + j)*n_pnt] = dataBuff.template get<prp_val>(dataBlockPos)[i][j][offset];
			}
		}
	}

	template<typename dataBuffer_type>
	__device__ __host__ static void copy_inv(arr_arr_ptr<1,nprp> & data_ptr, dataBuffer_type & dataBuff, unsigned int ppos, unsigned int dataBlockPos, unsigned int offset, unsigned int n_pnt)
	{
		for (int i = 0 ; i < N1 ; i++)
		{
			for (int j = 0 ; j < N2 ; j++)
			{
				dataBuff.template get<prp_val>(dataBlockPos)[i][j][offset] = ((copy_type *)data_ptr.ptr[0][prp_id])[ppos + (i*N2 + j)*n_pnt];
			}
		}
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to calculate the size to pack a point
 *
 * \tparam prp set for properties
 *
 */
template<typename AggregateT, typename dataBuffer_type, int ... prp>
struct sparsegridgpu_pack_impl
{
	typedef typename to_boost_vmpl<prp...>::type vprp;

	//! position of the block
	unsigned int dataBlockPos;

	//! offset
	unsigned int offset;

	//! data buffer
	dataBuffer_type & dataBuff;

	//! point
	unsigned int ppos;

	//! data pointer
	void * (& data_ptr)[sizeof...(prp)+1];

	//! Number of points to pack
	unsigned int n_pnt;

	/*! \brief constructor
	 *
	 */
	__device__ __host__ inline sparsegridgpu_pack_impl(unsigned int dataBlockPos,
								   unsigned int offset,
								   dataBuffer_type & dataBuff,
								   unsigned int ppos,
								   void * (& data_ptr)[sizeof...(prp)+1],
								   unsigned int n_pnt)
	:dataBlockPos(dataBlockPos),offset(offset),dataBuff(dataBuff),ppos(ppos),data_ptr(data_ptr),n_pnt(n_pnt)
	{};

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<vprp,T>::type prp_cp;

		// Remove the reference from the type to copy
		typedef typename boost::mpl::at<typename AggregateT::type,prp_cp>::type pack_type;

		meta_copy_block<pack_type,sizeof...(prp)+1,prp_cp::value,T::value>::copy(data_ptr,dataBuff,ppos,dataBlockPos,offset,n_pnt);
	}
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to calculate the size to pack a point
 *
 * \tparam prp set for properties
 *
 */
template<typename AggregateT, typename dataBuffer_type, int ... prp>
struct sparsegridgpu_unpack_impl
{
	typedef typename to_boost_vmpl<prp...>::type vprp;

	//! position of the block
	unsigned int dataBlockPos;

	//! offset
	unsigned int offset;

	//! data buffer
	dataBuffer_type & dataBuff;

	//! point
	unsigned int ppos;

	//! data pointer
	arr_arr_ptr<1,sizeof...(prp)> & data_ptr;

	//! Number of points to pack
	unsigned int n_pnt;

	/*! \brief constructor
	 *
	 */
	__device__ __host__ inline sparsegridgpu_unpack_impl(unsigned int dataBlockPos,
								   unsigned int offset,
								   dataBuffer_type & dataBuff,
								   unsigned int ppos,
								   arr_arr_ptr<1,sizeof...(prp)> & data_ptr,
								   unsigned int n_pnt)
	:dataBlockPos(dataBlockPos),offset(offset),dataBuff(dataBuff),ppos(ppos),data_ptr(data_ptr),n_pnt(n_pnt)
	{};

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<vprp,T>::type prp_cp;

		// Remove the reference from the type to copy
		typedef typename boost::mpl::at<typename AggregateT::type,prp_cp>::type pack_type;

		meta_copy_block<pack_type,sizeof...(prp),prp_cp::value,T::value>::copy_inv(data_ptr,dataBuff,ppos,dataBlockPos,offset,n_pnt);
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to calculate the size to pack a point
 *
 * \tparam prp set for properties
 *
 */
template<typename AggregateT, int ... prp>
struct sparsegridgpu_pack_request
{
	typedef typename to_boost_vmpl<prp...>::type vprp;

	//! point size
	size_t point_size;

	/*! \brief constructor
	 *
	 */
	inline sparsegridgpu_pack_request()
	:point_size(0)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<vprp,T>::type prp_cp;

		// Remove the reference from the type to copy
		typedef typename boost::mpl::at<typename AggregateT::type,prp_cp>::type pack_type;

		point_size += sizeof(pack_type);
	}
};

template<unsigned int edgeSize, unsigned int dim>
inline __device__ unsigned int coordToLin(const unsigned int (&coord)[dim], const unsigned int paddingSize = 0)
{
	unsigned int linId = coord[dim - 1];
	for (int d = dim - 2; d >= 0; --d)
	{
		linId *= edgeSize + 2 * paddingSize;
		linId += coord[d];
	}
	return linId;
}


template<unsigned int edgeSize, typename CoordT, unsigned int dim>
inline __device__ unsigned int coordToLin(const grid_key_dx<dim, CoordT> &coord, const unsigned int paddingSize = 0)
{
	unsigned int linId = coord.get(dim - 1);
	for (int d = dim - 2; d >= 0; --d)
	{
		linId *= edgeSize + 2 * paddingSize;
		linId += coord.get(d);
	}
	return linId;
}

template <typename CoordT,unsigned int dim>
inline __device__ unsigned int coordToLin(const grid_key_dx<dim, CoordT> & coord, grid_key_dx<dim, int> & blockDimensions)
{
	unsigned int linId = coord.get(dim - 1);
	for (int d = dim - 2; d >= 0; --d)
	{
		linId *= blockDimensions.get(d);
		linId += coord.get(d);
	}
	return linId;
}



template<unsigned int edgeSize, unsigned int dim>
inline __device__ void linToCoordWithOffset(const unsigned int linId, const unsigned int offset, unsigned int (&coord)[dim])
{
	unsigned int linIdTmp = linId;
	for (unsigned int d = 0; d < dim; ++d)
	{
		coord[d] = linIdTmp % edgeSize;
		coord[d] += offset;
		linIdTmp /= edgeSize;
	}
}

template<unsigned int edgeSize, unsigned int dim>
inline __device__ void linToCoord(const unsigned int linId, unsigned int (&coord)[dim])
{
	unsigned int linIdTmp = linId;
	for (unsigned int d = 0; d < dim; ++d)
	{
		coord[d] = linIdTmp % edgeSize;
		linIdTmp /= edgeSize;
	}
}

constexpr int gt = 0;
constexpr int nt = 1;

template<unsigned int nLoop, unsigned int dim, typename AggregateBlockT, unsigned int pMask , unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl
{
	template<typename AggrWrapperT,
	         typename SharedPtrT,
	         typename vector_type,
	         typename vector_type2,
	         typename blockMapType,
	         typename AggrBck>
	__device__ static inline void load(const AggrWrapperT &block,
							SharedPtrT * sharedRegionPtr,
							const vector_type & ghostLayerToThreadsMapping,
							const vector_type2 & nn_blocks,
							const blockMapType & blockMap,
							unsigned int stencilSupportRadius,
							unsigned int ghostLayerSize,
							const unsigned int blockId,
							AggrBck & bck)
	{
		printf("Error to implement loadGhostBlock_impl with nLoop=%d \n",nLoop);
	}
};

template<unsigned int dim, typename AggregateBlockT, unsigned int pMask , unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl<1,dim,AggregateBlockT,pMask,p,ct_params,blockEdgeSize>
{
	template<typename AggrWrapperT,
			 typename SharedPtrT,
			 typename vector_type,
			 typename vector_type2,
			 typename blockMapType,
			 typename AggrBck>
	__device__ static inline void load(const AggrWrapperT &block,
							SharedPtrT * sharedRegionPtr,
							const vector_type & ghostLayerToThreadsMapping,
							const vector_type2 & nn_blocks,
							const blockMapType & blockMap,
							unsigned int stencilSupportRadius,
							unsigned int ghostLayerSize,
							const unsigned int blockIdPos,
							AggrBck & bck)
	{
			typedef ScalarTypeOf<AggregateBlockT, p> ScalarT;

			const int pos = threadIdx.x % ghostLayerSize;

			__shared__ int neighboursPos[ct_params::nNN];

			const unsigned int edge = blockEdgeSize + 2*stencilSupportRadius;
			short int neighbourNum = ghostLayerToThreadsMapping.template get<nt>(pos);

			// Convert pos into a linear id accounting for the inner domain offsets
			const unsigned int linId = ghostLayerToThreadsMapping.template get<gt>(pos);
			// Now get linear offset wrt the first element of the block

			int ctr = linId;
			unsigned int acc = 1;
			unsigned int offset = 0;
			for (int i = 0; i < dim; ++i)
			{
				int v = (ctr % edge) - stencilSupportRadius;
				v = (v < 0)?(v + blockEdgeSize):v;
				v = (v >= blockEdgeSize)?v-blockEdgeSize:v;
				offset += v*acc;
				ctr /= edge;
				acc *= blockEdgeSize;
			}

			// Convert pos into a linear id accounting for the ghost offsets
			unsigned int coord[dim];
			linToCoordWithOffset<blockEdgeSize>(threadIdx.x, stencilSupportRadius, coord);
			const unsigned int linId2 = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);

			unsigned int nnb = nn_blocks.template get<0>(blockIdPos*ct_params::nNN + (threadIdx.x % ct_params::nNN));

			if (threadIdx.x < ct_params::nNN)
			{
				neighboursPos[threadIdx.x] = nnb;
			}

			__syncthreads();

			// Actually load the data into the shared region
			auto nPos = neighboursPos[neighbourNum];

			auto gdata = blockMap.template get_ele<p>(nPos)[offset];

			// Actually load the data into the shared region
			//ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
			auto bdata = block.template get<p>()[threadIdx.x];

			auto bmask = block.template get<pMask>()[threadIdx.x];
			auto gmask = blockMap.template get_ele<pMask>(nPos)[offset];

			if (bmask == 0)	{ set_compile_condition<pMask != p>::template set<p>(bdata,bck);}
			if (gmask == 0)	{ set_compile_condition<pMask != p>::template set<p>(gdata,bck);}

			sharedRegionPtr[linId] = gdata;
			sharedRegionPtr[linId2] = bdata;
	}
};

template<unsigned int dim, typename AggregateBlockT, unsigned int pMask , unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl<2,dim,AggregateBlockT,pMask,p,ct_params,blockEdgeSize>
{
    template<typename AggrWrapperT,
    		 typename SharedPtrT,
    		 typename vector_type,
    		 typename vector_type2,
    		 typename blockMapType,
    		 typename AggrBck>
    __device__ static inline void load(const AggrWrapperT &block,
                                       SharedPtrT * sharedRegionPtr,
                                       const vector_type & ghostLayerToThreadsMapping,
                                       const vector_type2 & nn_blocks,
                                       const blockMapType & blockMap,
                                       unsigned int stencilSupportRadius,
                                       unsigned int ghostLayerSize,
                                       const unsigned int blockIdPos,
                                       AggrBck & bck)
    {
        typedef ScalarTypeOf<AggregateBlockT, p> ScalarT;

        const int pos = threadIdx.x % ghostLayerSize;
        const int pos_d1 = (threadIdx.x + blockDim.x) % ghostLayerSize;

        __shared__ int neighboursPos[ct_params::nNN];

        const unsigned int edge = blockEdgeSize + 2*stencilSupportRadius;
        short int neighbourNum = ghostLayerToThreadsMapping.template get<nt>(pos);
        short int neighbourNum2 = ghostLayerToThreadsMapping.template get<nt>(pos_d1);

        // Convert pos into a linear id accounting for the inner domain offsets
        const unsigned int linId = ghostLayerToThreadsMapping.template get<gt>(pos);
        const unsigned int linId2 = ghostLayerToThreadsMapping.template get<gt>(pos_d1);
        // Now get linear offset wrt the first element of the block

        int ctr = linId;
        int ctr2 = linId2;
        unsigned int acc = 1;
        unsigned int offset = 0;
        unsigned int offset2 = 0;
        for (int i = 0; i < dim; ++i)
        {
            int v = (ctr % edge) - stencilSupportRadius;
            int v2 = (ctr2 % edge) - stencilSupportRadius;
            v = (v < 0)?(v + blockEdgeSize):v;
            v2 = (v2 < 0)?(v2 + blockEdgeSize):v2;
            v = (v >= blockEdgeSize)?v-blockEdgeSize:v;
            v2 = (v2 >= blockEdgeSize)?v2-blockEdgeSize:v2;
            offset += v*acc;
            offset2 += v2*acc;
            ctr /= edge;
            ctr2 /= edge;
            acc *= blockEdgeSize;
        }

        // Convert pos into a linear id accounting for the ghost offsets
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(threadIdx.x, stencilSupportRadius, coord);
        const int linId_b = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
//        const unsigned int linId_b = shift_position<dim,blockEdgeSize>::shift(threadIdx.x,stencilSupportRadius);

//        printf("AAA %d %d \n",linId_b,linId_b_test);

        unsigned int nnb = nn_blocks.template get<0>(blockIdPos*ct_params::nNN + (threadIdx.x % ct_params::nNN));

        if (threadIdx.x < ct_params::nNN)
        {
            neighboursPos[threadIdx.x] = nnb;
        }

        __syncthreads();

        // Actually load the data into the shared region
        auto nPos = neighboursPos[neighbourNum];
        auto nPos2 = neighboursPos[neighbourNum2];

        auto gdata = blockMap.template get_ele<p>(nPos)[offset];
        auto gdata2 = blockMap.template get_ele<p>(nPos2)[offset2];

        // Actually load the data into the shared region
        //ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
        auto bdata = block.template get<p>()[threadIdx.x];

        auto gmask = blockMap.template get_ele<pMask>(nPos)[offset];
        auto gmask2 = blockMap.template get_ele<pMask>(nPos2)[offset2];

        // Actually load the data into the shared region
        //ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
        auto bmask = block.template get<pMask>()[threadIdx.x];

		if (bmask == 0)	{set_compile_condition<pMask != p>::template set<p>(bdata,bck);}
		if (gmask == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata,bck);}
		if (gmask2 == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata2,bck);}

        sharedRegionPtr[linId] = gdata;
        sharedRegionPtr[linId2] = gdata2;
        sharedRegionPtr[linId_b] = bdata;
    }
};

template<unsigned int dim, typename AggregateBlockT, unsigned int pMask , unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl<3,dim,AggregateBlockT,pMask,p,ct_params,blockEdgeSize>
{
    template<typename AggrWrapperT,
             typename SharedPtrT,
             typename vector_type,
             typename vector_type2,
             typename blockMapType,
             typename AggrBck>
    __device__ static inline void load(const AggrWrapperT &block,
                                       SharedPtrT * sharedRegionPtr,
                                       const vector_type & ghostLayerToThreadsMapping,
                                       const vector_type2 & nn_blocks,
                                       const blockMapType & blockMap,
                                       unsigned int stencilSupportRadius,
                                       unsigned int ghostLayerSize,
                                       const unsigned int blockIdPos,
                                       AggrBck & bck)
    {
        typedef ScalarTypeOf<AggregateBlockT, p> ScalarT;

        const int pos = threadIdx.x % ghostLayerSize;
        const int pos_d1 = (threadIdx.x + 2*blockDim.x) % ghostLayerSize;

        __shared__ int neighboursPos[ct_params::nNN];

        const unsigned int edge = blockEdgeSize + 2*stencilSupportRadius;
        short int neighbourNum = ghostLayerToThreadsMapping.template get<nt>(pos);
        short int neighbourNum2 = ghostLayerToThreadsMapping.template get<nt>(pos + blockDim.x);
        short int neighbourNum3 = ghostLayerToThreadsMapping.template get<nt>(pos_d1);

        // Convert pos into a linear id accounting for the inner domain offsets
        const unsigned int linId = ghostLayerToThreadsMapping.template get<gt>(pos);
        const unsigned int linId2 = ghostLayerToThreadsMapping.template get<gt>(pos + blockDim.x);
        const unsigned int linId3 = ghostLayerToThreadsMapping.template get<gt>(pos_d1);
        // Now get linear offset wrt the first element of the block

        int ctr = linId;
        int ctr2 = linId2;
        int ctr3 = linId3;
        unsigned int acc = 1;
        unsigned int offset = 0;
        unsigned int offset2 = 0;
        unsigned int offset3 = 0;
        for (int i = 0; i < dim; ++i)
        {
            int v = (ctr % edge) - stencilSupportRadius;
            int v2 = (ctr2 % edge) - stencilSupportRadius;
            int v3 = (ctr3 % edge) - stencilSupportRadius;
            v = (v < 0)?(v + blockEdgeSize):v;
            v2 = (v2 < 0)?(v2 + blockEdgeSize):v2;
            v3 = (v3 < 0)?(v3 + blockEdgeSize):v3;
            v = (v >= blockEdgeSize)?v-blockEdgeSize:v;
            v2 = (v2 >= blockEdgeSize)?v2-blockEdgeSize:v2;
            v3 = (v3 >= blockEdgeSize)?v3-blockEdgeSize:v3;
            offset += v*acc;
            offset2 += v2*acc;
            offset3 += v3*acc;
            ctr /= edge;
            ctr2 /= edge;
            ctr3 /= edge;
            acc *= blockEdgeSize;
        }

        // Convert pos into a linear id accounting for the ghost offsets
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(threadIdx.x, stencilSupportRadius, coord);
        const int linId_b = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
//        const unsigned int linId_b = shift_position<dim,blockEdgeSize>::shift(threadIdx.x,stencilSupportRadius);

//        printf("AAA %d %d \n",linId_b,linId_b_test);

        unsigned int nnb = nn_blocks.template get<0>(blockIdPos*ct_params::nNN + (threadIdx.x % ct_params::nNN));

        if (threadIdx.x < ct_params::nNN)
        {
            neighboursPos[threadIdx.x] = nnb;
        }

        __syncthreads();

        // Actually load the data into the shared region
        auto nPos = neighboursPos[neighbourNum];
        auto nPos2 = neighboursPos[neighbourNum2];
        auto nPos3 = neighboursPos[neighbourNum3];

        auto gdata = blockMap.template get_ele<p>(nPos)[offset];
        auto gdata2 = blockMap.template get_ele<p>(nPos2)[offset2];
        auto gdata3 = blockMap.template get_ele<p>(nPos3)[offset3];

        // Actually load the data into the shared region
        //ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
        auto bdata = block.template get<p>()[threadIdx.x];

        auto gmask = blockMap.template get_ele<pMask>(nPos)[offset];
        auto gmask2 = blockMap.template get_ele<pMask>(nPos2)[offset2];
        auto gmask3 = blockMap.template get_ele<pMask>(nPos3)[offset3];

        // Actually load the data into the shared region
        //ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
        auto bmask = block.template get<pMask>()[threadIdx.x];

		if (bmask == 0)	{set_compile_condition<pMask != p>::template set<p>(bdata,bck);}
		if (gmask == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata,bck);}
		if (gmask2 == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata2,bck);}
		if (gmask3 == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata3,bck);}

        sharedRegionPtr[linId] = gdata;
        sharedRegionPtr[linId2] = gdata2;
        sharedRegionPtr[linId3] = gdata3;
        sharedRegionPtr[linId_b] = bdata;
    }
};

template<unsigned int dim, typename AggregateBlockT, unsigned int pMask , unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl<7,dim,AggregateBlockT,pMask,p,ct_params,blockEdgeSize>
{
    template<typename AggrWrapperT,
             typename SharedPtrT,
             typename vector_type,
             typename vector_type2,
             typename blockMapType,
             typename AggrBck>
    __device__ static inline void load(const AggrWrapperT &block,
                                       SharedPtrT * sharedRegionPtr,
                                       const vector_type & ghostLayerToThreadsMapping,
                                       const vector_type2 & nn_blocks,
                                       const blockMapType & blockMap,
                                       unsigned int stencilSupportRadius,
                                       unsigned int ghostLayerSize,
                                       const unsigned int blockIdPos,
                                       AggrBck & bck)
    {
        typedef ScalarTypeOf<AggregateBlockT, p> ScalarT;

        const int pos = threadIdx.x % ghostLayerSize;
        const int pos_d1 = (threadIdx.x + 6*blockDim.x) % ghostLayerSize;

        __shared__ int neighboursPos[ct_params::nNN];

        const unsigned int edge = blockEdgeSize + 2*stencilSupportRadius;
        short int neighbourNum = ghostLayerToThreadsMapping.template get<nt>(pos);
        short int neighbourNum2 = ghostLayerToThreadsMapping.template get<nt>(pos + blockDim.x);
        short int neighbourNum3 = ghostLayerToThreadsMapping.template get<nt>(pos + 2*blockDim.x);
        short int neighbourNum4 = ghostLayerToThreadsMapping.template get<nt>(pos + 3*blockDim.x);
        short int neighbourNum5 = ghostLayerToThreadsMapping.template get<nt>(pos + 4*blockDim.x);
        short int neighbourNum6 = ghostLayerToThreadsMapping.template get<nt>(pos + 5*blockDim.x);
        short int neighbourNum7 = ghostLayerToThreadsMapping.template get<nt>(pos_d1);

        // Convert pos into a linear id accounting for the inner domain offsets
        const unsigned int linId = ghostLayerToThreadsMapping.template get<gt>(pos);
        const unsigned int linId2 = ghostLayerToThreadsMapping.template get<gt>(pos + blockDim.x);
        const unsigned int linId3 = ghostLayerToThreadsMapping.template get<gt>(pos + 2*blockDim.x);
        const unsigned int linId4 = ghostLayerToThreadsMapping.template get<gt>(pos + 3*blockDim.x);
        const unsigned int linId5 = ghostLayerToThreadsMapping.template get<gt>(pos + 4*blockDim.x);
        const unsigned int linId6 = ghostLayerToThreadsMapping.template get<gt>(pos + 5*blockDim.x);
        const unsigned int linId7 = ghostLayerToThreadsMapping.template get<gt>(pos_d1);
        // Now get linear offset wrt the first element of the block

        int ctr = linId;
        int ctr2 = linId2;
        int ctr3 = linId3;
        int ctr4 = linId4;
        int ctr5 = linId5;
        int ctr6 = linId6;
        int ctr7 = linId7;
        unsigned int acc = 1;
        unsigned int offset = 0;
        unsigned int offset2 = 0;
        unsigned int offset3 = 0;
        unsigned int offset4 = 0;
        unsigned int offset5 = 0;
        unsigned int offset6 = 0;
        unsigned int offset7 = 0;
        for (int i = 0; i < dim; ++i)
        {
            int v = (ctr % edge) - stencilSupportRadius;
            int v2 = (ctr2 % edge) - stencilSupportRadius;
            int v3 = (ctr3 % edge) - stencilSupportRadius;
            int v4 = (ctr4 % edge) - stencilSupportRadius;
            int v5 = (ctr5 % edge) - stencilSupportRadius;
            int v6 = (ctr6 % edge) - stencilSupportRadius;
            int v7 = (ctr7 % edge) - stencilSupportRadius;
            v = (v < 0)?(v + blockEdgeSize):v;
            v2 = (v2 < 0)?(v2 + blockEdgeSize):v2;
            v3 = (v3 < 0)?(v3 + blockEdgeSize):v3;
            v4 = (v4 < 0)?(v4 + blockEdgeSize):v4;
            v5 = (v5 < 0)?(v5 + blockEdgeSize):v5;
            v6 = (v6 < 0)?(v6 + blockEdgeSize):v6;
            v7 = (v7 < 0)?(v7 + blockEdgeSize):v7;
            v = (v >= blockEdgeSize)?v-blockEdgeSize:v;
            v2 = (v2 >= blockEdgeSize)?v2-blockEdgeSize:v2;
            v3 = (v3 >= blockEdgeSize)?v3-blockEdgeSize:v3;
            v4 = (v4 >= blockEdgeSize)?v4-blockEdgeSize:v4;
            v5 = (v5 >= blockEdgeSize)?v5-blockEdgeSize:v5;
            v6 = (v6 >= blockEdgeSize)?v6-blockEdgeSize:v6;
            v7 = (v7 >= blockEdgeSize)?v7-blockEdgeSize:v7;
            offset += v*acc;
            offset2 += v2*acc;
            offset3 += v3*acc;
            offset4 += v4*acc;
            offset5 += v5*acc;
            offset6 += v6*acc;
            offset7 += v7*acc;
            ctr /= edge;
            ctr2 /= edge;
            ctr3 /= edge;
            ctr4 /= edge;
            ctr5 /= edge;
            ctr6 /= edge;
            ctr7 /= edge;
            acc *= blockEdgeSize;
        }

        // Convert pos into a linear id accounting for the ghost offsets
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(threadIdx.x, stencilSupportRadius, coord);
        const int linId_b = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
//        const unsigned int linId_b = shift_position<dim,blockEdgeSize>::shift(threadIdx.x,stencilSupportRadius);

//        printf("AAA %d %d \n",linId_b,linId_b_test);

        unsigned int nnb = nn_blocks.template get<0>(blockIdPos*ct_params::nNN + (threadIdx.x % ct_params::nNN));

        if (threadIdx.x < ct_params::nNN)
        {
            neighboursPos[threadIdx.x] = nnb;
        }

        __syncthreads();

        // Actually load the data into the shared region
        auto nPos = neighboursPos[neighbourNum];
        auto nPos2 = neighboursPos[neighbourNum2];
        auto nPos3 = neighboursPos[neighbourNum3];
        auto nPos4 = neighboursPos[neighbourNum4];
        auto nPos5 = neighboursPos[neighbourNum5];
        auto nPos6 = neighboursPos[neighbourNum6];
        auto nPos7 = neighboursPos[neighbourNum7];

        auto gdata = blockMap.template get_ele<p>(nPos)[offset];
        auto gdata2 = blockMap.template get_ele<p>(nPos2)[offset2];
        auto gdata3 = blockMap.template get_ele<p>(nPos3)[offset3];
        auto gdata4 = blockMap.template get_ele<p>(nPos4)[offset4];
        auto gdata5 = blockMap.template get_ele<p>(nPos5)[offset5];
        auto gdata6 = blockMap.template get_ele<p>(nPos6)[offset6];
        auto gdata7 = blockMap.template get_ele<p>(nPos7)[offset7];

        // Actually load the data into the shared region
        //ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
        auto bdata = block.template get<p>()[threadIdx.x];

        auto bmask = block.template get<pMask>()[threadIdx.x];
        auto gmask = blockMap.template get_ele<pMask>(nPos)[offset];
        auto gmask2 = blockMap.template get_ele<pMask>(nPos2)[offset2];
        auto gmask3 = blockMap.template get_ele<pMask>(nPos3)[offset3];
        auto gmask4 = blockMap.template get_ele<pMask>(nPos4)[offset4];
        auto gmask5 = blockMap.template get_ele<pMask>(nPos5)[offset5];
        auto gmask6 = blockMap.template get_ele<pMask>(nPos6)[offset6];
        auto gmask7 = blockMap.template get_ele<pMask>(nPos7)[offset7];

		if (bmask == 0)	{set_compile_condition<pMask != p>::template set<p>(bdata,bck);}
		if (gmask == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata,bck);}
		if (gmask2 == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata2,bck);}
		if (gmask3 == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata3,bck);}
		if (gmask4 == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata4,bck);}
		if (gmask5 == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata5,bck);}
		if (gmask6 == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata6,bck);}
		if (gmask7 == 0)	{set_compile_condition<pMask != p>::template set<p>(gdata7,bck);}

        sharedRegionPtr[linId] = gdata;
        sharedRegionPtr[linId2] = gdata2;
        sharedRegionPtr[linId3] = gdata3;
        sharedRegionPtr[linId4] = gdata4;
        sharedRegionPtr[linId5] = gdata5;
        sharedRegionPtr[linId6] = gdata6;
        sharedRegionPtr[linId7] = gdata7;
        sharedRegionPtr[linId_b] = bdata;
    }
};

#endif /* SPARSEGRIDGPU_KER_UTIL_HPP_ */
