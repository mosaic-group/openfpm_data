/*
 * SparseGridGpu_ker_util.hpp
 *
 *  Created on: Aug 7, 2019
 *      Author: i-bird
 */

#ifndef SPARSEGRIDGPU_KER_UTIL_HPP_
#define SPARSEGRIDGPU_KER_UTIL_HPP_

#include "util/variadic_to_vmpl.hpp"

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

constexpr int gt = 0;
constexpr int nt = 1;

template<unsigned int nLoop, unsigned int dim, typename AggregateBlockT, unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl
{
	template<typename AggrWrapperT, typename SharedPtrT, typename vector_type, typename vector_type2, typename blockMapType>
	__device__ static inline void load(const AggrWrapperT &block,
							SharedPtrT * sharedRegionPtr,
							const vector_type & ghostLayerToThreadsMapping,
							const vector_type2 & nn_blocks,
							const blockMapType & blockMap,
							unsigned int stencilSupportRadius,
							unsigned int ghostLayerSize,
							const unsigned int blockId)
	{
		printf("Error to implement loadGhostBlock_impl with nLoop=%d \n",nLoop);
	}
};

template<unsigned int dim, typename AggregateBlockT, unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl<1,dim,AggregateBlockT,p,ct_params,blockEdgeSize>
{
	template<typename AggrWrapperT, typename SharedPtrT, typename vector_type, typename vector_type2, typename blockMapType>
	__device__ static inline void load(const AggrWrapperT &block,
							SharedPtrT * sharedRegionPtr,
							const vector_type & ghostLayerToThreadsMapping,
							const vector_type2 & nn_blocks,
							const blockMapType & blockMap,
							unsigned int stencilSupportRadius,
							unsigned int ghostLayerSize,
							const unsigned int blockIdPos)
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
			const int linId2 = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);

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

			sharedRegionPtr[linId] = gdata;
			sharedRegionPtr[linId2] = bdata;
	}
};

template<unsigned int dim, typename AggregateBlockT, unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl<2,dim,AggregateBlockT,p,ct_params,blockEdgeSize>
{
    template<typename AggrWrapperT, typename SharedPtrT, typename vector_type, typename vector_type2, typename blockMapType>
    __device__ static inline void load(const AggrWrapperT &block,
                                       SharedPtrT * sharedRegionPtr,
                                       const vector_type & ghostLayerToThreadsMapping,
                                       const vector_type2 & nn_blocks,
                                       const blockMapType & blockMap,
                                       unsigned int stencilSupportRadius,
                                       unsigned int ghostLayerSize,
                                       const unsigned int blockIdPos)
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

        sharedRegionPtr[linId] = gdata;
        sharedRegionPtr[linId2] = gdata2;
        sharedRegionPtr[linId_b] = bdata;
    }
};

template<unsigned int dim, typename AggregateBlockT, unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl<3,dim,AggregateBlockT,p,ct_params,blockEdgeSize>
{
    template<typename AggrWrapperT, typename SharedPtrT, typename vector_type, typename vector_type2, typename blockMapType>
    __device__ static inline void load(const AggrWrapperT &block,
                                       SharedPtrT * sharedRegionPtr,
                                       const vector_type & ghostLayerToThreadsMapping,
                                       const vector_type2 & nn_blocks,
                                       const blockMapType & blockMap,
                                       unsigned int stencilSupportRadius,
                                       unsigned int ghostLayerSize,
                                       const unsigned int blockIdPos)
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

        sharedRegionPtr[linId] = gdata;
        sharedRegionPtr[linId2] = gdata2;
        sharedRegionPtr[linId3] = gdata3;
        sharedRegionPtr[linId_b] = bdata;
    }
};

template<unsigned int dim, typename AggregateBlockT, unsigned int p, typename ct_params, unsigned int blockEdgeSize>
struct loadGhostBlock_impl<7,dim,AggregateBlockT,p,ct_params,blockEdgeSize>
{
    template<typename AggrWrapperT, typename SharedPtrT, typename vector_type, typename vector_type2, typename blockMapType>
    __device__ static inline void load(const AggrWrapperT &block,
                                       SharedPtrT * sharedRegionPtr,
                                       const vector_type & ghostLayerToThreadsMapping,
                                       const vector_type2 & nn_blocks,
                                       const blockMapType & blockMap,
                                       unsigned int stencilSupportRadius,
                                       unsigned int ghostLayerSize,
                                       const unsigned int blockIdPos)
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
