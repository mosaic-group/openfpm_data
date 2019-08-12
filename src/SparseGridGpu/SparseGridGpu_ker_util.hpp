/*
 * SparseGridGpu_ker_util.hpp
 *
 *  Created on: Aug 7, 2019
 *      Author: i-bird
 */

#ifndef SPARSEGRIDGPU_KER_UTIL_HPP_
#define SPARSEGRIDGPU_KER_UTIL_HPP_

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


#endif /* SPARSEGRIDGPU_KER_UTIL_HPP_ */
