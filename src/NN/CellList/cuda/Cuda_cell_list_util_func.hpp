/*
 * Cuda_cell_list_util_func.hpp
 *
 *  Created on: Jun 17, 2018
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CUDA_CUDA_CELL_LIST_UTIL_FUNC_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CUDA_CUDA_CELL_LIST_UTIL_FUNC_HPP_

#include <boost/integer/integer_mask.hpp>
#include <Vector/map_vector_sparse.hpp>


template<unsigned int dim, typename ids_type, typename transform_type>
struct cid_
{
	static inline __device__ __host__ unsigned int get_cid(
		openfpm::array<ids_type,1> & numCellDiv,
		ids_type * e)
	{
		unsigned int id = e[dim-1];

#pragma unroll
		for (int i = 1; i >= 0 ; i-- )
		{id = e[i] + numCellDiv[i]*id;}

		return id;
	}

	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,dim> & numCellDiv,
		const grid_key_dx<1,unsigned int> & e)
	{
		unsigned int id = e.get(dim-1);

#pragma unroll
		for (int i = 1; i >= 0 ; i-- )
		{id = e.get(i) + numCellDiv[i]*id;}

		return id;
	}

	template<typename T>
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,dim> & numCellDiv,
		const openfpm::array<T,dim> & unitCellP2,
		const transform_type & pointTransform,
		const Point<dim,T> & p)
	{
		unsigned int id = p.get(dim-1) / unitCellP2[dim-1];

#pragma unroll
		for (int i = 1; i >= 0 ; i-- )
		{id = pointTransform.transform(p.get(i),i) / unitCellP2[i] + numCellDiv[i]*id;}

		return id;
	}
};

template<typename ids_type, typename transform_type>
struct cid_<1,ids_type, transform_type>
{
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,1> & numCellDiv,
		ids_type * e)
	{
		return e[0];
	}

	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,1> & numCellDiv,
		const grid_key_dx<1,ids_type> & e)
	{
		return e.get(0);
	}

	template<typename T>
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,1> & numCellDiv,
		const openfpm::array<T,1> & unitCellP2,
		const openfpm::array<ids_type,1> & cellPadDim,
		const transform_type & pointTransform,
		const Point<1,T> & p)
	{
		return pointTransform.transform(p.get(0),0) / unitCellP2[0] + cellPadDim[0];
	}

	template<typename T>
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,1> & numCellDiv,
		const openfpm::array<T,1> & unitCellP2,
		const openfpm::array<ids_type,1> & cellPadDim,
		const transform_type & pointTransform,
		const T * p,
		ids_type * e)
	{
		e[0] = openfpm::math::uint_floor(pointTransform.transform(p,0)/unitCellP2[0]) + cellPadDim[0];

		return e[0];
	}

	template<typename T>
	static inline __device__ __host__ grid_key_dx<2,ids_type> get_cid_key(
		const openfpm::array<T,1> & unitCellP2,
		const openfpm::array<ids_type,1> & cellPadDim,
		const transform_type & pointTransform,
		const Point<2,T> & p)
	{
		grid_key_dx<1,ids_type> e;

		e.set_d(0,openfpm::math::uint_floor(pointTransform.transform(p,0)/unitCellP2[0]) + cellPadDim[0]);

		return e;
	}

	template <typename U = unsigned int, typename sfinae=typename std::enable_if<std::is_same<ids_type,U>::value >::type >
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,1> & numCellDiv,
		const grid_key_dx<1,unsigned int> & e)
	{
		return e.get(0);
	}
};

template<typename ids_type, typename transform_type>
struct cid_<2,ids_type,transform_type>
{
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,2> & numCellDiv,
		ids_type * e)
	{
		return e[0] + numCellDiv[0] * e[1];
	}

	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,2> & numCellDiv,
		const grid_key_dx<2,ids_type> & e)
	{
		return e.get(0) + numCellDiv[0] * e.get(1);
	}

	template<typename T>
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,2> & numCellDiv,
		const openfpm::array<T,2> & unitCellP2,
		const openfpm::array<ids_type,2> & cellPadDim,
		const transform_type & pointTransform,
		const Point<2,T> & p)
	{

		return openfpm::math::uint_floor(pointTransform.transform(p,0)/unitCellP2[0]) + cellPadDim[0] +
				  (openfpm::math::uint_floor(pointTransform.transform(p,1)/unitCellP2[1]) + cellPadDim[1])*numCellDiv[0];
	}

	template<typename T>
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,2> & numCellDiv,
		const openfpm::array<T,2> & unitCellP2,
		const openfpm::array<ids_type,2> & cellPadDim,
		const transform_type & pointTransform,
		const T * p,
		ids_type * e)
	{
		e[0] = openfpm::math::uint_floor(pointTransform.transform(p,0)/unitCellP2[0]) + cellPadDim[0];
		e[1] = openfpm::math::uint_floor(pointTransform.transform(p,1)/unitCellP2[1]) + cellPadDim[1];

		return e[0] + e[1]*numCellDiv[0];
	}

	template<typename T>
	static inline __device__ __host__ grid_key_dx<2,ids_type> get_cid_key(
		const openfpm::array<T,2> & unitCellP2,
		const openfpm::array<ids_type,2> & cellPadDim,
		const transform_type & pointTransform,
		const Point<2,T> & p)
	{
		grid_key_dx<2,ids_type> e;

		e.set_d(0,openfpm::math::uint_floor(pointTransform.transform(p,0)/unitCellP2[0]) + cellPadDim[0]);
		e.set_d(1,openfpm::math::uint_floor(pointTransform.transform(p,1)/unitCellP2[1]) + cellPadDim[1]);

		return e;
	}

	template <typename U = unsigned int, typename sfinae=typename std::enable_if<std::is_same<ids_type,U>::value >::type >
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,2> & numCellDiv,
		const grid_key_dx<2,unsigned int> & e)
	{
		return e.get(0) + e.get(1)*numCellDiv[0];
	}
};


template<typename ids_type,typename transform_type>
struct cid_<3,ids_type,transform_type>
{

	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,3> & numCellDiv,
		const ids_type * e)
	{
		return e[0] + (e[1] + e[2]*numCellDiv[1])*numCellDiv[0];
	}

	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,3> & numCellDiv,
		const grid_key_dx<3,ids_type> & e)
	{
		return e.get(0) + (e.get(1) + e.get(2)*numCellDiv[1])*numCellDiv[0];
	}

	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,3> & numCellDiv,
		const openfpm::array<ids_type,3> & cellPadDim,
		const grid_key_dx<3,ids_type> & e)
	{
		return (e.get(0) + cellPadDim[0]) + ((e.get(1) + cellPadDim[1]) + (e.get(2) + cellPadDim[2])*numCellDiv[1])*numCellDiv[0];
	}

	template<typename T> static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,3> & numCellDiv,
		const openfpm::array<T,3> & unitCellP2,
		const openfpm::array<ids_type,3> & cellPadDim,
		const transform_type & pointTransform,
		const Point<3,T> & p)
	{
		return openfpm::math::uint_floor(pointTransform.transform(p,0)/unitCellP2[0]) + cellPadDim[0] +
			(openfpm::math::uint_floor(pointTransform.transform(p,1)/unitCellP2[1]) + cellPadDim[1] +
			(openfpm::math::uint_floor(pointTransform.transform(p,2)/unitCellP2[2]) + cellPadDim[2])*numCellDiv[1])*numCellDiv[0];
	}

	template<typename T>
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,3> & numCellDiv,
		const openfpm::array<T,3> & unitCellP2,
		const openfpm::array<ids_type,3> & cellPadDim,
		const transform_type & pointTransform,
		const T * p,
		ids_type * e)
	{
		e[0] = openfpm::math::uint_floor(pointTransform.transform(p,0)/unitCellP2[0]) + cellPadDim[0];
		e[1] = openfpm::math::uint_floor(pointTransform.transform(p,1)/unitCellP2[1]) + cellPadDim[1];
		e[2] = openfpm::math::uint_floor(pointTransform.transform(p,2)/unitCellP2[2]) + cellPadDim[2];

		return e[0] + (e[1] + e[2]*numCellDiv[1])*numCellDiv[0];
	}

	template<typename T>
	static inline __device__ __host__ grid_key_dx<3,ids_type> get_cid_key(
		const openfpm::array<T,3> & unitCellP2,
		const openfpm::array<ids_type,3> & cellPadDim,
		const transform_type & pointTransform,
		const Point<3,T> & p)
	{
		grid_key_dx<3,ids_type> e;

		e.set_d(0,openfpm::math::uint_floor(pointTransform.transform(p,0)/unitCellP2[0]) + cellPadDim[0]);
		e.set_d(1,openfpm::math::uint_floor(pointTransform.transform(p,1)/unitCellP2[1]) + cellPadDim[1]);
		e.set_d(2,openfpm::math::uint_floor(pointTransform.transform(p,2)/unitCellP2[2]) + cellPadDim[2]);

		return e;
	}

	template <typename U = unsigned int, typename sfinae=typename std::enable_if<std::is_same<ids_type,U>::value >::type >
	static inline __device__ __host__ unsigned int get_cid(
		const openfpm::array<ids_type,3> & numCellDiv,
		const grid_key_dx<3,unsigned int> & e)
	{
		return e.get(0) + (e.get(1) + e.get(2)*numCellDiv[1])*numCellDiv[0];
	}
};

/////////////////////////// THIS ONLY WORK IF NVCC IS COMPILING THIS //////////////////////////////////////////////////////

#ifdef __NVCC__

template<unsigned int dim, typename pos_type, typename ids_type, typename transform_type,
typename vector_pos_type, typename vector_cnt_type, typename vector_pids_type>
__global__ void fill_cellIndex_LocalIndex(
	openfpm::array<ids_type,dim> numCellDiv,
	openfpm::array<pos_type,dim> unitCellP2,
	openfpm::array<ids_type,dim> cellPadDim,
	transform_type pointTransform,
	size_t n_part,
	size_t start,
	vector_pos_type vPos,
	vector_cnt_type numPartInCell,
	vector_pids_type cellIndex_LocalIndex)
{
	unsigned int i, cid, ins;
	ids_type e[dim+1];

	i = threadIdx.x + blockIdx.x * blockDim.x + start;
	ins = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= n_part) return;

	pos_type p[dim];

	for (size_t k = 0 ; k < dim ; k++)
		p[k] = vPos.template get<0>(i)[k];

	cid = cid_<dim,ids_type,transform_type>::get_cid(numCellDiv,unitCellP2,cellPadDim,pointTransform,p,e);

	e[dim] = atomicAdd(&numPartInCell.template get<0>(cid), 1);
	cellIndex_LocalIndex.template get<0>(ins)[0] = cid;
	cellIndex_LocalIndex.template get<0>(ins)[1] = e[dim];
}

template<unsigned int dim, typename pos_type, typename ids_type, typename transform_type,
typename vector_pos_type, typename vector_cnt_type>
__global__ void fill_cellIndex(
	openfpm::array<ids_type,dim> numCellDiv,
	openfpm::array<pos_type,dim> unitCellP2,
	openfpm::array<ids_type,dim> cellPadDim,
	transform_type pointTransform,
	size_t n_part,
	size_t start,
	vector_pos_type vPos,
	vector_cnt_type cellIndex)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x + start;
	int ins = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= n_part) return;

	Point<dim,pos_type> p;

	for (size_t k = 0 ; k < dim ; k++)
		p[k] = vPos.template get<0>(i)[k];

	cellIndex.template get<0>(ins) = cid_<dim,ids_type,transform_type>::get_cid(numCellDiv,unitCellP2,cellPadDim,pointTransform,p);
}

template<typename vector_sparse, typename vector_cell>
__global__ void fill_vsCellIndex_PartIndex(
	vector_sparse vecSparseCellIndex_PartIndex,
	vector_cell cellIndex)
{
	vecSparseCellIndex_PartIndex.init();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	if (p < cellIndex.size())
	{
		int c = cellIndex.template get<0>(p);
		vecSparseCellIndex_PartIndex.template insert<0>(c) = p;
	}

	vecSparseCellIndex_PartIndex.flush_block_insert();
}

template<typename vector_starts_type, typename vector_pids_type, typename vector_cells_type>
__global__ void fill_cells(
	vector_starts_type numPartInCellPrefixSum,
	vector_pids_type cellIndex_LocalIndex,
	vector_cells_type cellIndexLocalIndexToUnsorted,
	size_t startParticle=0)
{
	unsigned int cid, id, cellStart;

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= cellIndex_LocalIndex.size()) return;

	cid = cellIndex_LocalIndex.template get<0>(tid)[0];

	cellStart = numPartInCellPrefixSum.template get<0>(cid);
	id = cellStart + cellIndex_LocalIndex.template get<0>(tid)[1];

	cellIndexLocalIndexToUnsorted.template get<0>(id) = tid + startParticle;
}

template <typename vector_map_type, typename vector_cells_type>
__global__ void constructSortUnsortBidirectMap(
	vector_map_type sortedToUnsortedIndex,
	vector_map_type unsortedToSortedIndex,
	const vector_cells_type cellIndexLocalIndexToUnsorted)
{
	unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= sortedToUnsortedIndex.size())	{return;}

	unsigned int pid = cellIndexLocalIndexToUnsorted.template get<0>(tid);

	sortedToUnsortedIndex.template get<0>(tid) = pid;
	unsortedToSortedIndex.template get<0>(pid) = tid;
}

template <typename vector_type, typename vector_map_type>
__global__ void reorderParticlesPos(
	const vector_type vectorIn,
	vector_type vectorOut,
	const vector_map_type indexMap,
	size_t start = 0)
{
	int keyIn = start + threadIdx.x + blockIdx.x * blockDim.x;
	if (keyIn >= indexMap.size())	{return;}

	unsigned int keyOut = indexMap.template get<0>(keyIn);

	vectorOut.set(keyOut,vectorIn,keyIn);
}

template <typename vector_type, typename vector_map_type, unsigned int ... prp>
__global__ void reorderParticlesPrp(
	const vector_type vectorIn,
	vector_type vectorOut,
	vector_map_type indexMap,
	size_t start = 0)
{
	int keyIn = start + threadIdx.x + blockIdx.x * blockDim.x;
	if (keyIn >= indexMap.size())	{return;}

	unsigned int keyOut = indexMap.template get<0>(keyIn);

	vectorOut.template set<prp ...>(keyOut,vectorIn,keyIn);
}

template <typename vector_type, typename vector_map_type>
__global__ void reorderParticlesPosCoalWrite(
	const vector_type vectorIn,
	vector_type vectorOut,
	const vector_map_type indexMap,
	size_t start = 0)
{
	int keyWrite = start + threadIdx.x + blockIdx.x * blockDim.x;
	if (keyWrite >= indexMap.size())	{return;}

	unsigned int keyRead = indexMap.template get<0>(keyWrite);

	vectorOut.set(keyWrite,vectorIn,keyRead);
}

template <typename vector_type, typename vector_map_type, unsigned int ... prp>
__global__ void reorderParticlesPrpCoalWrite(
	const vector_type vectorIn,
	vector_type vectorOut,
	vector_map_type indexMap,
	size_t start = 0)
{
	int keyWrite = start + threadIdx.x + blockIdx.x * blockDim.x;
	if (keyWrite >= indexMap.size())	{return;}

	unsigned int keyRead = indexMap.template get<0>(keyWrite);

	vectorOut.template set<prp ...>(keyWrite,vectorIn,keyRead);
}


template<typename vector_sort_index, typename vector_out_type>
__global__ void mark_domain_particles(
	vector_sort_index sortedToUnsortedIndex,
	vector_out_type isSortedDomainOrGhost,
	size_t ghostMarker)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i >= sortedToUnsortedIndex.size()) return;

	isSortedDomainOrGhost.template get<0>(i) = (sortedToUnsortedIndex.template get<0>(i) < ghostMarker)?1:0;
}

template<typename scan_type, typename vector_out_type>
__global__ void collect_domain_ghost_ids(
	scan_type isUnsortedDomainOrGhostPrefixSum,
	vector_out_type sortedToSortedIndexNoGhost)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i >= isUnsortedDomainOrGhostPrefixSum.size()-1) return;

	auto pp = isUnsortedDomainOrGhostPrefixSum.template get<0>(i+1);
	auto p = isUnsortedDomainOrGhostPrefixSum.template get<0>(i);

	if (pp != p)
		sortedToSortedIndexNoGhost.template get<0>(isUnsortedDomainOrGhostPrefixSum.template get<0>(i)) = i;
}

template<typename cl_sparse_type, typename vector_type, typename vector_type2>
__global__ void countNonEmptyNeighborCells(
	cl_sparse_type vecSparseCellIndex_PartIndex,
	vector_type neighborCellCount,
	vector_type2 neighborCellOffset)
{
	typedef typename cl_sparse_type::index_type index_type;

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= vecSparseCellIndex_PartIndex.size()) {return;}

	openfpm::sparse_index<index_type> id;
	id.id = tid;

	index_type cell = vecSparseCellIndex_PartIndex.get_index(id);

	for (int i = 0 ; i < neighborCellOffset.size() ; i++)
	{
		index_type neighborCellIndex = cell + neighborCellOffset.template get<0>(i);
		index_type start = vecSparseCellIndex_PartIndex.template get<0>(neighborCellIndex);

		if (start != (index_type)-1)
		{
			// Cell exist
			neighborCellCount.template get<0>(tid) += 1;
		}
	}
};

template<typename cl_sparse_type, typename vector_type, typename vector_type2, typename vector_type3>
__global__ void fillNeighborCellList(
	cl_sparse_type vecSparseCellIndex_PartIndex,
	vector_type neighborCellCountPrefixSum,
	vector_type2 neighborCellOffset,
	vector_type3 neighborPartIndexFrom_To,
	typename cl_sparse_type::index_type stop)
{
	typedef typename cl_sparse_type::index_type index_type;

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= vecSparseCellIndex_PartIndex.size())	{return;}

	openfpm::sparse_index<index_type> id; id.id = tid;
	index_type cell = vecSparseCellIndex_PartIndex.get_index(id);

	for (int i = 0, cnt = 0; i < neighborCellOffset.size(); i++)
	{
		index_type neighborCellIndex = cell + neighborCellOffset.template get<0>(i);
		auto sid = vecSparseCellIndex_PartIndex.get_sparse(neighborCellIndex);

		if (sid.id != vecSparseCellIndex_PartIndex.size())
		{
			neighborPartIndexFrom_To.template get<0>(neighborCellCountPrefixSum.template get<0>(tid) + cnt) = vecSparseCellIndex_PartIndex.template get<0>(sid);

			if (++sid.id != vecSparseCellIndex_PartIndex.size())
				neighborPartIndexFrom_To.template get<1>(neighborCellCountPrefixSum.template get<0>(tid) + cnt) = vecSparseCellIndex_PartIndex.template get<0>(sid);
			else
				neighborPartIndexFrom_To.template get<1>(neighborCellCountPrefixSum.template get<0>(tid) + cnt) = stop;

			++cnt;
		}
	}
};


/////////////////////////// THIS ONLY WORK IF NVCC IS COMPILING THIS //////////////////////////////////////////////////////

#endif

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CUDA_CUDA_CELL_LIST_UTIL_FUNC_HPP_ */
