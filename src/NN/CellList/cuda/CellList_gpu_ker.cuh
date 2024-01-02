/*
 * CellList_gpu_ker.cuh
 *
 *  Created on: Jul 30, 2018
 *      Author: i-bird
 */

#ifndef CELLLIST_GPU_KER_CUH_
#define CELLLIST_GPU_KER_CUH_

#include "NN/CellList/CellList_def.hpp"
#include "NN/CellList/cuda/CellDecomposer_gpu_ker.cuh"

// #ifdef USE_LOW_REGISTER_ITERATOR

// #ifdef __NVCC__
// __constant__ int cells_striding[126];
// #endif

// pr_int=openfpm::math::pow(2*radius_int+1,dim)
// template<unsigned int dim, int radius_int, unsigned int pr_int, typename ids_type>
// struct NN_gpu_int_base_lr_impl
// {
// 	unsigned int ca_pnt;
// 	unsigned int ca_lin;
// 	unsigned int c_id;

// 	__device__ inline void init(const grid_key_dx<dim,ids_type> & cellPosition, const openfpm::array<ids_type,dim> & numCellDim)
// 	{
// #ifdef __NVCC__
// 		ca_pnt = 0;

// 		ca_lin = cid_<dim,ids_type,int>::get_cid(numCellDim,cellPosition);
// 		c_id = ca_lin + cells_striding[ca_pnt];
// #endif
// 	}

// 	__device__ inline void nextCell(const openfpm::array<ids_type,dim> & numCellDim)
// 	{
// #ifdef __NVCC__
// 		++ca_pnt;

// 		c_id = ca_lin + cells_striding[ca_pnt];
// #endif
// 	}

// 	__device__ inline bool isNextCell()
// 	{
// 		return ca_pnt < pr_int;
// 	}
// };

// #endif

template<unsigned int dim, int radius_int, typename ids_type>
struct NN_gpu_int_base_hr_impl
{
	grid_key_dx<dim,ids_type> cellPosAct;
	grid_key_dx<dim,ids_type> cellPosStart;
	grid_key_dx<dim,ids_type> cellPosStop;

	unsigned int cellIndexAct;

	__device__ __host__ inline void init(
		const grid_key_dx<dim,ids_type> & cellPosition,
		const openfpm::array<ids_type,dim> & numCellDim)
	{
		for (unsigned int i = 0 ; i < dim ; i++)
		{
			cellPosStart.set_d(i,cellPosition.get(i) - radius_int);
			cellPosStop.set_d(i,cellPosition.get(i) + radius_int);
			cellPosAct.set_d(i,cellPosition.get(i) - radius_int);
		}

		cellIndexAct = cid_<dim,ids_type,int>::get_cid(numCellDim,cellPosStart);
	}

	__device__ __host__ inline void nextCell(const openfpm::array<ids_type,dim> & numCellDim)
	{
		cellPosAct.set_d(0,cellPosAct.get(0)+1);

		//! check the overflow of all the index with exception of the last dimensionality

		for (int i = 0 ; i < dim-1 ; i++)
		{
			unsigned int id = cellPosAct.get(i);
			if ((int)id > cellPosStop.get(i))
			{
				// ! overflow, increment the next index

				cellPosAct.set_d(i,cellPosStart.get(i));
				id = cellPosAct.get(i+1);
				cellPosAct.set_d(i+1,id+1);
			}
			else break;
		}

		cellIndexAct = cid_<dim,ids_type,int>::get_cid(numCellDim,cellPosAct);
	}

	__device__ __host__ inline bool isNextCell()
	{
		return cellPosAct.get(dim-1) <= cellPosStop.get(dim-1);
	}
};

// #ifdef USE_LOW_REGISTER_ITERATOR

// template<unsigned int dim, int radius_int, unsigned int pr_int, typename ids_type>
// struct NN_gpu_int_base: public NN_gpu_int_base_hr_impl<dim,radius_int,pr_int,ids_type>
// {};

// template<int radius_int, unsigned int pr_int, typename ids_type>
// struct NN_gpu_int_base<2,radius_int,pr_int,ids_type>: public NN_gpu_int_base_lr_impl<2,radius_int,pr_int,ids_type>
// {};

// template<int radius_int, unsigned int pr_int, typename ids_type>
// struct NN_gpu_int_base<3,radius_int,pr_int,ids_type>: public NN_gpu_int_base_lr_impl<3,radius_int,pr_int,ids_type>
// {};

// #else

template<unsigned int dim, int radius_int, typename ids_type>
struct NN_gpu_int_base: public NN_gpu_int_base_hr_impl<dim,radius_int,ids_type>
{};

// #endif

template<unsigned int dim, typename ids_type, unsigned int radius_int, bool is_sparse>
class NN_gpu_it: public NN_gpu_int_base<dim,radius_int,ids_type>
{
	const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & numPartInCellPrefixSum;
	const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & sortedToUnsortedIndex;
	const openfpm::array<ids_type,dim> & numCellDim;
	const openfpm::array<ids_type,dim> & cellPadDim;

	unsigned int neighborPartIndexStart;
	unsigned int neighborPartIndexStop;

	inline __device__ __host__ void SelectValid()
	{
		while (neighborPartIndexStart == neighborPartIndexStop && isNext())
		{
			this->nextCell(numCellDim);

			if (isNext() == false) {break;}

			neighborPartIndexStart = numPartInCellPrefixSum.template get<0>(this->cellIndexAct);
			neighborPartIndexStop = numPartInCellPrefixSum.template get<0>(this->cellIndexAct+1);
		}
	}


public:

	inline __device__ __host__ NN_gpu_it(
		const grid_key_dx<dim,ids_type> & cellPosition,
		const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & numPartInCellPrefixSum,
		const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & sortedToUnsortedIndex,
		const openfpm::array<ids_type,dim> & numCellDim,
		const openfpm::array<ids_type,dim> & cellPadDim)
	: numPartInCellPrefixSum(numPartInCellPrefixSum),
	sortedToUnsortedIndex(sortedToUnsortedIndex),
	numCellDim(numCellDim),
	cellPadDim(cellPadDim)
	{
		this->init(cellPosition,numCellDim);

		neighborPartIndexStart = numPartInCellPrefixSum.template get<0>(this->cellIndexAct);
		neighborPartIndexStop = numPartInCellPrefixSum.template get<0>(this->cellIndexAct+1);

		SelectValid();
	}

	inline __device__ __host__ unsigned int get_sort()
	{
		return neighborPartIndexStart;
	}

	inline __device__ __host__ unsigned int get()
	{
		return sortedToUnsortedIndex.template get<0>(neighborPartIndexStart);
	}

	inline __device__ __host__ NN_gpu_it<dim,ids_type,radius_int,is_sparse> & operator++()
	{
		++neighborPartIndexStart; SelectValid();

		return *this;
	}

	inline __device__ unsigned int get_start(unsigned int ce_id)
	{
		return numPartInCellPrefixSum.template get<0>(ce_id);
	}

	inline __device__ unsigned int get_cid()
	{
		return this->cellIndexAct;
	}

	inline __device__ __host__ bool isNext()
	{
		return this->isNextCell();
	}
};

template<unsigned int dim, typename ids_type, unsigned int radius_int>
class NN_gpu_it<dim,ids_type,radius_int,true>
{
	unsigned int neighborPartIndexStart;
	unsigned int neighborPartIndexStop;
	unsigned int neighborCellIndexStart;
	unsigned int neighborCellIndexStop;

	const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & sortedToUnsortedIndex;
	const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & neighborCellCountPrefixSum;
	const openfpm::vector_gpu_ker<aggregate<unsigned int,unsigned int>,memory_traits_inte> & neighborPartIndexStart_To;

	__device__ __host__ void SelectValid()
	{
		while ((neighborPartIndexStart == neighborPartIndexStop) && isNext())
		{
			++neighborCellIndexStart;

			if (neighborCellIndexStart < neighborCellIndexStop)
			{
				neighborPartIndexStart = neighborPartIndexStart_To.template get<0>(neighborCellIndexStart);
				neighborPartIndexStop = neighborPartIndexStart_To.template get<1>(neighborCellIndexStart);
			}
		}
	}

public:

	__device__ NN_gpu_it(
		unsigned int cellIndex,
		const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & neighborCellCountPrefixSum,
		const openfpm::vector_gpu_ker<aggregate<unsigned int,unsigned int>,memory_traits_inte> & neighborPartIndexStart_To,
		const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & sortedToUnsortedIndex)
	: sortedToUnsortedIndex(sortedToUnsortedIndex),
	neighborCellCountPrefixSum(neighborCellCountPrefixSum),
	neighborPartIndexStart_To(neighborPartIndexStart_To)
	{
		if (cellIndex == (unsigned int)-1)
		{
			neighborCellIndexStop = neighborCellIndexStart;
			return;
		}

		neighborCellIndexStart = neighborCellCountPrefixSum.template get<0>(cellIndex);
		neighborCellIndexStop = neighborCellCountPrefixSum.template get<0>(cellIndex + 1);

		neighborPartIndexStart = neighborPartIndexStart_To.template get<0>(neighborCellIndexStart);
		neighborPartIndexStop = neighborPartIndexStart_To.template get<1>(neighborCellIndexStart);

		SelectValid();
	}

	__device__ unsigned int get_sort()
	{
		return neighborPartIndexStart;
	}

	__device__ unsigned int get()
	{
		return sortedToUnsortedIndex.template get<0>(neighborPartIndexStart);
	}

	__device__ __host__ NN_gpu_it<dim,ids_type,radius_int,true> & operator++()
	{
		++neighborPartIndexStart; SelectValid();

		return *this;
	}

	inline __device__ __host__ bool isNext()
	{
		return neighborCellIndexStart < neighborCellIndexStop;
	}
};


template<unsigned int dim, typename ids_type>
class NN_gpu_it_radius
{
	unsigned int cellIndexAct;
	unsigned int radNeighborCellIndex_i;
	unsigned int neighborPartIndexStart;
	unsigned int neighborCellIndex;

	const openfpm::vector_gpu_ker<aggregate<int>,memory_traits_inte> & radNeighborCellIndex;
	const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & neighborCellCountPrefixSum;
	const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & sortedToUnsortedIndex;
	const openfpm::array<ids_type,dim> & numCellDim;
	const openfpm::array<ids_type,dim> & cellPadDim;

	__device__ __host__ inline void SelectValid()
	{
		// here it's not neighborCellCountPrefixSum.template get<0>(cellIndexAct + radNeighborCellIndex.template get<0>(radNeighborCellIndex_i)+1)
		// as we check whether all the particles in cell neighborCellIndex were iterated, until the next cell neighborCellIndex+1
		while (isNext() && neighborPartIndexStart == neighborCellCountPrefixSum.template get<0>(neighborCellIndex+1))
		{
			radNeighborCellIndex_i++;

			if (radNeighborCellIndex_i >= radNeighborCellIndex.size())
				break;

			neighborCellIndex = cellIndexAct + radNeighborCellIndex.template get<0>(radNeighborCellIndex_i);
			neighborPartIndexStart = neighborCellCountPrefixSum.template get<0>(neighborCellIndex);
		}
	}

public:

	__device__ __host__ inline NN_gpu_it_radius(
		const grid_key_dx<dim,ids_type> & cellPosition,
		const openfpm::vector_gpu_ker<aggregate<int>,memory_traits_inte> & radNeighborCellIndex,
		const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & neighborCellCountPrefixSum,
		const openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & sortedToUnsortedIndex,
		const openfpm::array<ids_type,dim> & numCellDim,
		const openfpm::array<ids_type,dim> & cellPadDim)
	: radNeighborCellIndex_i(0),
	radNeighborCellIndex(radNeighborCellIndex),
	neighborCellCountPrefixSum(neighborCellCountPrefixSum),
	sortedToUnsortedIndex(sortedToUnsortedIndex),
	numCellDim(numCellDim),
	cellPadDim(cellPadDim)
	{
		cellIndexAct = cid_<dim,ids_type,int>::get_cid(numCellDim,cellPosition);
		neighborCellIndex = cellIndexAct + radNeighborCellIndex.template get<0>(radNeighborCellIndex_i);
		neighborPartIndexStart = neighborCellCountPrefixSum.template get<0>(neighborCellIndex);
		SelectValid();
	}

	__device__ unsigned int get_sort()
	{
		return neighborPartIndexStart;
	}

	__device__ unsigned int get()
	{
		return sortedToUnsortedIndex.template get<0>(neighborPartIndexStart);
	}

	__device__ __host__ NN_gpu_it_radius<dim,ids_type> & operator++()
	{
		++neighborPartIndexStart; SelectValid();

		return *this;
	}

	__device__ unsigned int get_start(unsigned int ce_id)
	{
		return neighborCellCountPrefixSum.template get<0>(ce_id);
	}

	__device__ unsigned int get_cid()
	{
		return neighborCellIndex;
	}

	__device__ __host__ bool isNext()
	{
		return radNeighborCellIndex_i < radNeighborCellIndex.size();
	}
};

template<unsigned int dim, typename T, typename ids_type, typename transform_type, bool is_sparse>
class CellList_gpu_ker: public CellDecomposer_gpu_ker<dim,T,ids_type,transform_type>
{
	//! starting point for each cell
	openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> numPartInCellPrefixSum;

	//! Sorted to non sorted ids conversion
	openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> sortedToUnsortedIndex;

	//! Domain particles ids
	openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> indexSorted;

	//! radius cells
	openfpm::vector_gpu_ker<aggregate<int>,memory_traits_inte> radNeighborCellIndex;

	//! Ghost particle marker
	unsigned int ghostMarker;

public:

	typedef int yes_is_gpu_ker_celllist;

	//! Indicate this structure has a function to check the device pointer
	typedef int yes_has_check_device_pointer;

	__host__ __device__ inline CellList_gpu_ker()
	:ghostMarker(0)
	{}

	__host__ __device__ inline CellList_gpu_ker(
		openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> numPartInCellPrefixSum,
		openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> sortedToUnsortedIndex,
		openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> indexSorted,
		openfpm::vector_gpu_ker<aggregate<int>,memory_traits_inte> radNeighborCellIndex,
		openfpm::array<T,dim> unitCellP2,
		openfpm::array<ids_type,dim> & numCellDim,
		openfpm::array<ids_type,dim> & cellPadDim,
		const transform_type & pointTransform,
		unsigned int ghostMarker,
		SpaceBox<dim,T> cellListSpaceBox,
		grid_sm<dim,void> cellListGrid,
		Point<dim,long int> cellShift)
	: CellDecomposer_gpu_ker<dim,T,ids_type,transform_type>(unitCellP2,numCellDim,cellPadDim,pointTransform,cellListSpaceBox,cellListGrid,cellShift),
	numPartInCellPrefixSum(numPartInCellPrefixSum),
	sortedToUnsortedIndex(sortedToUnsortedIndex),
	indexSorted(indexSorted),
	radNeighborCellIndex(radNeighborCellIndex),
	ghostMarker(ghostMarker)
	{}


	template<unsigned int stub = NO_CHECK>
	inline __device__ __host__ NN_gpu_it<dim,ids_type,1,is_sparse> getNNIterator(
		const grid_key_dx<dim,ids_type> & cellPosition)
	{
		NN_gpu_it<dim,ids_type,1,is_sparse> ngi(cellPosition,numPartInCellPrefixSum,sortedToUnsortedIndex,this->get_div_c(),this->get_off());

		return ngi;
	}

	inline __device__ __host__ NN_gpu_it_radius<dim,ids_type> getNNIteratorRadius(
		const grid_key_dx<dim,ids_type> & cellPosition)
	{
		NN_gpu_it_radius<dim,ids_type> ngi(cellPosition,radNeighborCellIndex,numPartInCellPrefixSum,sortedToUnsortedIndex,this->get_div_c(),this->get_off());

		return ngi;
	}

	template<unsigned int radius_int = 2> 
	inline __device__ NN_gpu_it<dim,ids_type,radius_int,is_sparse> getNNIteratorBox(
		const grid_key_dx<dim,ids_type> & cellPosition)
	{
		NN_gpu_it<dim,ids_type,radius_int,is_sparse> ngi(cellPosition,numPartInCellPrefixSum,sortedToUnsortedIndex,this->get_div_c(),this->get_off());

		return ngi;
	}

	inline __device__ openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & getDomainSortIds()
	{
		return indexSorted;
	}

	inline __device__ openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & getSortToNonSort()
	{
		return sortedToUnsortedIndex;
	}

	/*! \brief Get the number of cells this cell-list contain
	 *
	 * \return number of cells
	 */
	inline __device__ unsigned int getNCells() const
	{
		return numPartInCellPrefixSum.size() - 1;
	}

	/*! \brief Return the number of elements in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return number of elements in the cell
	 *
	 */
	inline __device__ unsigned int getNelements(unsigned int cell_id) const
	{
		return numPartInCellPrefixSum.template get<0>(cell_id+1) - numPartInCellPrefixSum.template get<0>(cell_id);
	}

	/*! \brief Get an element in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline __device__ unsigned int get(
		unsigned int cell, 
		unsigned int ele)
	{
		unsigned int p_id = numPartInCellPrefixSum.template get<0>(cell) + ele;
		return sortedToUnsortedIndex.template get<0>(p_id);
	}


	inline __device__ unsigned int get_g_m()
	{
		return ghostMarker;
	}

#ifdef SE_CLASS1

	/*! \brief Check if the device pointer is owned by this structure
	 *
	 * \return a structure pointer check with information about the match
	 *
	 */
	pointer_check check_device_pointer(void * ptr)
	{
		pointer_check pc;

		pc = numPartInCellPrefixSum.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Cell index overflow (numPartInCellPrefixSum): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = sortedToUnsortedIndex.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (str): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = indexSorted.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (indexSorted): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = radNeighborCellIndex.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (indexSorted): ") + "\n" + pc.match_str;
			return pc;
		}

		return pc;
	}

#endif
};


template<unsigned int dim, typename T, typename ids_type, typename transform_type>
class CellList_gpu_ker<dim,T,ids_type,transform_type,true>: public CellDecomposer_gpu_ker<dim,T,ids_type,transform_type>
{
	//! starting point for each cell
	openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> neighborCellCountPrefixSum;

	//! starting point for each cell
	openfpm::vector_gpu_ker<aggregate<unsigned int, unsigned int>,memory_traits_inte> neighborCellVsIndexFrom_To;

	//! Sorted to non sorted ids conversion
	openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> sortedToUnsortedIndex;

	//! Domain particles ids
	openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> indexSorted;

	//! Set of cells sparse
	openfpm::vector_sparse_gpu_ker<aggregate<unsigned int>,int,memory_traits_inte> vecSparseCellIndex_PartIndex;

	//! Ghost particle marker
	unsigned int ghostMarker;

public:

	//! Indicate this structure has a function to check the device pointer
	typedef int yes_has_check_device_pointer;

	__host__ __device__ inline CellList_gpu_ker(
		openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> neighborCellCountPrefixSum,
		openfpm::vector_gpu_ker<aggregate<unsigned int, unsigned int>,memory_traits_inte> neighborCellVsIndexFrom_To,
		openfpm::vector_sparse_gpu_ker<aggregate<unsigned int>,int,memory_traits_inte> vecSparseCellIndex_PartIndex,
		openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> sortedToUnsortedIndex,
		openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> indexSorted,
		openfpm::array<T,dim> & unitCellP2,
		openfpm::array<ids_type,dim> & numCellDim,
		openfpm::array<ids_type,dim> & cellPadDim,
		const transform_type & pointTransform,
		unsigned int ghostMarker,
		SpaceBox<dim,T> cellListSpaceBox,
		grid_sm<dim,void> cellListGrid,
		Point<dim,long int> cellShift
	)
	: CellDecomposer_gpu_ker<dim,T,ids_type,transform_type>(unitCellP2,numCellDim,cellPadDim,pointTransform,cellListSpaceBox,cellListGrid,cellShift),
	neighborCellCountPrefixSum(neighborCellCountPrefixSum),
	neighborCellVsIndexFrom_To(neighborCellVsIndexFrom_To),
	sortedToUnsortedIndex(sortedToUnsortedIndex),
	indexSorted(indexSorted),
	vecSparseCellIndex_PartIndex(vecSparseCellIndex_PartIndex),
	ghostMarker(ghostMarker)
	{}

	inline __device__ auto getCell(
		const Point<dim,T> & xp) 
	const -> decltype(vecSparseCellIndex_PartIndex.get_sparse(0))
	{
		unsigned int cell = cid_<dim,ids_type,transform_type>::get_cid(this->get_div_c(),this->get_spacing_c(),this->get_off(),this->get_t(),xp);

		return vecSparseCellIndex_PartIndex.get_sparse(cell);
	}


	template<unsigned int stub = NO_CHECK>
	inline __device__ NN_gpu_it<dim,ids_type,1,true> getNNIterator(
		decltype(vecSparseCellIndex_PartIndex.get_sparse(0)) cId)
	{
		NN_gpu_it<dim,ids_type,1,true> ngi(cId.id,neighborCellCountPrefixSum,neighborCellVsIndexFrom_To,sortedToUnsortedIndex);

		return ngi;
	}

	template<unsigned int radius_int = 2>
	inline __device__ NN_gpu_it<dim,ids_type,radius_int,true> getNNIteratorBox(
		decltype(vecSparseCellIndex_PartIndex.get_sparse(0)) cId)
	{
		NN_gpu_it<dim,ids_type,radius_int,true> ngi(cId.id,neighborCellCountPrefixSum,neighborCellVsIndexFrom_To,sortedToUnsortedIndex);

		return ngi;
	}


	inline __device__ openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & getDomainSortIds()
	{
		return indexSorted;
	}

	inline __device__ openfpm::vector_gpu_ker<aggregate<unsigned int>,memory_traits_inte> & getSortToNonSort()
	{
		return sortedToUnsortedIndex;
	}


	inline __device__ unsigned int get_g_m()
	{
		return ghostMarker;
	}

#ifdef SE_CLASS1

	/*! \brief Check if the device pointer is owned by this structure
	 *
	 * \return a structure pointer check with information about the match
	 *
	 */
	pointer_check check_device_pointer(void * ptr)
	{
		pointer_check pc;

		pc = neighborCellCountPrefixSum.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Cell index overflow (starts): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = neighborCellVsIndexFrom_To.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Cell particle buffer overflow (neighborCellVsIndexFrom_To): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = sortedToUnsortedIndex.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (str): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = indexSorted.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (indexSorted): ") + "\n" + pc.match_str;
			return pc;
		}

		return pc;
	}

#endif
};


#endif /* CELLLIST_GPU_KER_CUH_ */
