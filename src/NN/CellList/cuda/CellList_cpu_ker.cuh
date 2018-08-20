/*
 * CellList_gpu_ker.cuh
 *
 *  Created on: Jul 30, 2018
 *      Author: i-bird
 */

#ifndef CELLLIST_CPU_KER_CUH_
#define CELLLIST_CPU_KER_CUH_

#include "Cuda_cell_list_util_func.hpp"

template<unsigned int dim, typename T, typename Mem_type, typename transform>
class CellList_cpu_ker: Mem_type
{
	typedef typename Mem_type::local_index_type cnt_type;

	typedef typename Mem_type::local_index_type ids_type;

	//! Spacing
	openfpm::array<T,dim,cnt_type> spacing_c;

	//! \brief number of sub-divisions in each direction
	openfpm::array<ids_type,dim,cnt_type> div_c;

	//! \brief cell padding
	openfpm::array<ids_type,dim,cnt_type> off;

	//! transformation
	transform t;

public:

	CellList_cpu_ker(const Mem_type & mt,
			 	 	 openfpm::array<T,dim,cnt_type> & spacing_c,
			 	 	 openfpm::array<ids_type,dim,cnt_type> & div_c,
			 	 	 openfpm::array<ids_type,dim,cnt_type> & off,
			 	 	 const transform & t)
	:Mem_type(mt),spacing_c(spacing_c),div_c(div_c),off(off),t(t)
	{}

	inline __device__ unsigned int getCell(const Point<dim,T> & xp)
	{
		return cid_<dim,cnt_type,ids_type,transform>::get_cid(div_c,spacing_c,off,t,xp);
	}

	/*! \brief Return the number of elements in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return number of elements in the cell
	 *
	 */
	inline __device__ int getNelements(unsigned int cell) const
	{
		return Mem_type::getNelements(cell);
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
	inline __device__  unsigned int get(unsigned int cell, unsigned int ele)
	{
		return Mem_type::get(cell,ele);
	}

};


#endif /* CELLLIST_GPU_KER_CUH_ */
