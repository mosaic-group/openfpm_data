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

	//! Grid structure of the Cell list
	grid_sm<dim,void> gr_cell;

	//! cell_shift
	Point<dim,long int> cell_shift;

	//! Unit box of the Cell list
	Box<dim,T> box_unit;

public:

	CellList_cpu_ker(const Mem_type & mt,
			 	 	 openfpm::array<T,dim,cnt_type> & spacing_c,
			 	 	 openfpm::array<ids_type,dim,cnt_type> & div_c,
			 	 	 openfpm::array<ids_type,dim,cnt_type> & off,
					 grid_sm<dim,void> & gr_cell,
					 Point<dim,long int> & cell_shift,
					 Box<dim,T> & box_unit,
			 	 	 const transform & t)
	:Mem_type(mt),spacing_c(spacing_c),div_c(div_c),off(off),t(t),gr_cell(gr_cell),cell_shift(cell_shift),box_unit(box_unit)
	{}

	inline __device__ unsigned int getCell(const Point<dim,T> & xp) const
	{
		return cid_<dim,cnt_type,ids_type,transform>::get_cid(div_c,spacing_c,off,t,xp);
	}

	/*! \brief Return the underlying grid information of the cell list
	 *
	 * \return the grid infos
	 *
	 */
	const __device__ grid_sm<dim,void> & getGrid() const
	{
		return gr_cell;
	}

		/*! \brief Get the cell-ids
	 *
	 * Convert the point coordinates into the cell ids (Careful it include padding)
	 *
	 * \param pos Point position
	 *
	 * \return the cell-ids ad a grid_key_dx<dim>
	 *
	 */
	inline __device__ grid_key_dx<dim> getCellGrid(const Point<dim,T> & pos) const
	{
		grid_key_dx<dim> key;
		key.set_d(0,ConvertToID(pos,0));

		for (size_t s = 1 ; s < dim ; s++)
		{
			key.set_d(s,ConvertToID(pos,s));
		}

		return key;
	}

	/*! \brief Get the cell-ids
	 *
	 * Convert the point coordinates into the cell ids
	 *
	 * \param pos Point position
	 *
	 * \return the cell-ids ad a grid_key_dx<dim>
	 *
	 */
	inline __device__ grid_key_dx<dim> getCellGrid(const T (& pos)[dim]) const
	{
		grid_key_dx<dim> key;
		key.set_d(0,ConvertToID(pos,0));

		for (size_t s = 1 ; s < dim ; s++)
		{
			key.set_d(s,ConvertToID(pos,s));
		}

		return key;
	}

	/*! \brief Convert the coordinates into id
	 *
	 * \param x coordinate
	 * \param s dimension
	 *
	 */
	inline __device__ size_t ConvertToID(const T (&x)[dim], size_t s) const
	{
		size_t id = openfpm::math::size_t_floor(t.transform(x,s) / box_unit.getHigh(s)) + off[s];
		id = (id >= gr_cell.size(s))?(gr_cell.size(s)-1-cell_shift.get(s)):id-cell_shift.get(s);
		return id;
	}

	/*! \brief Convert the coordinates into id
	 *
	 * \param x point
	 * \param s dimension
	 *
	 */
	inline __device__ size_t ConvertToID(const Point<dim,T> & x, size_t s, size_t sc = 0) const
	{
		size_t id = openfpm::math::size_t_floor(t.transform(x,s) / box_unit.getHigh(s)) + off[s];
		id = (id >= gr_cell.size(s))?(gr_cell.size(s)-1-cell_shift.get(s)):id-cell_shift.get(s);
		return id;
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
