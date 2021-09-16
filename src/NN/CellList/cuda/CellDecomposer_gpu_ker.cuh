/*
 * CellDecomposer_gpu_ker.hpp
 *
 *  Created on: Apr 28, 2019
 *      Author: i-bird
 */

#ifndef CELLDECOMPOSER_GPU_KER_HPP_
#define CELLDECOMPOSER_GPU_KER_HPP_

#include "util/cuda_launch.hpp"
#include "util/multi_array_openfpm/array_openfpm.hpp"
#include "Grid/grid_sm.hpp"
#include "NN/CellList/cuda/Cuda_cell_list_util_func.hpp"
#include "NN/CellList/CellDecomposer.hpp"

template <unsigned int dim, typename T, typename cnt_type, typename ids_type, typename transform>
class CellDecomposer_gpu_ker
{
	//! Spacing
	openfpm::array<T,dim,cnt_type> spacing_c;

	//! \brief number of sub-divisions in each direction
	openfpm::array<ids_type,dim,cnt_type> div_c;

	//! \brief cell offset
	openfpm::array<ids_type,dim,cnt_type> off;

	//! transformation
	transform t;

	//! Unit box of the Cell list
	SpaceBox<dim,T> box_unit;

	//! Grid structure of the Cell list
	grid_sm<dim,void> gr_cell;

	//! cell_shift
	Point<dim,long int> cell_shift;
public:

	__device__ __host__ CellDecomposer_gpu_ker()
	{}

	__device__ __host__ CellDecomposer_gpu_ker(
		openfpm::array<T,dim,cnt_type> & spacing_c,
		openfpm::array<ids_type,dim,cnt_type> & div_c,
		openfpm::array<ids_type,dim,cnt_type> & off,
		const transform & t)
	: spacing_c(spacing_c),div_c(div_c),off(off),t(t)
	{}

	__device__ __host__ CellDecomposer_gpu_ker(
		openfpm::array<T,dim,cnt_type> & spacing_c,
		openfpm::array<ids_type,dim,cnt_type> & div_c,
		openfpm::array<ids_type,dim,cnt_type> & off,
		const transform & t,
		SpaceBox<dim,T> box_unit,
		grid_sm<dim,void> gr_cell,
		Point<dim,long int> cell_shift)
	: spacing_c(spacing_c),div_c(div_c),off(off),t(t),
	box_unit(box_unit), gr_cell(gr_cell), cell_shift(cell_shift)
	{}

	__device__ __host__ grid_sm<dim,void> getGrid()
	{
		size_t sz[dim];

		for (size_t i = 0 ; i < dim ; i++)
		{
			sz[i] = div_c[i] + 2*off[i];
		}

		return grid_sm<dim,void> (sz);
	}

	__device__ __host__ void getGridSize(size_t (& sz)[dim]) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			sz[i] = div_c[i] + 2*off[i];
		}
	}

	template<typename ids_type2>
	__device__ __host__ mem_id getGridLinId(const grid_key_dx<dim,ids_type2> & gk) const
	{
		mem_id lid = gk.get(0);
		for (mem_id i = 1 ; i < dim ; i++)
		{
			lid += gk.get(i) * (div_c[i-1] + 2*off[i-1]);
		}

		return lid;
	}

	__device__ __host__ inline grid_key_dx<dim,ids_type> getCell(const Point<dim,T> & xp) const
	{
		return cid_<dim,cnt_type,ids_type,transform>::get_cid_key(spacing_c,off,t,xp);
	}

	__device__ __host__ inline cnt_type LinId(const grid_key_dx<dim,ids_type> & k) const
	{
		return cid_<dim,cnt_type,ids_type,transform>::get_cid(div_c,k);
	}

	__device__ inline const openfpm::array<T,dim,cnt_type> & get_spacing_c() const
	{
		return spacing_c;
	}

	__device__ __host__ inline const openfpm::array<ids_type,dim,cnt_type> & get_div_c() const
	{
		return div_c;
	}

	__device__ __host__ inline const openfpm::array<ids_type,dim,cnt_type> & get_off() const
	{
		return off;
	}

	__device__ __host__ inline const transform & get_t() const
	{
		return t;
	}

	__device__ __host__ inline grid_key_dx<dim> getCellGrid(const T (& pos)[dim]) const
	{
		grid_key_dx<dim> key;
		key.set_d(0,ConvertToID(pos,0));

		for (size_t s = 1 ; s < dim ; s++)
		{
			key.set_d(s,ConvertToID(pos,s));
		}

		return key;
	}

	__device__ __host__ inline grid_key_dx<dim> getCellGrid(const Point<dim,T> & pos) const
	{
		grid_key_dx<dim> key;
		key.set_d(0,ConvertToID(pos,0));

		for (size_t s = 1 ; s < dim ; s++)
		{
			key.set_d(s,ConvertToID(pos,s));
		}

		return key;
	}

	__device__ __host__ inline size_t ConvertToID(const T (&x)[dim] ,size_t s) const
	{
		size_t id = openfpm::math::size_t_floor(t.transform(x,s) / box_unit.getHigh(s)) + off[s];
		id = (id >= gr_cell.size(s))?(gr_cell.size(s)-1-cell_shift.get(s)):id-cell_shift.get(s);
		return id;
	}

	__device__ __host__ inline size_t ConvertToID(const Point<dim,T> & x ,size_t s, size_t sc = 0) const
	{
		size_t id = openfpm::math::size_t_floor(t.transform(x,s) / box_unit.getHigh(s)) + off[s];
		id = (id >= gr_cell.size(s))?(gr_cell.size(s)-1-cell_shift.get(s)):id-cell_shift.get(s);
		return id;
	}
};

#endif /* CELLDECOMPOSER_GPU_KER_HPP_ */
