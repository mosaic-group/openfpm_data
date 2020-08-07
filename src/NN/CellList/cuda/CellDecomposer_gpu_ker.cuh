/*
 * CellDecomposer_gpu_ker.hpp
 *
 *  Created on: Apr 28, 2019
 *      Author: i-bird
 */

#ifndef CELLDECOMPOSER_GPU_KER_HPP_
#define CELLDECOMPOSER_GPU_KER_HPP_

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

public:

	__device__ __host__ CellDecomposer_gpu_ker()
	{}

	__device__ __host__ CellDecomposer_gpu_ker(openfpm::array<T,dim,cnt_type> & spacing_c,
	         	 	 	  openfpm::array<ids_type,dim,cnt_type> & div_c,
	         	 	 	  openfpm::array<ids_type,dim,cnt_type> & off,
	         	 	 	  const transform & t)
	:spacing_c(spacing_c),div_c(div_c),off(off),t(t)
	{}

	__host__ grid_sm<dim,void> getGrid()
	{
		size_t sz[dim];

		for (size_t i = 0 ; i < dim ; i++)
		{
			sz[i] = div_c[i] + 2*off[i];
		}

		return grid_sm<dim,void> (sz);
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
};

#endif /* CELLDECOMPOSER_GPU_KER_HPP_ */
