/*
 * CellDecomposer_gpu_ker.hpp
 *
 *  Created on: Apr 28, 2019
 *      Author: i-bird
 */

#ifndef CELLDECOMPOSER_GPU_KER_HPP_
#define CELLDECOMPOSER_GPU_KER_HPP_

#include "util/cuda_util.hpp"
#include "util/multi_array_openfpm/array_openfpm.hpp"
#include "Grid/grid_sm.hpp"
#include "NN/CellList/cuda/Cuda_cell_list_util_func.hpp"
#include "NN/CellList/CellDecomposer.hpp"

template <unsigned int dim, typename T, typename ids_type, typename transform_type>
class CellDecomposer_gpu_ker
{
	//! Spacing
	openfpm::array<T,dim> unitCellP2;

	//! \brief number of sub-divisions in each direction
	openfpm::array<ids_type,dim> numCellDiv;

	//! \brief cell offset
	openfpm::array<ids_type,dim> cellPadDim;

	//! transformation
	transform_type pointTransform;

	//! Unit box of the Cell list
	SpaceBox<dim,T> cellListSpaceBox;

	//! Grid structure of the Cell list
	grid_sm<dim,void> cellListGrid;

	//! cellShift
	Point<dim,long int> cellShift;
public:

	__device__ __host__ CellDecomposer_gpu_ker()
	{}

	__device__ __host__ CellDecomposer_gpu_ker(
		openfpm::array<T,dim> & unitCellP2,
		openfpm::array<ids_type,dim> & numCellDiv,
		openfpm::array<ids_type,dim> & cellPadDim,
		const transform_type & pointTransform)
	: unitCellP2(unitCellP2),
	numCellDiv(numCellDiv),
	cellPadDim(cellPadDim),
	pointTransform(pointTransform)
	{}

	__device__ __host__ CellDecomposer_gpu_ker(
		openfpm::array<T,dim> & unitCellP2,
		openfpm::array<ids_type,dim> & numCellDiv,
		openfpm::array<ids_type,dim> & cellPadDim,
		const transform_type & pointTransform,
		SpaceBox<dim,T> cellListSpaceBox,
		grid_sm<dim,void> cellListGrid,
		Point<dim,long int> cellShift)
	: unitCellP2(unitCellP2),
	numCellDiv(numCellDiv),
	cellPadDim(cellPadDim),
	pointTransform(pointTransform),
	cellListSpaceBox(cellListSpaceBox),
	cellListGrid(cellListGrid),
	cellShift(cellShift)
	{}

	__device__ __host__ grid_sm<dim,void> getGrid()
	{
		size_t sz[dim];

		for (size_t i = 0 ; i < dim ; i++)
			sz[i] = numCellDiv[i];

		return grid_sm<dim,void> (sz);
	}

	__device__ __host__ void getGridSize(size_t (& sz)[dim]) const
	{
		for (size_t i = 0 ; i < dim ; i++)
			sz[i] = numCellDiv[i] + 2*cellPadDim[i];
	}

	template<typename ids_type2>
	__device__ __host__ mem_id getGridLinId(const grid_key_dx<dim,ids_type2> & gk) const
	{
		mem_id lid = gk.get(0);
		for (mem_id i = 1 ; i < dim ; i++)
			lid += gk.get(i) * (numCellDiv[i-1] + 2*cellPadDim[i-1]);

		return lid;
	}

	__device__ __host__ inline grid_key_dx<dim,ids_type> getCell(const Point<dim,T> & xp) const
	{
		return cid_<dim,ids_type,transform_type>::get_cid_key(unitCellP2,cellPadDim,pointTransform,xp);
	}

	__device__ __host__ inline unsigned int LinId(const grid_key_dx<dim,ids_type> & k) const
	{
		return cid_<dim,ids_type,transform_type>::get_cid(numCellDiv,k);
	}

	__device__ inline const openfpm::array<T,dim> & get_spacing_c() const
	{
		return unitCellP2;
	}

	__device__ __host__ inline const openfpm::array<ids_type,dim> & get_div_c() const
	{
		return numCellDiv;
	}

	__device__ __host__ inline const openfpm::array<ids_type,dim> & get_off() const
	{
		return cellPadDim;
	}

	__device__ __host__ inline const transform_type & get_t() const
	{
		return pointTransform;
	}

	__device__ __host__ inline grid_key_dx<dim> getCellGrid(const T (& pos)[dim]) const
	{
		grid_key_dx<dim> key;
		key.set_d(0,ConvertToID(pos,0));

		for (size_t s = 1 ; s < dim ; s++)
			key.set_d(s,ConvertToID(pos,s));

		return key;
	}

	__device__ __host__ inline grid_key_dx<dim> getCellGrid(const Point<dim,T> & pos) const
	{
		grid_key_dx<dim> key;
		key.set_d(0,ConvertToID(pos,0));

		for (size_t s = 1 ; s < dim ; s++)
			key.set_d(s,ConvertToID(pos,s));

		return key;
	}

	__device__ __host__ inline size_t ConvertToID(const T (&x)[dim] ,size_t s) const
	{
		size_t id = openfpm::math::size_t_floor(pointTransform.transform(x,s) / cellListSpaceBox.getHigh(s)) + cellPadDim[s];
		id = (id >= cellListGrid.size(s))?(cellListGrid.size(s)-1-cellShift.get(s)):id-cellShift.get(s);
		return id;
	}

	__device__ __host__ inline size_t ConvertToID(const Point<dim,T> & x ,size_t s, size_t sc = 0) const
	{
		size_t id = openfpm::math::size_t_floor(pointTransform.transform(x,s) / cellListSpaceBox.getHigh(s)) + cellPadDim[s];
		id = (id >= cellListGrid.size(s))?(cellListGrid.size(s)-1-cellShift.get(s)):id-cellShift.get(s);
		return id;
	}
};

#endif /* CELLDECOMPOSER_GPU_KER_HPP_ */
