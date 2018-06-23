/*
 * CellList_gpu.hpp
 *
 *  Created on: Jun 11, 2018
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_GPU_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_GPU_HPP_

#include "config.h"

#ifdef CUDA_GPU

#include <cuda_runtime_api.h>
#include "CellDecomposer.hpp"
#include "Vector/map_vector.hpp"
#include "cuda/Cuda_cell_list_util_func.hpp"

constexpr int count = 0;
constexpr int start = 1;

template<unsigned int dim, typename T,  typename Memory, typename cnt_type = int, typename ids_type = short int, typename transform = no_transform<dim,T>>
class CellList_gpu : public CellDecomposer_sm<dim,T,transform>
{
	//! \brief Cell information
	openfpm::vector<aggregate<cnt_type,cnt_type>,Memory,typename memory_traits_inte<aggregate<cnt_type,cnt_type>>::type,memory_traits_inte> cl_n;

	//! \brief particle information
	openfpm::vector<aggregate<ids_type[dim+1]>,Memory,typename memory_traits_inte<aggregate<ids_type[dim+1]>>::type,memory_traits_inte> part_ids;

	//! Spacing
	CudaMemory spacing;

	//! \brief number of sub-divisions in each direction
	CudaMemory div;

	//! Initialize the structures of the data structure
	void InitializeStructures(const size_t (& div)[dim], size_t tot_n_cell)
	{
		spacing.allocate(sizeof(T)*dim);
		this->div.allocate(dim*sizeof(ids_type));

		T (& div_p)[dim] = *static_cast<T (*)[dim]>(this->div.getPointer());
		T (& spacing_p)[dim] = *static_cast<T (*)[dim]>(this->spacing.getPointer());

		for (size_t i = 0 ; i < dim ; i++)
		{
			div_p[i] = div[i];
			spacing_p[i] = this->getCellBox().getP2().get(i);
		}

		// Force to copy into device
		this->spacing.getDevicePointer();
		this->div.getDevicePointer();

		cl_n.resize(tot_n_cell);
	}

public:

	CellList_gpu(const Box<dim,T> & box, const size_t (&div)[dim], const size_t pad = 1)
	{
		Initialize(box,div,pad);
	}

	/*! Initialize the cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param pad padding cell
	 * \param slot maximum number of slot
	 *
	 */
	void Initialize(const Box<dim,T> & box, const size_t (&div)[dim], const size_t pad = 1)
	{
		SpaceBox<dim,T> sbox(box);

		// Initialize point transformation

		Initialize(sbox,div,pad);
	}

	/*! Initialize the cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param pad padding cell
	 * \param slot maximum number of slot
	 *
	 */
	void Initialize(const SpaceBox<dim,T> & box, const size_t (&div)[dim], const size_t pad = 1)
	{
		Matrix<dim,T> mat;
		CellDecomposer_sm<dim,T,transform>::setDimensions(box,div, mat, pad);

		// create the array that store the number of particle on each cell and se it to 0
		InitializeStructures(this->gr_cell.getSize(),this->gr_cell.size());
	}


	/*! \brief construct from a list of particles
	 *
	 * \param pl Particles list
	 *
	 */
	template<typename vector> void construct(vector & pl)
	{
		// First we set the count memory to zero

		CUDA_SAFE(cudaMemset(cl_n.template getDeviceBuffer<0>(),0,cl_n.size()*sizeof(cnt_type)));

		part_ids.resize(pl.size());

		// Than we construct the ids

		auto ite_gpu = pl.getGPUIterator();

		subindex<dim,T,cnt_type,ids_type><<<ite_gpu.wthr,ite_gpu.thr>>>(*static_cast<ids_type (*)[dim]>(div.getDevicePointer()),
																		*static_cast<T (*)[dim]>(spacing.getDevicePointer()),
																		pl.capacity(),
																		pl.size(),
																		static_cast<T *>(pl.template getDeviceBuffer<0>()),
																		static_cast<cnt_type *>(cl_n.template getDeviceBuffer<0>()),
																		static_cast<ids_type *>(part_ids.template getDeviceBuffer<0>()));
	}

};

#endif

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_GPU_HPP_ */
