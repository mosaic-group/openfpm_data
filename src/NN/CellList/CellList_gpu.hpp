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
#include "util/cuda/scan_cuda.cuh"

constexpr int count = 0;
constexpr int start = 1;

template<unsigned int dim, typename cnt_type, typename ids_type>
class NN_gpu_it
{
	grid_key_dx<dim,ids_type> cell_act;

	grid_key_dx<dim,ids_type> cell_start;
	grid_key_dx<dim,ids_type> cell_stop;

	const openfpm::vector_gpu_ker<aggregate<cnt_type>> & starts;

	const openfpm::array<ids_type,dim,cnt_type> & div_c;

	const openfpm::array<ids_type,dim,cnt_type> & off;

	cnt_type p_id;
	cnt_type c_id;

	__device__ void SelectValid()
	{
		while (p_id >= starts.template get<0>(c_id+1) && isNext())
		{
			cnt_type id = cell_act.get(0);
			cell_act.set_d(0,id+1);

			//! check the overflow of all the index with exception of the last dimensionality

			int i = 0;
			for ( ; i < dim-1 ; i++)
			{
				size_t id = cell_act.get(i);
				if ((int)id > cell_stop.get(i))
				{
					// ! overflow, increment the next index

					cell_act.set_d(i,cell_start.get(i));
					id = cell_act.get(i+1);
					cell_act.set_d(i+1,id+1);
				}
				else
				{
					break;
				}
			}

			c_id = cid_<dim,cnt_type,ids_type,int>::get_cid(div_c,cell_act);
			p_id = starts.template get<0>(c_id);
		}
	}

public:

	__device__ NN_gpu_it(const grid_key_dx<dim,ids_type> & cell_pos,
			             const openfpm::vector_gpu_ker<aggregate<cnt_type>> & starts,
			             const openfpm::array<ids_type,dim,cnt_type> & div_c,
			             const openfpm::array<ids_type,dim,cnt_type> & off)
	:starts(starts),div_c(div_c),off(off)
	{
		// calculate start and stop

		for (size_t i = 0 ; i < dim ; i++)
		{
			cell_start.set_d(i,cell_pos.get(i) - 1);
			cell_stop.set_d(i,cell_pos.get(i) + 1);
			cell_act.set_d(i,cell_pos.get(i) - 1);
		}

		c_id = cid_<dim,cnt_type,ids_type,int>::get_cid(div_c,cell_start);
		p_id = starts.template get<0>(c_id);

		SelectValid();
	}

	__device__ cnt_type get()
	{
		return p_id;
	}

	__device__ NN_gpu_it<dim,cnt_type,ids_type> & operator++()
	{
		++p_id;

		SelectValid();

		return *this;
	}

	__device__ bool isNext()
	{
		return cell_act.get(dim-1) <= cell_stop.get(dim-1);
	}
};

template<unsigned int dim, typename T, typename cnt_type, typename ids_type, typename transform>
class CellList_gpu_ker
{
	openfpm::vector_gpu_ker<aggregate<cnt_type>> starts;

	//! Spacing
	openfpm::array<T,dim,cnt_type> spacing_c;

	//! \brief number of sub-divisions in each direction
	openfpm::array<ids_type,dim,cnt_type> div_c;

	//! \brief cell offset
	openfpm::array<ids_type,dim,cnt_type> off;

	//! transformation
	transform t;

public:

	CellList_gpu_ker(openfpm::vector_gpu_ker<aggregate<cnt_type>> starts,
					 openfpm::array<T,dim,cnt_type> & spacing_c,
			         openfpm::array<ids_type,dim,cnt_type> & div_c,
			         openfpm::array<ids_type,dim,cnt_type> & off,
			         const transform & t)
	:starts(starts),spacing_c(spacing_c),div_c(div_c),off(off),t(t)
	{}

	inline __device__ grid_key_dx<dim,ids_type> getCell(const Point<dim,T> & xp)
	{
		return cid_<3,cnt_type,ids_type,transform>::get_cid_key(spacing_c,off,t,xp);
	}

	__device__ NN_gpu_it<dim,cnt_type,ids_type> getNNIterator(const grid_key_dx<dim,ids_type> & cid)
	{
		NN_gpu_it<dim,cnt_type,ids_type> ngi(cid,starts,div_c,off);

		return ngi;
	}
};

template<unsigned int dim, typename T,  typename Memory, typename transform = no_transform_only<dim,T>, typename cnt_type = unsigned int, typename ids_type = unsigned char>
class CellList_gpu : public CellDecomposer_sm<dim,T,transform>
{
	typedef openfpm::vector<aggregate<cnt_type>,Memory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> vector_cnt_type;

	//! \brief Number of particles in each cell
	vector_cnt_type cl_n;

	//! \brief for each cell the particles id in it
	vector_cnt_type cells;

	//! \brief Cell scan with + operation of cl_n
	vector_cnt_type starts;

	//! \brief particle ids information the first "dim" is the cell-id in grid coordinates, the last is the local-id inside the cell
	openfpm::vector<aggregate<ids_type[dim+1]>,Memory,typename memory_traits_inte<aggregate<ids_type[dim+1]>>::type,memory_traits_inte> part_ids;

	//! \brief for each sorted index it show the index in the unordered
	vector_cnt_type sorted_to_not_sorted;

	//! Spacing
	openfpm::array<T,dim,cnt_type> spacing_c;

	//! \brief number of sub-divisions in each direction
	openfpm::array<ids_type,dim,cnt_type> div_c;

	//! \brief cell padding
	openfpm::array<ids_type,dim,cnt_type> off;

	//! Initialize the structures of the data structure
	void InitializeStructures(const size_t (& div)[dim], size_t tot_n_cell, size_t pad)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			div_c[i] = div[i];
			spacing_c[i] = this->getCellBox().getP2().get(i);
			off[i] = pad;
		}

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
		InitializeStructures(this->gr_cell.getSize(),this->gr_cell.size(),pad);
	}

	vector_cnt_type & getSortToNonSort()
	{
		return sorted_to_not_sorted;
	}

	/*! \brief construct from a list of particles
	 *
	 * \warning pl is assumed to be already be in device memory
	 *
	 * \param pl Particles list
	 *
	 */
	template<typename vector, typename vector_prp, unsigned int ... prp> void construct(vector & pl, vector_prp & pl_prp, vector_prp & pl_prp_out)
	{
		part_ids.resize(pl.size());

		// Than we construct the ids

		auto ite_gpu = pl.getGPUIterator();

		cl_n.resize(this->gr_cell.size()+1);
		CUDA_SAFE(cudaMemset(cl_n.template getDeviceBuffer<0>(),0,cl_n.size()*sizeof(cnt_type)));

		part_ids.resize(pl.size());

		subindex<dim,T,cnt_type,ids_type><<<ite_gpu.wthr,ite_gpu.thr>>>(div_c,
																		spacing_c,
																		off,
																		this->getTransform(),
																		pl.capacity(),
																		pl.size(),
																		static_cast<T *>(pl.template getDeviceBuffer<0>()),
																		static_cast<cnt_type *>(cl_n.template getDeviceBuffer<0>()),
																		static_cast<ids_type *>(part_ids.template getDeviceBuffer<0>()));

		// now we scan

		scan<cnt_type,ids_type>(cl_n,starts);

		// now we construct the cells

		cells.resize(pl.size());
		auto itgg = part_ids.getGPUIterator();

		fill_cells<dim,cnt_type,ids_type,shift_ph<0,cnt_type>><<<itgg.wthr,itgg.thr>>>(0,
																					   div_c,
																					   off,
																					   part_ids.size(),
																					   part_ids.capacity(),
																					   static_cast<cnt_type *>(starts.template getDeviceBuffer<0>()),
																					   static_cast<ids_type *>(part_ids.template getDeviceBuffer<0>()),
																					   static_cast<cnt_type *>(cells.template getDeviceBuffer<0>()) );

		sorted_to_not_sorted.resize(pl.size());
		auto ite = pl.getGPUIterator();

		// Here we test fill cell
		reorder_parts<decltype(pl_prp.template toGPU<prp...>()),
				      decltype(sorted_to_not_sorted.template toGPU()),
				      cnt_type,shift_ph<0,cnt_type>><<<ite.wthr,ite.thr>>>(pl.size(),
				                                                           pl_prp.template toGPU<prp...>(),
				                                                           pl_prp_out.template toGPU<>(),
				                                                           sorted_to_not_sorted.template toGPU<>(),
				                                                           static_cast<cnt_type *>(cells.template getDeviceBuffer<0>()));


	}

	CellList_gpu_ker<dim,T,cnt_type,ids_type,transform> toGPU()
	{
		return CellList_gpu_ker<dim,T,cnt_type,ids_type,transform>
		       (starts.template toGPU<>(),
		        spacing_c,
		        div_c,
		        off,
		        this->getTransform());
	}

};

#endif

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_GPU_HPP_ */
