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

#include "NN/CellList/CellDecomposer.hpp"
#include "Vector/map_vector.hpp"
#include "Cuda_cell_list_util_func.hpp"
#include "util/cuda/scan_cuda.cuh"
#include "NN/CellList/cuda/CellList_gpu_ker.cuh"
#include "util/cuda_util.hpp"
#include "NN/CellList/CellList_util.hpp"

constexpr int count = 0;
constexpr int start = 1;


template<unsigned int dim, typename T,  typename Memory, typename transform = no_transform_only<dim,T>, typename cnt_type = unsigned int, typename ids_type = unsigned short>
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

	//! Sorted domain particles domain or ghost
	vector_cnt_type sorted_domain_particles_dg;

	//! \brief the index of all the domain particles in the sorted vector
	vector_cnt_type sorted_domain_particles_ids;

	//! \brief for each non sorted index it show the index in the ordered vector
	vector_cnt_type non_sorted_to_sorted;

	//! Spacing
	openfpm::array<T,dim,cnt_type> spacing_c;

	//! \brief number of sub-divisions in each direction
	openfpm::array<ids_type,dim,cnt_type> div_c;

	//! \brief cell padding
	openfpm::array<ids_type,dim,cnt_type> off;

	//! scan object
	scan<cnt_type,ids_type> sc;

	//! Additional information in general (used to understand if the cell-list)
	//! has been constructed from an old decomposition
	size_t n_dec;

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

	//! Indicate that this cell list is a gpu type cell-list
	typedef int yes_is_gpu_celllist;

	/*! \brief Copy constructor
	 *
	 *
	 *
	 */
	CellList_gpu(const CellList_gpu<dim,T,Memory,transform,cnt_type,ids_type> & clg)
	{
		cl_n = clg.cl_n;
		cells = clg.cells;
		starts = clg.starts;
		part_ids = clg.part_ids;
		sorted_to_not_sorted = clg.sorted_to_not_sorted;

		spacing_c = clg.spacing_c;
		div_c = clg.div_c;
		off = clg.off;
		g_m = clg.g_m;
		n_dec = clg.n_dec;
	}

	/*! \brief Copy constructor from temporal
	 *
	 *
	 *
	 */
	CellList_gpu(CellList_gpu<dim,T,Memory,transform,cnt_type,ids_type> && clg)
	{
		cl_n.swap(clg.cl_n);
		cells.swap(clg.cells);
		starts.swap(clg.starts);
		part_ids.swap(clg.part_ids);
		sorted_to_not_sorted.swap(clg.sorted_to_not_sorted);

		spacing_c = clg.spacing_c;
		div_c = clg.div_c;
		off = clg.off;
		g_m = clg.g_m;
		n_dec = clg.n_dec;
	}

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

	vector_cnt_type & getNonSortToSort()
	{
		return non_sorted_to_sorted;
	}

	vector_cnt_type & getDomainSortIds()
	{
		return sorted_domain_particles_ids;
	}

	/*! \brief construct from a list of particles
	 *
	 * \warning pl is assumed to be already be in device memory
	 *
	 * \param pl Particles list
	 *
	 */
	template<typename vector, typename vector_prp> void construct(vector & pl, vector & pl_out, vector_prp & pl_prp, vector_prp & pl_prp_out, mgpu::standard_context_t & mgpuContext, size_t g_m = 0)
	{
#ifdef __NVCC__

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
		sc.scan_(cl_n,starts);

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
		non_sorted_to_sorted.resize(pl.size());

		sorted_domain_particles_ids.resize(pl.size());
		sorted_domain_particles_dg.resize(pl.size());

		auto ite = pl.getGPUIterator();

		// Here we test fill cell
		reorder_parts<decltype(pl_prp.toKernel()),
				      decltype(pl.toKernel()),
				      decltype(sorted_to_not_sorted.toKernel()),
				      cnt_type,shift_ph<0,cnt_type>><<<ite.wthr,ite.thr>>>(pl.size(),
				                                                           pl_prp.toKernel(),
				                                                           pl_prp_out.toKernel(),
				                                                           pl.toKernel(),
				                                                           pl_out.toKernel(),
				                                                           sorted_to_not_sorted.toKernel(),
				                                                           non_sorted_to_sorted.toKernel(),
				                                                           static_cast<cnt_type *>(cells.template getDeviceBuffer<0>()));


		ite = sorted_domain_particles_ids.getGPUIterator();

		mark_domain_particles<<<ite.wthr,ite.thr>>>(sorted_to_not_sorted.toKernel(),sorted_domain_particles_ids.toKernel(),sorted_domain_particles_dg.toKernel(),g_m);


		// now we sort the particles
		mergesort((int *)sorted_domain_particles_dg.template getDeviceBuffer<0>(),(int *)sorted_domain_particles_ids.template getDeviceBuffer<0>(),
				         sorted_domain_particles_dg.size(), mgpu::template less_t<int>(), mgpuContext);

#else

		std::cout << "Error: " <<  __FILE__ << ":" << __LINE__ << " you are calling CellList_gpu.construct() this function is suppose must be compiled with NVCC compiler, but it look like has been compiled by the standard system compiler" << std::endl;

#endif

	}

	CellList_gpu_ker<dim,T,cnt_type,ids_type,transform> toKernel()
	{
		return CellList_gpu_ker<dim,T,cnt_type,ids_type,transform>
		       (starts.toKernel(),
		    	sorted_to_not_sorted.toKernel(),
		    	sorted_domain_particles_ids.toKernel(),
		        spacing_c,
		        div_c,
		        off,
		        this->getTransform());
	}

	/*! \brief Clear the structure
	 *
	 *
	 */
	void clear()
	{
		cl_n.clear();
		cells.clear();
		starts.clear();
		part_ids.clear();
		sorted_to_not_sorted.clear();
	}

	/////////////////////////////////////

	//! Ghost marker
	size_t g_m = 0;

	/*! \brief return the ghost marker
	 *
	 * \return ghost marker
	 *
	 */
	inline size_t get_gm()
	{
		return g_m;
	}

	/*! \brief Set the ghost marker
	 *
	 * \param g_m marker
	 *
	 */
	inline void set_gm(size_t g_m)
	{
		this->g_m = g_m;
	}

	/////////////////////////////////////

	/*! \brief Set the n_dec number
	 *
	 * \param n_dec
	 *
	 */
	void set_ndec(size_t n_dec)
	{
		this->n_dec = n_dec;
	}

	/*! \brief Set the n_dec number
	 *
	 * \return n_dec
	 *
	 */
	size_t get_ndec() const
	{
		return n_dec;
	}

	/////////////////////////////////////
};


#endif

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_GPU_HPP_ */
