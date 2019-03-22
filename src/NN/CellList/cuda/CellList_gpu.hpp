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
#include "NN/CellList/cuda/CellList_gpu_ker.cuh"
#include "util/cuda_util.hpp"
#include "NN/CellList/CellList_util.hpp"
#include "NN/CellList/CellList.hpp"
#include "Vector/map_vector_sparse.hpp"

#ifdef __NVCC__
#include "util/cuda/moderngpu/kernel_scan.hxx"
#endif

constexpr int count = 0;
constexpr int start = 1;

template<unsigned int dim,
		 typename T,
		 typename Memory,
		 typename transform = no_transform_only<dim,T>,
		 typename cnt_type = unsigned int,
		 typename ids_type = int,
		 bool is_sparse = false>
class CellList_gpu : public CellDecomposer_sm<dim,T,transform>
{
	typedef openfpm::vector<aggregate<cnt_type>,Memory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> vector_cnt_type;

	//! \brief Number of particles in each cell
	vector_cnt_type cl_n;

	//! \brief for each cell the particles id in it
	vector_cnt_type cells;

	//! \brief Cell scan with + operation of cl_n
	vector_cnt_type starts;

	//! \brief sparse vector in case of sparse Cell-list
	openfpm::vector_sparse_gpu<aggregate<cnt_type>> cl_sparse;

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

	//! Radius neighborhood
	openfpm::vector<aggregate<int>,Memory,typename memory_traits_inte<aggregate<int>>::type,memory_traits_inte> nnc_rad;

	//! scan object
	//scan<cnt_type,ids_type> sc;

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

	/*! \brief This function construct a sparse cell-list
	 *
	 *
	 */
	template<typename vector, typename vector_prp>
	void construct_sparse(vector & pl,
  	   	   	 	 	 	  vector & pl_out,
  	   	   	 	 	 	  vector_prp & pl_prp,
  	   	   	 	 	 	  vector_prp & pl_prp_out,
  	   	   	 	 	 	  mgpu::ofp_context_t & mgpuContext,
  	   	   	 	 	 	  size_t g_m = 0,
  	   	   	 	 	 	  cl_construct_opt opt = cl_construct_opt::Full)
	{
#ifdef __NVCC__

		part_ids.resize(pl.size());
		cl_n.resize(pl.size());

		// Than we construct the ids

		auto ite_gpu = pl.getGPUIterator();

		part_ids.resize(pl.size());

		if (ite_gpu.wthr.x == 0)
		{
			// no particles
			starts.resize(cl_n.size());
			starts.template fill<0>(0);
			return;
		}

		CUDA_LAUNCH((subindex<true,dim,T,cnt_type,ids_type>),ite_gpu,div_c,
																		spacing_c,
																		off,
																		this->getTransform(),
																		pl.capacity(),
																		pl.size(),
																		part_ids.capacity(),
																		static_cast<T *>(pl.template getDeviceBuffer<0>()),
																		static_cast<cnt_type *>(cl_n.template getDeviceBuffer<0>()),
																		static_cast<ids_type *>(part_ids.template getDeviceBuffer<0>()));

		// now we construct the cells

		cells.resize(pl.size());

		// Here we fill the sparse vector
		CUDA_LAUNCH((fill_cells_sparse),ite_gpu,cl_sparse.toKernel(),cl_n.toKernel());

		auto itgg = part_ids.getGPUIterator();

		CUDA_LAUNCH((fill_cells<dim,cnt_type,ids_type,shift_ph<0,cnt_type>>),itgg,0,
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

		auto ite = pl.getGPUIterator(64);

		// Here we reorder the particles to improve coalescing access
		CUDA_LAUNCH((reorder_parts<decltype(pl_prp.toKernel()),
				      decltype(pl.toKernel()),
				      decltype(sorted_to_not_sorted.toKernel()),
				      cnt_type,shift_ph<0,cnt_type>>),ite,pl.size(),
				                                                           pl_prp.toKernel(),
				                                                           pl_prp_out.toKernel(),
				                                                           pl.toKernel(),
				                                                           pl_out.toKernel(),
				                                                           sorted_to_not_sorted.toKernel(),
				                                                           non_sorted_to_sorted.toKernel(),
				                                                           static_cast<cnt_type *>(cells.template getDeviceBuffer<0>()));


		if (opt == cl_construct_opt::Full)
		{
			ite = sorted_domain_particles_ids.getGPUIterator();

			CUDA_LAUNCH((mark_domain_particles),ite,sorted_to_not_sorted.toKernel(),sorted_domain_particles_ids.toKernel(),sorted_domain_particles_dg.toKernel(),g_m);


			// now we sort the particles
			mergesort((int *)sorted_domain_particles_dg.template getDeviceBuffer<0>(),(int *)sorted_domain_particles_ids.template getDeviceBuffer<0>(),
							 sorted_domain_particles_dg.size(), mgpu::template less_t<int>(), mgpuContext);
		}

	#else

			std::cout << "Error: " <<  __FILE__ << ":" << __LINE__ << " you are calling CellList_gpu.construct() this function is suppose must be compiled with NVCC compiler, but it look like has been compiled by the standard system compiler" << std::endl;

	#endif
	}

	/*! \brief This function construct a dense cell-list
	 *
	 *
	 */
	template<typename vector, typename vector_prp>
	void construct_dense(vector & pl,
			   	   	   	 vector & pl_out,
			   	   	   	 vector_prp & pl_prp,
			   	   	   	 vector_prp & pl_prp_out,
			   	   	   	 mgpu::ofp_context_t & mgpuContext,
			   	   	   	 size_t g_m = 0,
			   	   	   	 cl_construct_opt opt = cl_construct_opt::Full)
	{
#ifdef __NVCC__

		part_ids.resize(pl.size());

		// Than we construct the ids

		auto ite_gpu = pl.getGPUIterator();

		cl_n.resize(this->gr_cell.size()+1);
		CUDA_SAFE(cudaMemset(cl_n.template getDeviceBuffer<0>(),0,cl_n.size()*sizeof(cnt_type)));

		part_ids.resize(pl.size());

		if (ite_gpu.wthr.x == 0)
		{
			// no particles
			starts.resize(cl_n.size());
			starts.template fill<0>(0);
			return;
		}

		CUDA_LAUNCH((subindex<false,dim,T,cnt_type,ids_type>),ite_gpu,div_c,
																		spacing_c,
																		off,
																		this->getTransform(),
																		pl.capacity(),
																		pl.size(),
																		part_ids.capacity(),
																		static_cast<T *>(pl.template getDeviceBuffer<0>()),
																		static_cast<cnt_type *>(cl_n.template getDeviceBuffer<0>()),
																		static_cast<ids_type *>(part_ids.template getDeviceBuffer<0>()));

		// now we scan
		starts.resize(cl_n.size());
		mgpu::scan((cnt_type *)cl_n.template getDeviceBuffer<0>(), cl_n.size(), (cnt_type *)starts.template getDeviceBuffer<0>() , mgpuContext);

		// now we construct the cells

		cells.resize(pl.size());
		auto itgg = part_ids.getGPUIterator();

		CUDA_LAUNCH((fill_cells<dim,cnt_type,ids_type,shift_ph<0,cnt_type>>),itgg,0,
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

		auto ite = pl.getGPUIterator(64);

		// Here we reorder the particles to improve coalescing access
		CUDA_LAUNCH((reorder_parts<decltype(pl_prp.toKernel()),
				      decltype(pl.toKernel()),
				      decltype(sorted_to_not_sorted.toKernel()),
				      cnt_type,shift_ph<0,cnt_type>>),ite,pl.size(),
				                                                           pl_prp.toKernel(),
				                                                           pl_prp_out.toKernel(),
				                                                           pl.toKernel(),
				                                                           pl_out.toKernel(),
				                                                           sorted_to_not_sorted.toKernel(),
				                                                           non_sorted_to_sorted.toKernel(),
				                                                           static_cast<cnt_type *>(cells.template getDeviceBuffer<0>()));


		if (opt == cl_construct_opt::Full)
		{
			ite = sorted_domain_particles_ids.getGPUIterator();

			CUDA_LAUNCH((mark_domain_particles),ite,sorted_to_not_sorted.toKernel(),sorted_domain_particles_ids.toKernel(),sorted_domain_particles_dg.toKernel(),g_m);


			// now we sort the particles
			mergesort((int *)sorted_domain_particles_dg.template getDeviceBuffer<0>(),(int *)sorted_domain_particles_ids.template getDeviceBuffer<0>(),
							 sorted_domain_particles_dg.size(), mgpu::template less_t<int>(), mgpuContext);
		}

	#else

			std::cout << "Error: " <<  __FILE__ << ":" << __LINE__ << " you are calling CellList_gpu.construct() this function is suppose must be compiled with NVCC compiler, but it look like has been compiled by the standard system compiler" << std::endl;

	#endif
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
	:CellDecomposer_sm<dim,T,transform>(clg)
	{
		cl_n = clg.cl_n;
		cells = clg.cells;
		starts = clg.starts;
		part_ids = clg.part_ids;
		sorted_to_not_sorted = clg.sorted_to_not_sorted;
		sorted_domain_particles_dg = clg.sorted_domain_particles_dg;
		sorted_domain_particles_ids = clg.sorted_domain_particles_ids;
		non_sorted_to_sorted = clg.non_sorted_to_sorted;

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
	:CellDecomposer_sm<dim,T,transform>(clg)
	{
		cl_n.swap(clg.cl_n);
		cells.swap(clg.cells);
		starts.swap(clg.starts);
		part_ids.swap(clg.part_ids);
		sorted_to_not_sorted.swap(clg.sorted_to_not_sorted);
		sorted_domain_particles_dg.swap(clg.sorted_domain_particles_dg);
		sorted_domain_particles_ids.swap(clg.sorted_domain_particles_ids);
		non_sorted_to_sorted.swap(clg.non_sorted_to_sorted);

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


	/*! \brief Set the radius for the getNNIteratorRadius
	 *
	 * \param radius
	 *
	 */
	void setRadius(T radius)
	{
		openfpm::vector<long int> nnc_rad_;

		NNcalc_rad(radius,nnc_rad_,this->getCellBox(),this->getGrid());

		nnc_rad.resize(nnc_rad_.size());

		// copy to nnc_rad

		for (unsigned int i = 0 ; i < nnc_rad_.size() ; i++)
		{nnc_rad.template get<0>(i) = nnc_rad_.template get<0>(i);}

		nnc_rad.template hostToDevice<0>();
	}

	template<typename vector, typename vector_prp>
	void construct_sparse()
	{
		///////////////////////////// In case of sparse ///////////////////////////

		// 1 sub index does add elements to linear_sparse
		// flush will construct the cells structure + offsets + internal index
		//
		//
		// we construct for each cell the NN cells the number of neigborhood cells
		//
		//
		// toKernelSparse()

		///////////////////////////////////////////////////////////////////////////
	}

	/*! \brief construct from a list of particles
	 *
	 * \warning pl is assumed to be already be in device memory
	 *
	 * \param pl Particles list
	 *
	 */
	template<typename vector, typename vector_prp>
	void construct(vector & pl,
				   vector & pl_out,
				   vector_prp & pl_prp,
				   vector_prp & pl_prp_out,
				   mgpu::ofp_context_t & mgpuContext,
				   size_t g_m = 0,
				   cl_construct_opt opt = cl_construct_opt::Full)
	{
		if (is_sparse == false) {construct_dense(pl,pl_out,pl_prp,pl_prp_out,mgpuContext,g_m,opt);}
		else {construct_sparse(pl,pl_out,pl_prp,pl_prp_out,mgpuContext,g_m,opt);}
	}

	CellList_gpu_ker<dim,T,cnt_type,ids_type,transform> toKernel()
	{
		if (nnc_rad.size() == 0)
		{
			// set the radius equal the cell spacing on direction X
			// (must be initialized to something to avoid warnings)
			setRadius(this->getCellBox().getHigh(0));
		}

		return CellList_gpu_ker<dim,T,cnt_type,ids_type,transform>
		       (starts.toKernel(),
		    	sorted_to_not_sorted.toKernel(),
		    	sorted_domain_particles_ids.toKernel(),
		    	nnc_rad.toKernel(),
		        spacing_c,
		        div_c,
		        off,
		        this->getTransform(),
		        g_m);
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

	/*! \brief Transfer the information computed on gpu to construct the cell-list on gpu
	 *
	 */
	void debug_deviceToHost()
	{
		cl_n.template deviceToHost<0>();
		cells.template deviceToHost<0>();
		starts.template deviceToHost<0>();
	}

	/*! \brief Return the numbers of cells contained in this cell-list
	 *
	 * \return the number of cells
	 *
	 */
	size_t getNCells()
	{
		return cl_n.size();
	}

	/*! \brief Return the numbers of elements in the cell
	 *
	 * \return the number of elements in the cell
	 *
	 */
	size_t getNelements(size_t i)
	{
		return cl_n.template get<0>(i);
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
	inline auto get(size_t cell, size_t ele) -> decltype(cells.template get<0>(starts.template get<0>(cell)+ele))
	{
		return cells.template get<0>(starts.template get<0>(cell)+ele);
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
	inline auto get(size_t cell, size_t ele) const -> decltype(cells.template get<0>(starts.template get<0>(cell)+ele))
	{
		return cells.template get<0>(starts.template get<0>(cell)+ele);
	}

	/*! \brief swap the information of the two cell-lists
	 *
	 *
	 *
	 */
	void swap(CellList_gpu<dim,T,Memory,transform,cnt_type,ids_type> & clg)
	{
		((CellDecomposer_sm<dim,T,transform> *)this)->swap(clg);
		cl_n.swap(clg.cl_n);
		cells.swap(clg.cells);
		starts.swap(clg.starts);
		part_ids.swap(clg.part_ids);
		sorted_to_not_sorted.swap(clg.sorted_to_not_sorted);
		sorted_domain_particles_dg.swap(clg.sorted_domain_particles_dg);
		sorted_domain_particles_ids.swap(clg.sorted_domain_particles_ids);
		non_sorted_to_sorted.swap(clg.non_sorted_to_sorted);

		spacing_c.swap(clg.spacing_c);
		div_c.swap(clg.div_c);
		off.swap(clg.off);

		size_t g_m_tmp = g_m;
		g_m = clg.g_m;
		clg.g_m = g_m_tmp;

		size_t n_dec_tmp = n_dec;
		n_dec = clg.n_dec;
		clg.n_dec = n_dec_tmp;
	}
};


#endif

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_GPU_HPP_ */
