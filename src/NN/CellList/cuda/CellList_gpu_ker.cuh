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

#ifdef USE_LOW_REGISTER_ITERATOR

#ifdef __NVCC__
__constant__ int cells_striding[126];
#endif

template<unsigned int dim, int r_int, unsigned int pr_int, typename ids_type, typename cnt_type>
struct NN_gpu_int_base_lr_impl
{
	cnt_type ca_pnt;

	cnt_type ca_lin;

	cnt_type c_id;

	__device__ inline void init_impl(const grid_key_dx<dim,ids_type> & cell_pos, const openfpm::array<ids_type,dim,cnt_type> & div_c)
	{
#ifdef __NVCC__
		ca_pnt = 0;

		ca_lin = cid_<dim,cnt_type,ids_type,int>::get_cid(div_c,cell_pos);
		c_id = ca_lin + cells_striding[ca_pnt];
#endif
	}

	__device__ inline void SelectValid_impl(const openfpm::array<ids_type,dim,cnt_type> & div_c)
	{
#ifdef __NVCC__
		++ca_pnt;

		c_id = ca_lin + cells_striding[ca_pnt];
#endif
	}

	__device__ inline bool isNext_impl()
	{
		return ca_pnt < pr_int;
	}
};

#endif

template<unsigned int dim, int r_int, unsigned int pr_int, typename ids_type, typename cnt_type>
struct NN_gpu_int_base_hr_impl
{
	grid_key_dx<dim,ids_type> cell_act;

	grid_key_dx<dim,ids_type> cell_start;
	grid_key_dx<dim,ids_type> cell_stop;

	cnt_type c_id;

	__device__ inline void init_impl(const grid_key_dx<dim,ids_type> & cell_pos, const openfpm::array<ids_type,dim,cnt_type> & div_c)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			cell_start.set_d(i,cell_pos.get(i) - r_int);
			cell_stop.set_d(i,cell_pos.get(i) + r_int);
			cell_act.set_d(i,cell_pos.get(i) - r_int);
		}

		c_id = cid_<dim,cnt_type,ids_type,int>::get_cid(div_c,cell_start);
	}

	__device__ inline void SelectValid_impl(const openfpm::array<ids_type,dim,cnt_type> & div_c)
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
	}

	__device__ inline bool isNext_impl()
	{
		return cell_act.get(dim-1) <= cell_stop.get(dim-1);
	}
};

#ifdef USE_LOW_REGISTER_ITERATOR

template<unsigned int dim, int r_int, unsigned int pr_int, typename ids_type, typename cnt_type>
struct NN_gpu_int_base: public NN_gpu_int_base_hr_impl<dim,r_int,pr_int,ids_type,cnt_type>
{};

template<int r_int, unsigned int pr_int, typename ids_type, typename cnt_type>
struct NN_gpu_int_base<2,r_int,pr_int,ids_type,cnt_type>: public NN_gpu_int_base_lr_impl<2,r_int,pr_int,ids_type,cnt_type>
{};

template<int r_int, unsigned int pr_int, typename ids_type, typename cnt_type>
struct NN_gpu_int_base<3,r_int,pr_int,ids_type,cnt_type>: public NN_gpu_int_base_lr_impl<3,r_int,pr_int,ids_type,cnt_type>
{};

#else

template<unsigned int dim, int r_int, unsigned int pr_int, typename ids_type, typename cnt_type>
struct NN_gpu_int_base: public NN_gpu_int_base_hr_impl<dim,r_int,pr_int,ids_type,cnt_type>
{};

#endif

template<unsigned int dim, typename cnt_type, typename ids_type, unsigned int r_int, bool is_sparse>
class NN_gpu_it: public NN_gpu_int_base<dim,r_int,openfpm::math::pow(2*r_int+1,dim),ids_type,cnt_type>
{
	const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & starts;

	const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt;

	const openfpm::array<ids_type,dim,cnt_type> & div_c;

	const openfpm::array<ids_type,dim,cnt_type> & off;

	cnt_type p_id;
	cnt_type p_id_end;

	inline __device__ void SelectValid()
	{
		while (p_id >= p_id_end && isNext())
		{
			this->SelectValid_impl(div_c);

			if (isNext() == false) {break;}

			p_id = starts.template get<0>(this->c_id);
			p_id_end = starts.template get<0>(this->c_id+1);
		}
	}


public:

	inline __device__ NN_gpu_it(const grid_key_dx<dim,ids_type> & cell_pos,
			             const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & starts,
			             const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt,
			             const openfpm::array<ids_type,dim,cnt_type> & div_c,
			             const openfpm::array<ids_type,dim,cnt_type> & off)
	:starts(starts),srt(srt),div_c(div_c),off(off)
	{
		// calculate start and stop

		this->init_impl(cell_pos,div_c);

		p_id = starts.template get<0>(this->c_id);
		p_id_end = starts.template get<0>(this->c_id+1);

		SelectValid();
	}

	inline __device__ cnt_type get_sort()
	{
		return p_id;
	}

	inline __device__ cnt_type get()
	{
		return srt.template get<0>(p_id);
	}

	inline __device__ NN_gpu_it<dim,cnt_type,ids_type,r_int,is_sparse> & operator++()
	{
		++p_id;

		SelectValid();

		return *this;
	}

	inline __device__ cnt_type get_start(unsigned int ce_id)
	{
		return starts.template get<0>(ce_id);
	}

	inline __device__ cnt_type get_cid()
	{
		return this->c_id;
	}

	inline __device__ bool isNext()
	{
		return this->isNext_impl();
	}
};

template<unsigned int dim, typename cnt_type, typename ids_type, unsigned int r_int>
class NN_gpu_it<dim,cnt_type,ids_type,r_int,true>
{
	cnt_type p_id;
	cnt_type p_id_end;

	cnt_type cells_list_start;
	cnt_type cells_list_stop;

	const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt;

	const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & cells_nn;

	const openfpm::vector_gpu_ker<aggregate<cnt_type,cnt_type>,memory_traits_inte> & cell_nn_list;

	__device__ void SelectValid()
	{
		while (p_id >= p_id_end && isNext())
		{
			++cells_list_start;

			if (cells_list_start < cells_list_stop)
			{
				// calculate start and stop
				p_id = cell_nn_list.template get<0>(cells_list_start);
				p_id_end = cell_nn_list.template get<1>(cells_list_start);
			}
		}
	}


public:

	__device__ NN_gpu_it(cnt_type c_id_sparse,
            			const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & cells_nn,
            			const openfpm::vector_gpu_ker<aggregate<cnt_type,cnt_type>,memory_traits_inte> & cell_nn_list,
            			const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt)
	:srt(srt),cells_nn(cells_nn),cell_nn_list(cell_nn_list)
	{
		if (c_id_sparse == (cnt_type)-1)
		{
			cells_list_stop = cells_list_start;
			return;
		}

		cells_list_start = cells_nn.template get<0>(c_id_sparse);
		cells_list_stop = cells_nn.template get<0>(c_id_sparse + 1);

		// calculate start and stop
		p_id = cell_nn_list.template get<0>(cells_list_start);
		p_id_end = cell_nn_list.template get<1>(cells_list_start);

		SelectValid();
	}

	__device__ cnt_type get_sort()
	{
		return p_id;
	}

	__device__ cnt_type get()
	{
		return srt.template get<0>(p_id);
	}

	__device__ NN_gpu_it<dim,cnt_type,ids_type,r_int,true> & operator++()
	{
		++p_id;

		SelectValid();

		return *this;
	}

	__device__ bool isNext()
	{
		return cells_list_start < cells_list_stop;
	}
};

template<unsigned int dim, typename cnt_type, typename ids_type>
class NN_gpu_it_radius
{
	cnt_type pos;

	cnt_type act;

	const openfpm::vector_gpu_ker<aggregate<int>,memory_traits_inte> & cells;

	const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & starts;

	const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt;

	const openfpm::array<ids_type,dim,cnt_type> & div_c;

	const openfpm::array<ids_type,dim,cnt_type> & off;

	cnt_type p_id;
	cnt_type c_id;

	__device__ inline void SelectValid()
	{
		while (isNext() && p_id >= starts.template get<0>(c_id+1))
		{
			act++;

			if (act >= cells.size())
			{break;}

			c_id = pos + cells.template get<0>(act);
			p_id = starts.template get<0>(c_id);
		}
	}


public:

	__device__ inline NN_gpu_it_radius(const grid_key_dx<dim,ids_type> & cell_pos,
						 const openfpm::vector_gpu_ker<aggregate<int>,memory_traits_inte> & cells,
			             const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & starts,
			             const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt,
			             const openfpm::array<ids_type,dim,cnt_type> & div_c,
			             const openfpm::array<ids_type,dim,cnt_type> & off)
	:act(0),cells(cells),starts(starts),srt(srt),div_c(div_c),off(off)
	{
		// calculate start and stop

		pos = cid_<dim,cnt_type,ids_type,int>::get_cid(div_c,cell_pos);
		c_id = pos + cells.template get<0>(act);
		p_id = starts.template get<0>(c_id);

		SelectValid();
	}

	__device__ cnt_type get_sort()
	{
		return p_id;
	}

	__device__ cnt_type get()
	{
		return srt.template get<0>(p_id);
	}

	__device__ NN_gpu_it_radius<dim,cnt_type,ids_type> & operator++()
	{
		++p_id;

		SelectValid();

		return *this;
	}

	__device__ cnt_type get_start(unsigned int ce_id)
	{
		return starts.template get<0>(ce_id);
	}

	__device__ cnt_type get_cid()
	{
		return c_id;
	}

	__device__ bool isNext()
	{
		return act < cells.size();
	}
};

template<unsigned int dim,typename cnt_type,typename ids_type,bool is_sparse>
struct NN_gpu_selector
{
	static NN_gpu_it<dim,cnt_type,ids_type,1,is_sparse> get(grid_key_dx<dim,ids_type> & cid,
															openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & starts,
															openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt,
															openfpm::array<ids_type,dim,cnt_type> div_c,
															openfpm::array<ids_type,dim,cnt_type> off)
	{
		NN_gpu_it<dim,cnt_type,ids_type,1,is_sparse> ngi(cid,starts,srt,div_c,off);

		return ngi;
	}
};

template<unsigned int dim,typename cnt_type,typename ids_type>
struct NN_gpu_selector<dim,cnt_type,ids_type,true>
{
	static NN_gpu_it<dim,cnt_type,ids_type,1,true> get(grid_key_dx<dim,ids_type> & cid,
															cnt_type c_id_sparse,
															openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & starts,
															openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt,
															openfpm::array<ids_type,dim,cnt_type> div_c,
															openfpm::array<ids_type,dim,cnt_type> off,
															const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & cells_nn,
															const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & cell_nn_list)
	{
		NN_gpu_it<dim,cnt_type,ids_type,1,true> ngi(c_id_sparse,cells_nn,cell_nn_list,srt);

		return ngi;
	}
};

template<unsigned int dim, typename T, typename cnt_type, typename ids_type, typename transform, bool is_sparse>
class CellList_gpu_ker: public CellDecomposer_gpu_ker<dim,T,cnt_type,ids_type,transform>
{
	//! starting point for each cell
	openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> starts;

	//! Sorted to non sorted ids conversion
	openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> srt;

	//! Domain particles ids
	openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> dprt;

	//! radius cells
	openfpm::vector_gpu_ker<aggregate<int>,memory_traits_inte> rad_cells;

	//! Ghost particle marker
	unsigned int g_m;

public:

	typedef int yes_is_gpu_ker_celllist;

	//! Indicate this structure has a function to check the device pointer
	typedef int yes_has_check_device_pointer;

	__device__ inline CellList_gpu_ker(openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> starts,
					 openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> srt,
					 openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> dprt,
					 openfpm::vector_gpu_ker<aggregate<int>,memory_traits_inte> rad_cells,
					 openfpm::array<T,dim,cnt_type> & spacing_c,
			         openfpm::array<ids_type,dim,cnt_type> & div_c,
			         openfpm::array<ids_type,dim,cnt_type> & off,
			         const transform & t,
			         unsigned int g_m)
	:CellDecomposer_gpu_ker<dim,T,cnt_type,ids_type,transform>(spacing_c,div_c,off,t),
	 starts(starts),srt(srt),dprt(dprt),rad_cells(rad_cells),g_m(g_m)
	{
	}


	template<unsigned int stub = NO_CHECK>
	inline __device__ NN_gpu_it<dim,cnt_type,ids_type,1,is_sparse> getNNIterator(const grid_key_dx<dim,ids_type> & cid)
	{
		NN_gpu_it<dim,cnt_type,ids_type,1,is_sparse> ngi(cid,starts,srt,this->get_div_c(),this->get_off());

		return ngi;
	}

	inline __device__ NN_gpu_it_radius<dim,cnt_type,ids_type> getNNIteratorRadius(const grid_key_dx<dim,ids_type> & cid)
	{
		NN_gpu_it_radius<dim,cnt_type,ids_type> ngi(cid,rad_cells,starts,srt,this->get_div_c(),this->get_off());

		return ngi;
	}

	template<unsigned int r_int = 2> inline __device__ NN_gpu_it<dim,cnt_type,ids_type,r_int,is_sparse> getNNIteratorBox(const grid_key_dx<dim,ids_type> & cid)
	{
		NN_gpu_it<dim,cnt_type,ids_type,r_int,is_sparse> ngi(cid,starts,srt,this->get_div_c(),this->get_off());

		return ngi;
	}

	inline __device__ openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & getDomainSortIds()
	{
		return dprt;
	}

	inline __device__ openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & getSortToNonSort()
	{
		return srt;
	}

	/*! \brief Get the number of cells this cell-list contain
	 *
	 * \return number of cells
	 */
	inline __device__ unsigned int getNCells() const
	{
		return starts.size() - 1;
	}

	/*! \brief Return the number of elements in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return number of elements in the cell
	 *
	 */
	inline __device__ cnt_type getNelements(const cnt_type cell_id) const
	{
		return starts.template get<0>(cell_id+1) - starts.template get<0>(cell_id);
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
	inline __device__ cnt_type get(size_t cell, size_t ele)
	{
		cnt_type p_id = starts.template get<0>(cell) + ele;
		return srt.template get<0>(p_id);
	}


	inline __device__ unsigned int get_g_m()
	{
		return g_m;
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

		pc = starts.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Cell index overflow (starts): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = srt.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (str): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = dprt.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (dprt): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = rad_cells.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (dprt): ") + "\n" + pc.match_str;
			return pc;
		}

		return pc;
	}

#endif
};


template<unsigned int dim, typename T, typename cnt_type, typename ids_type, typename transform>
class CellList_gpu_ker<dim,T,cnt_type,ids_type,transform,true>: public CellDecomposer_gpu_ker<dim,T,cnt_type,ids_type,transform>
{
	//! starting point for each cell
	openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> cell_nn;

	//! starting point for each cell
	openfpm::vector_gpu_ker<aggregate<cnt_type,cnt_type>,memory_traits_inte> cell_nn_list;

	//! Sorted to non sorted ids conversion
	openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> srt;

	//! Domain particles ids
	openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> dprt;

	//! Set of cells sparse
	openfpm::vector_sparse_gpu_ker<aggregate<cnt_type>,cnt_type,memory_traits_inte> cl_sparse;

	//! Ghost particle marker
	unsigned int g_m;

public:

	//! Indicate this structure has a function to check the device pointer
	typedef int yes_has_check_device_pointer;

	__device__ inline CellList_gpu_ker(openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> cell_nn,
					 openfpm::vector_gpu_ker<aggregate<cnt_type,cnt_type>,memory_traits_inte> cell_nn_list,
					 openfpm::vector_sparse_gpu_ker<aggregate<cnt_type>,cnt_type,memory_traits_inte> cl_sparse,
					 openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> srt,
					 openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> dprt,
					 openfpm::array<T,dim,cnt_type> & spacing_c,
			         openfpm::array<ids_type,dim,cnt_type> & div_c,
			         openfpm::array<ids_type,dim,cnt_type> & off,
			         const transform & t,
			         unsigned int g_m)
	:CellDecomposer_gpu_ker<dim,T,cnt_type,ids_type,transform>(spacing_c,div_c,off,t),cell_nn(cell_nn),cell_nn_list(cell_nn_list),srt(srt),dprt(dprt),
	 cl_sparse(cl_sparse),g_m(g_m)
	{
	}

	inline __device__ auto getCell(const Point<dim,T> & xp) const -> decltype(cl_sparse.get_sparse(0))
	{
		cnt_type cell = cid_<dim,cnt_type,ids_type,transform>::get_cid(this->get_div_c(),this->get_spacing_c(),this->get_off(),this->get_t(),xp);

		return cl_sparse.get_sparse(cell);
	}


	template<unsigned int stub = NO_CHECK>
	inline __device__ NN_gpu_it<dim,cnt_type,ids_type,1,true> getNNIterator(decltype(cl_sparse.get_sparse(0)) cid)
	{
		NN_gpu_it<dim,cnt_type,ids_type,1,true> ngi(cid.id,cell_nn,cell_nn_list,srt);

		return ngi;
	}

	template<unsigned int r_int = 2>
	inline __device__ NN_gpu_it<dim,cnt_type,ids_type,r_int,true>
	getNNIteratorBox(decltype(cl_sparse.get_sparse(0)) cid)
	{
		NN_gpu_it<dim,cnt_type,ids_type,r_int,true> ngi(cid.id,cell_nn,cell_nn_list,srt);

		return ngi;
	}


	inline __device__ openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & getDomainSortIds()
	{
		return dprt;
	}

	inline __device__ openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & getSortToNonSort()
	{
		return srt;
	}


	inline __device__ unsigned int get_g_m()
	{
		return g_m;
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

		pc = cell_nn.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Cell index overflow (starts): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = cell_nn_list.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Cell particle buffer overflow (cell_nn_list): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = srt.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (str): ") + "\n" + pc.match_str;
			return pc;
		}

		pc = dprt.check_device_pointer(ptr);

		if (pc.match == true)
		{
			pc.match_str = std::string("Particle index overflow (dprt): ") + "\n" + pc.match_str;
			return pc;
		}

		return pc;
	}

#endif
};


#endif /* CELLLIST_GPU_KER_CUH_ */
