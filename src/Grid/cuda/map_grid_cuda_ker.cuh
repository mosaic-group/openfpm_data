/*
 * map_grid_cuda_ker.hpp
 *
 *  Created on: Jun 28, 2018
 *      Author: i-bird
 */

#ifndef MAP_GRID_CUDA_KER_HPP_
#define MAP_GRID_CUDA_KER_HPP_

#include "config.h"
#include "Grid/grid_base_impl_layout.hpp"
#include "util/tokernel_transformation.hpp"
#ifdef CUDA_GPU
#include "memory/CudaMemory.cuh"
#endif
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */

template<typename T_type_src,typename T_type_dst>
struct copy_switch_memory_c_no_cpy
{
	//! encapsulated source object
	const T_type_src & src;
	//! encapsulated destination object
	T_type_dst & dst;


	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline copy_switch_memory_c_no_cpy(const T_type_src & src,
			                   	   	   	     T_type_dst & dst)
	:src(src),dst(dst)
	{
	};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		auto &grid_layout = boost::fusion::at_c<T::value>(dst);

		grid_layout.disable_manage_memory();
		grid_layout.mem = boost::fusion::at_c<T::value>(src).mem;
		grid_layout.mem_r.bind_ref(boost::fusion::at_c<T::value>(src).mem_r);
#ifdef CUDA_GPU
		if (grid_layout.mem)
		{grid_layout.mem_r.set_pointer(((CudaMemory*)grid_layout.mem)->getDevicePointer());}
#endif
	}
};

template<bool inte_or_lin,typename T>
struct grid_gpu_ker_constructor_impl
{
	template<typename ggk_type> static inline void construct(const ggk_type & cpy,ggk_type & this_)
	{
		copy_switch_memory_c_no_cpy<decltype(cpy.get_data_()),decltype(this_.get_data_())> bp_mc(cpy.get_data_(),this_.get_data_());

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(bp_mc);
	}
};

template<typename T>
struct grid_gpu_ker_constructor_impl<false,T>
{
	template<typename ggk_type> static inline void construct(const ggk_type & cpy,ggk_type & this_)
	{
		auto &grid_layout = this_.get_data_();

		grid_layout.disable_manage_memory();
		grid_layout.mem = cpy.get_data_().mem;
		grid_layout.mem_r.bind_ref(cpy.get_data_().mem_r);

#ifdef CUDA_GPU
		if (grid_layout.mem)
		{grid_layout.mem_r.set_pointer(((CudaMemory*)grid_layout.mem)->getDevicePointer());}
#endif
	}
};

template<unsigned int dim, int prp, typename ids_type>
__device__ void fill_grid_error_array_overflow(const void * sptr,grid_key_dx<dim,ids_type> key)
{
#ifdef CUDA_GPU

	int * ptr = (int *)&global_cuda_error_array[0];

	ptr[0] = 1;
    ptr[1] = ((size_t)sptr) & 0xFFFFFFFF;
    ptr[2] = (((size_t)sptr) & 0xFFFFFFFF00000000) >> 32;
	ptr[3] = prp;
	ptr[4] = dim;

	for (int i = 0 ; i < dim ; i++)
	{ptr[i+5] = key.get(i);}

#ifdef __NVCC__

	ptr[5+dim] = blockIdx.x;
	ptr[6+dim] = blockIdx.y;
	ptr[7+dim] = blockIdx.z;

	ptr[8+dim] = blockDim.x;
	ptr[9+dim] = blockDim.y;
	ptr[10+dim] = blockDim.z;

	ptr[11+dim] = threadIdx.x;
	ptr[12+dim] = threadIdx.y;
	ptr[13+dim] = threadIdx.z;

#endif

#endif
}

template<unsigned int dim>
__device__ void fill_grid_error_array(size_t lin_id)
{
#ifdef CUDA_GPU

	int * ptr = (int *)&global_cuda_error_array[0];

	ptr[0] = 1;
	ptr[1] = 1;
	ptr[2] = lin_id;

#endif
}

template<unsigned int dim, typename T, template <typename> class layout_base, typename linearizer>
class grid_gpu_ker_ref;

/*! \brief grid interface available when on gpu
 *
 * \tparam n_buf number of template buffers
 *
 */
template<unsigned int dim, typename T, template <typename> class layout_base, typename linearizer>
class grid_gpu_ker
{
	//! Type T
	typedef typename apply_transform<layout_base,T>::type T_;

	//! grid information
	linearizer g1;

	//! type of layout of the structure
	typedef typename layout_base<T_>::type layout;

	//! layout data
	mutable layout data_;



	/*! \brief Check that the key is inside the grid
	 *
	 * \param key
	 *
	 * \return
	 *
	 */
	template<typename ids_type> __device__ __host__ inline bool check_bound(const grid_key_dx<dim,ids_type> & v1) const
	{
		for (long int i = 0 ; i < dim ; i++)
		{
			if (v1.get(i) >= (long int)getGrid().size(i))
			{return false;}
			else if (v1.get(i) < 0)
			{return false;}
		}
		return true;
	}

	/*! \brief Check that the key is inside the grid
	 *
	 * \param key
	 *
	 * \return true if it is bound
	 *
	 */
	__device__ __host__ inline bool check_bound(size_t v1) const
	{
		return v1 < getGrid().size();
	}

public:

	//! it define that it is a grid
	typedef int yes_i_am_grid;

	//! Type of the value the vector is storing
	typedef T value_type;

	__device__ __host__ grid_gpu_ker()
	{}

	__device__ __host__ grid_gpu_ker(const linearizer & g1)
	:g1(g1)
	{
	}

	__device__ __host__ grid_gpu_ker(const grid_gpu_ker & cpy)
	:g1(cpy.g1)
	{
		grid_gpu_ker_constructor_impl<is_layout_inte<layout_base<T_>>::value,T_>::construct(cpy,*this);
	}

	__device__ __host__ void constructor_impl(const grid_gpu_ker & cpy)
	{
		g1 = cpy.g1;
		grid_gpu_ker_constructor_impl<is_layout_inte<layout_base<T_>>::value,T_>::construct(cpy,*this);
	}

	__device__ __host__ void constructor_impl(const grid_gpu_ker_ref<dim,T,layout_base,linearizer> & cpy)
	{
		g1 = cpy.ggk.g1;
		grid_gpu_ker_constructor_impl<is_layout_inte<layout_base<T_>>::value,T_>::construct(cpy.ggk,*this);
	}

	/*! \brief Return the internal grid information
	 *
	 * Return the internal grid information
	 *
	 * \return the internal grid
	 *
	 */
	__device__ __host__ const grid_sm<dim,void> & getGrid() const
	{
		return g1;
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename ids_type,typename r_type=decltype(layout_base<T_>::template get<p>(data_,g1,grid_key_dx<dim>()))>
	__device__ __host__ inline r_type get(const grid_key_dx<dim,ids_type> & v1)
	{
#ifdef SE_CLASS1
		if (check_bound(v1) == false)
		{fill_grid_error_array_overflow<dim,p>(this->template getPointer<p>(),v1);}
#endif

		return layout_base<T_>::template get<p>(data_,g1,v1);
	}

	/*! \brief Get the const reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the const reference of the element
	 *
	 */
	template <unsigned int p, typename ids_type, typename r_type=decltype(layout_base<T_>::template get<p>(data_,g1,grid_key_dx<dim>()))>
	__device__ __host__ inline r_type get_debug(const grid_key_dx<dim,ids_type> & v1) const
	{
#ifdef SE_CLASS1
		if (check_bound(v1) == false)
		{fill_grid_error_array_overflow<dim,p>(this->template getPointer<p>(),v1);}
#endif

		return layout_base<T_>::template get<p>(data_,g1,v1);
	}

	/*! \brief Get the const reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the const reference of the element
	 *
	 */
	template <unsigned int p, typename ids_type, typename r_type=decltype(layout_base<T_>::template get<p>(data_,g1,grid_key_dx<dim>()))>
	__device__ __host__ inline r_type get(const grid_key_dx<dim,ids_type> & v1) const
	{
#ifdef SE_CLASS1
		if (check_bound(v1) == false)
		{fill_grid_error_array_overflow<dim,p>(this->template getPointer<p>(),v1);}
#endif
		return layout_base<T_>::template get<p>(data_,g1,v1);
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param lin_id linearized element that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(layout_base<T_>::template get_lin<p>(data_,g1,0))>
	__device__ __host__ inline r_type get(const size_t lin_id)
	{
#ifdef SE_CLASS1
		if (check_bound(lin_id) == false)
		{fill_grid_error_array_overflow<p>(this->getPointer(),lin_id);}
#endif
		return layout_base<T_>::template get_lin<p>(data_,g1,lin_id);
	}

	/*! \brief Get the const reference of the selected element
	 *
	 * \param lin_id linearized element that identify the element in the grid
	 *
	 * \return the const reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(layout_base<T_>::template get_lin<p>(data_,g1,0))>
	__device__ __host__ inline const r_type get(size_t lin_id) const
	{
#ifdef SE_CLASS1
		if (check_bound(lin_id) == false)
		{fill_grid_error_array_overflow<p>(this->getPointer(),lin_id);}
#endif
		return layout_base<T_>::template get_lin<p>(data_,g1,lin_id);
	}

	/*! \brief Get the of the selected element as a boost::fusion::vector
	 *
	 * Get the selected element as a boost::fusion::vector
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \see encap_c
	 *
	 * \return an encap_c that is the representation of the object (careful is not the object)
	 *
	 */
	template<typename Tk>
	__device__ inline encapc<dim,T_,layout> get_o(const grid_key_dx<dim,Tk> & v1)
	{
#ifdef SE_CLASS1
		if (check_bound(v1) == false)
		{fill_grid_error_array_overflow<dim,-1>(this->template getPointer<0>(),v1);}
#endif
		return mem_geto<dim,T_,layout_base<T_>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(data_,g1,v1);
	}

	/*! \brief Get the of the selected element as a boost::fusion::vector
	 *
	 * Get the selected element as a boost::fusion::vector
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \see encap_c
	 *
	 * \return an encap_c that is the representation of the object (careful is not the object)
	 *
	 */
	template<typename Tk>
	__device__ inline const encapc<dim,T_,layout> get_o(const grid_key_dx<dim,Tk> & v1) const
	{
#ifdef SE_CLASS1
		if (check_bound(v1) == false)
		{fill_grid_error_array_overflow<dim,-1>(this->template getPointer<0>(),v1);}
#endif
		return mem_geto<dim,T,layout_base<T_>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(const_cast<decltype(this->data_) &>(data_),g1,v1);
	}


	__device__ inline void set(const grid_key_dx<dim> & key1,const grid_gpu_ker<dim,T_,layout_base, linearizer> & g, const grid_key_dx<dim> & key2)
	{
#ifdef SE_CLASS1
		if (check_bound(key1) == false)
		{fill_grid_error_array_overflow<dim,-1>(this->template getPointer<0>(),key1);}

		if (g.check_bound(key2) == false)
		{fill_grid_error_array_overflow<dim,-1>(g.template getPointer<0>(),key2);}

#endif

		this->get_o(key1) = g.get_o(key2);
	}

	template<unsigned int ... prp> __device__ inline void set(const grid_key_dx<dim> & key1,const grid_gpu_ker<dim,T_,layout_base, linearizer> & g, const grid_key_dx<dim> & key2)
	{
#ifdef SE_CLASS1
		if (check_bound(key1) == false)
		{fill_grid_error_array_overflow<dim,-1>(this->template getPointer<0>(),key1);}

		if (g.check_bound(key2) == false)
		{fill_grid_error_array_overflow<dim,-1>(g.template getPointer<0>(),key2);}

#endif

		auto edest = this->get_o(key1);
		auto esrc = g.get_o(key2);

		copy_cpu_encap_encap_prp<decltype(g.get_o(key2)),decltype(this->get_o(key1)),prp...> ec(esrc,edest);

		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prp)>>(ec);
	}

	/*! \brief set an element of the grid
	 *
	 * set an element of the grid
	 *
	 * \param dx is the grid key or the position to set
	 * \param obj value to set
	 *
	 */
	template<typename Memory> __device__ inline void set(grid_key_dx<dim> key1, const encapc<1,T,Memory> & obj)
	{
#ifdef SE_CLASS1
		if (check_bound(key1) == false)
		{fill_grid_error_array_overflow<dim,-1>(this->template getPointer<0>(),key1);}
#endif

		this->get_o(key1) = obj;
	}

	/*! \brief Get the pointer for the property p
	 *
	 * \tparam property p
	 *
	 */
	template<unsigned int p> __device__ __host__ void * getPointer()
	{
		return mem_getpointer<decltype(data_),layout_base<T>>::template getPointer<p>(data_);
	}

	/*! \brief Get the pointer for the property p
	 *
	 * \tparam property p
	 *
	 */
	template<unsigned int p> __device__ __host__ const void * getPointer() const
	{
		return mem_getpointer<decltype(data_),layout_base<T>>::template getPointer<p>(data_);
	}

	/*! \brief operator= this operator absorb the pointers, consider that this object wrap device pointers
	 *
	 * \param object to copy
	 *
	 */
	grid_gpu_ker<dim,T_,layout_base,linearizer> & operator=(const grid_gpu_ker<dim,T_,layout_base,linearizer> & g)
	{
		g1 = g.g1;

		grid_gpu_ker_constructor_impl<is_layout_inte<layout_base<T_>>::value,T_>::construct(g,*this);

		return *this;
	}

	/*! \brief Get an iterator for the GPU
	 *
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	struct ite_gpu<dim> getGPUIterator(grid_key_dx<dim> & key1, grid_key_dx<dim> & key2, size_t n_thr = default_kernel_wg_threads_) const
	{
		return getGPUIterator_impl<dim>(g1,key1,key2,n_thr);
	}

	/*! \brief Get the internal data_ structure
	 *
	 * \return the data_ structure
	 *
	 */
	__device__ __host__ inline layout & get_data_()
	{
		return data_;
	}

	/*! \brief Get the internal data_ structure
	 *
	 * \return the data_ structure
	 *
	 */
	__device__ __host__ inline const layout & get_data_() const
	{
		return data_;
	}
};

// This is an abstraction for reference type. It exist because the compiler by C++ starndard even if we return a reference deduce
// as value. To force as reference we have to create an object grid_gpu_ker_ref that emulate the reference concept 

template<unsigned int dim, typename T, template <typename> class layout_base, typename linearizer>
class grid_gpu_ker_ref
{
	grid_gpu_ker<dim,T,layout_base,linearizer> & ggk;

	typedef typename apply_transform<layout_base,T>::type T_;

	typedef typename layout_base<T_>::type layout;

public:

	//! it define that it is a grid
	typedef int yes_i_am_grid;

	//! Type of the value the vector is storing
	typedef T value_type;

	//! expose the dimansionality as a static const
	static constexpr unsigned int dims = dim;

	__device__ __host__ grid_gpu_ker_ref()
	{}

	__device__ __host__ grid_gpu_ker_ref(grid_gpu_ker<dim,T,layout_base,linearizer> & ggk)
	:ggk(ggk)
	{}


	__device__ __host__ const grid_sm<dim,void> & getGrid() const
	{
		return ggk.getGrid();
	}

	__device__ __host__ size_t size() const
	{
		return ggk.getGrid().size();
	}

	template <unsigned int p, typename ids_type>
	__device__ __host__ inline auto get(const grid_key_dx<dim,ids_type> & v1) -> decltype(ggk.template get<p>(v1))
	{
		return ggk.template get<p>(v1);
	}

	template <unsigned int p, typename ids_type>
	__device__ __host__ inline auto get(const grid_key_dx<dim,ids_type> & v1) const -> decltype(ggk.template get<p>(v1))
	{
		return ggk.template get<p>(v1);
	}

	template <unsigned int p>
	__device__ __host__ inline auto get(const size_t lin_id) -> decltype(ggk.template get<p>(lin_id))
	{
		return ggk.template get<p>(lin_id);
	}

	template <unsigned int p>
	__device__ __host__ inline auto get(size_t lin_id) const -> decltype(ggk.template get<p>(lin_id))
	{
		return ggk.template get<p>(lin_id);
	}

	template<typename Tk>
	__device__ inline auto get_o(const grid_key_dx<dim,Tk> & v1) -> decltype(ggk.get_o(v1))
	{
		return ggk.get_o(v1);
	}

	template<typename Tk>
	__device__ inline auto get_o(const grid_key_dx<dim,Tk> & v1) const  -> decltype(ggk.get_o(v1))
	{
		return ggk.get_o(v1);
	}


	__device__ inline void set(const grid_key_dx<dim> & key1,const grid_gpu_ker<dim,T_,layout_base, linearizer> & g, const grid_key_dx<dim> & key2)
	{
		ggk.set(key1,g,key2);
	}

	template<unsigned int ... prp> __device__ inline void set(const grid_key_dx<dim> & key1,const grid_gpu_ker<dim,T_,layout_base, linearizer> & g, const grid_key_dx<dim> & key2)
	{
		ggk.template set<prp ...>(key1,g,key2);
	}

	template<typename Memory> __device__ inline void set(grid_key_dx<dim> key1, const encapc<1,T,Memory> & obj)
	{
		ggk.set(key1,obj);
	}

	template<unsigned int p> __device__ __host__ void * getPointer()
	{
		return ggk.template getPointer<p>();
	}

	template<unsigned int p> __device__ __host__ const void * getPointer() const
	{
		return ggk.template getPointer<p>();
	}

	grid_gpu_ker<dim,T,layout_base,linearizer> & operator=(const grid_gpu_ker<dim,T,layout_base,linearizer> & g)
	{
		ggk.operator=(g);

		return *this;
	}

	struct ite_gpu<dim> getGPUIterator(grid_key_dx<dim> & key1, grid_key_dx<dim> & key2, size_t n_thr = default_kernel_wg_threads_) const
	{
		return ggk.getGPUIterator(key1,key2,n_thr);
	}

	__device__ __host__ inline layout & get_data_()
	{
		return ggk.get_data_();
	}

	__device__ __host__ inline const layout & get_data_() const
	{
		return ggk.get_data_();
	}

	const grid_gpu_ker_ref & toKernel() const
	{
		return *this;
	}

	friend class grid_gpu_ker<dim,T,layout_base,linearizer>;
};

#endif /* MAP_GRID_CUDA_KER_HPP_ */
