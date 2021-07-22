/*
 * map_vector_cuda.hpp
 *
 *  Created on: Jun 28, 2018
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_CUDA_HPP_
#define MAP_VECTOR_CUDA_HPP_

#ifdef __NVCC__

template<typename vector_src_type, typename vector_dst_type, unsigned int ... args>
__global__ void merge_add_prp_device_impl(vector_src_type v_src, vector_dst_type v_dst, unsigned int old_sz)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i >= v_src.size())
	{return;}

	// write the object in the last element
	object_s_di<decltype(v_src.get(i)),decltype(v_dst.get(old_sz+i)),OBJ_ENCAP,args...>(v_src.get(i),v_dst.get(old_sz+i));
}

template<typename vector_src_type, typename vector_dst_type>
__global__ void copy_two_vectors(vector_src_type v_dst, vector_dst_type v_src)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i >= v_src.size())
	{return;}

	v_dst.get(i) = v_src.get(i);
}


template<template<typename,typename> class op,
         typename vector_src_type,
		 typename vector_dst_type,
		 typename vector_opart_type,
		 unsigned int ... args>
__global__ void merge_add_prp_device_impl_src_dst_opar_offset(vector_src_type v_src, vector_dst_type v_dst, vector_opart_type opart, unsigned int start)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i >= v_src.size())
	{return;}

	// write the object in the last element
	object_s_di_op<op,decltype(v_src.get(0)),decltype(v_dst.get(0)),OBJ_ENCAP,args...>(v_src.get(i),v_dst.get(opart.template get<1>(start + i)));
}

template<template<typename,typename> class op,
         typename vector_src_type,
		 typename vector_dst_type,
		 typename vector_opart_type,
		 unsigned int ... args>
__global__ void merge_add_prp_device_impl_src_offset_dst_opar(vector_src_type v_src, vector_dst_type v_dst, vector_opart_type opart, unsigned int start)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i >= opart.size())
	{return;}

	// write the object in the last element
	object_si_di_op<op,decltype(v_src.get(0)),decltype(v_dst.get(0)),OBJ_ENCAP,args...>(v_src.get(start + i),v_dst.get(opart.template get<0>(i)));
}

#endif


template<int prp>
__device__ void fill_vector_error_array_overflow(const void * sptr,int key)
{
#ifdef CUDA_GPU

	int * ptr = (int *)&global_cuda_error_array[0];

	ptr[0] = 1;
    ptr[1] = ((size_t)sptr) & 0xFFFFFFFF;
    ptr[2] = (((size_t)sptr) & 0xFFFFFFFF00000000) >> 32;
	ptr[3] = prp;
	ptr[4] = 1;

	for (int i = 0 ; i < 1 ; i++)
	{ptr[i+5] = key;}

#ifdef __NVCC__

	ptr[5+1] = blockIdx.x;
	ptr[6+1] = blockIdx.y;
	ptr[7+1] = blockIdx.z;

	ptr[8+1] = blockDim.x;
	ptr[9+1] = blockDim.y;
	ptr[10+1] = blockDim.z;

	ptr[11+1] = threadIdx.x;
	ptr[12+1] = threadIdx.y;
	ptr[13+1] = threadIdx.z;

#endif

#endif
}


namespace openfpm
{

	/*! \brief grid interface available when on gpu
	 *
	 * \tparam n_buf number of template buffers
	 *
	 */

	template<typename T, template <typename> class layout_base>
	struct vector_gpu_ker
	{
		typedef vector_gpu_ker<T,layout_base> self_type;

		typedef typename apply_transform<layout_base,T>::type T_;

		//! Actual size of the vector, warning: it is not the space allocated in grid
		//! grid size increase by a fixed amount every time we need a vector bigger than
		//! the actually allocated space
		unsigned int v_size;

		//! 1-D static grid
		grid_gpu_ker<1,T_,layout_base> base;

		/*! \brief Check that the key is inside the grid
		 *
		 * \param key
		 *
		 * \return true if it is bound
		 *
		 */
		__device__ __host__ inline bool check_bound(size_t v1) const
		{
			return v1 < size();
		}

	public:

		//! it define that it is a vector
		typedef int yes_i_am_vector;

		//! Type of the encapsulation memory parameter
		typedef typename layout_base<T_>::type layout_type;

		//! Object container for T, it is the return type of get_o it return a object type trough
		// you can access all the properties of T
		typedef typename grid_base<1,T_,CudaMemory,typename layout_base<T_>::type>::container container;

		//! Type of the value the vector is storing
		typedef T_ value_type;

		//! Indicate this structure has a function to check the device pointer
		typedef int yes_has_check_device_pointer;

		/*! \brief Return the size of the vector
		 *
		 * \return the size
		 *
		 */
		__device__ __host__ unsigned int size() const
		{
			return v_size;
		}

		/*! \brief return the maximum capacity of the vector before reallocation
		 *
		 * \return the capacity of the vector
		 *
		 */

		__device__ __host__ unsigned int capacity() const
		{
			return base.size();
		}


		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \tparam p Property to get
		 * \param id Element to get
		 *
		 * \return the element value requested
		 *
		 */
		template <unsigned int p>
		__device__ __host__ inline auto get(unsigned int id) const -> decltype(base.template get<p>(grid_key_dx<1>(0)))
		{
#ifdef SE_CLASS1
			if (check_bound(id) == false)
			{fill_vector_error_array_overflow<p>(this->getPointer<p>(),id);}
#endif
			grid_key_dx<1> key(id);

			return base.template get<p>(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */
		__device__ __host__ inline auto get(unsigned int id) -> decltype(base.get_o(grid_key_dx<1>(id)))
		{
#ifdef SE_CLASS1
			if (check_bound(id) == false)
			{fill_vector_error_array_overflow<-1>(this->template getPointer<0>(),id);}
#endif

			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */
		inline __device__ __host__ auto get(unsigned int id) const -> const decltype(base.get_o(grid_key_dx<1>(id)))
		{
#ifdef SE_CLASS1
			if (check_bound(id) == false)
			{fill_vector_error_array_overflow<-1>(this->getPointer<0>(),id);}
#endif

			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * \deprecated
		 *
		 * exactly as get, exist to keep the compatibility with grid
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */

		inline __device__ __host__ auto get_o(unsigned int id) const -> decltype(base.get_o(id))
		{
#ifdef SE_CLASS1
			if (check_bound(id) == false)
			{fill_vector_error_array_overflow<-1>(this->template getPointer<0>(),id);}
#endif

			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * \deprecated
		 *
		 * exactly as get, exist to keep the compatibility with grid
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */

		inline __device__ __host__ auto get_o(unsigned int id) -> decltype(base.get_o(id))
		{
#ifdef SE_CLASS1
			if (check_bound(id) == false)
			{fill_vector_error_array_overflow<-1>(this->template getPointer<0>(),id);}
#endif

			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Get the last element of the vector
		 *
		 * \return the last element (encapsulated)
		 *
		 */
		inline auto last() const -> decltype(base.get_o(0))
		{
			grid_key_dx<1> key(size()-1);

			return base.get_o(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \tparam p Property to get
		 * \param id Element to get
		 *
		 * \return the element value requested
		 *
		 */
		template <unsigned int p>
		__device__ __host__ inline auto get(unsigned int id) -> decltype(base.template get<p>(grid_key_dx<1>(0)))
		{
#ifdef SE_CLASS1
			if (check_bound(id) == false)
			{fill_vector_error_array_overflow<p>(this->template getPointer<p>(),id);}
#endif

			grid_key_dx<1> key(id);

			return base.template get<p>(key);
		}

		/*! \brief Get the last element of the vector
		 *
		 * \return the element (encapsulated)
		 *
		 */
		inline auto last() -> decltype(base.get_o(0))
		{
			grid_key_dx<1> key(size()-1);

			return base.get_o(key);
		}

		vector_gpu_ker()
		:v_size(0)
		{}

		vector_gpu_ker(int v_size, const grid_gpu_ker<1,T_,layout_base> & cpy)
		:v_size(v_size),base(cpy)
		{}


		/*! \brief Set the object id to obj
		 *
		 * \param id element
		 * \param obj object (encapsulated)
		 *
		 */
		__device__ void set(int id, const container & obj)
		{
#ifdef SE_CLASS1
			if (check_bound(id) == false)
			{fill_vector_error_array_overflow<-1>(this->template getPointer<0>(),id);}
#endif

			//! copy the element
			base.set(id,obj);
		}

		/*! \brief Get the pointer for the property p
		 *
		 * \tparam property p
		 *
		 */
		template<unsigned int p> __device__ __host__ void * getPointer()
		{
			//! copy the element
			return base.template getPointer<p>();
		}

		/*! \brief Get the pointer for the property p
		 *
		 * \tparam property p
		 *
		 */
		template<unsigned int p> __device__ __host__ const void * getPointer() const
		{
			//! copy the element
			return base.template getPointer<p>();
		}

		/*! \brief It set an element of the vector from a object that is a subset of the vector properties
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam encap_S object that encapsulate the object
		 * \tparam args ids of the properties to map the object to
		 *
		 * \param i element to set
		 * \param obj object that encapsulate the object
		 *
		 * \param v source vector
		 *
		 */
		template <typename encap_S, unsigned int ...args> void set_o(unsigned int i, const encap_S & obj)
		{
#ifdef SE_CLASS1
			if (check_bound(i) == false)
			{fill_vector_error_array_overflow<-1>(this->template getPointer<0>(),i);}
#endif

			// write the object in the last element
			object_s_di<encap_S,decltype(get(i)),OBJ_ENCAP,args...>(obj,get(i));
		}

		/*! \brief Set the element of the vector v from another element of another vector
		 *
		 * \param id element id
		 * \param v vector source
		 * \param src source element
		 *
		 */
		__device__ void set(unsigned int id, const vector_gpu_ker<T_,layout_base> & v, unsigned int src)
		{
#ifdef SE_CLASS1
			if (check_bound(id) == false)
			{fill_vector_error_array_overflow<-1>(this->template getPointer<0>(),id);}
#endif

			base.set(id,v.base,src);
		}

		/*! \brief Set the element of the vector v from another element of another vector
		 *
		 * \param id element id
		 * \param v vector source
		 * \param src source element
		 *
		 */
		template<unsigned int ... prp>
		__device__ void set(unsigned int id, const vector_gpu_ker<T_,layout_base> & v, unsigned int src)
		{
#ifdef SE_CLASS1
			if (check_bound(id) == false)
			{fill_vector_error_array_overflow<-1>(this->template getPointer<0>(),id);}
#endif

			base.template set<prp...>(id,v.base,src);
		}

		/*! \brief Get an iterator for the GPU
		 *
		 *
		 */
		__host__ ite_gpu<1> getGPUIterator(size_t n_thr = default_kernel_wg_threads_) const
		{
			grid_key_dx<1> start(0);
			grid_key_dx<1> stop(size()-1);

			return base.getGPUIterator(start,stop,n_thr);
		}

		/*! \brief Get an iterator for the GPU
		 *
		 *
		 */
		ite_gpu<1> getGPUIteratorTo(size_t stop, size_t n_thr = default_kernel_wg_threads_) const
		{
			grid_key_dx<1> start(0);
			grid_key_dx<1> stop_(stop);

			return base.getGPUIterator(start,stop_,n_thr);
		}

		/*! \brief operator= this operator absorb the pointers, consider that this object wrap device pointers
		 *
		 * \param object to copy
		 *
		 */
		__host__ vector_gpu_ker<T,layout_base> & operator=(const vector_gpu_ker<T,layout_base> & v)
		{
			v_size = v.v_size;
			base = v.base;

			return *this;
		}

		/*! \brief Return the base
		 *
		 * \return the base
		 *
		 */
		__device__ grid_gpu_ker<1,T_,layout_base> & getBase()
		{
			return base;
		}

		void * internal_get_size_pointer()	{return &v_size;}

		void print_size()
		{
#ifndef DISABLE_ALL_RTTI
			std::cout << "the size of: " << demangle(typeid(self_type).name()) << " is " << sizeof(self_type) << std::endl;
			std::cout << "    " << demangle(typeid(decltype(v_size)).name()) << ":" << sizeof(decltype(v_size)) << std::endl;
			std::cout << "    " << demangle(typeid(decltype(base)).name()) << ":" << sizeof(decltype(base)) << std::endl;
#endif
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
			pc.match = false;

			check_device_ptr<self_type> ptr_chk(ptr,*this);

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::max_prop>>(ptr_chk);

			if (ptr_chk.result == true)
			{
				pc.match = true;
				pc.match_str += std::string("Property: ") + std::to_string(ptr_chk.prp) + "\n";
			}

			return pc;
		}

#endif
	};

}

#endif /* MAP_VECTOR_CUDA_HPP_ */
