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
__global__ void copy_two_vectors(vector_src_type v_src, vector_dst_type v_dst)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	v_src.get(i) = v_dst.get(i);
}

#endif

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

	public:

		//! it define that it is a vector
		typedef int yes_i_am_vector;

		//! Type of the encapsulation memory parameter
		typedef typename layout_base<T_>::type layout_type;

		//! Object container for T, it is the return type of get_o it return a object type trough
		// you can access all the properties of T
		typedef typename grid_cpu<1,T_,CudaMemory,typename layout_base<T_>::type>::container container;

		//! Type of the value the vector is storing
		typedef T_ value_type;

		/*! \brief Return the size of the vector
		 *
		 * \return the size
		 *
		 */
		__device__ unsigned int size() const
		{
			return v_size;
		}

		/*! \brief return the maximum capacity of the vector before reallocation
		 *
		 * \return the capacity of the vector
		 *
		 */

		__device__ unsigned int capacity() const
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
		__device__ inline auto get(size_t id) const -> decltype(base.template get<p>(grid_key_dx<1>(0)))
		{
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
		inline __device__ auto get(size_t id) -> decltype(base.get_o(grid_key_dx<1>(id)))
		{
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
		inline __device__ auto get(size_t id) const -> const decltype(base.get_o(grid_key_dx<1>(id)))
		{
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

		inline __device__ auto get_o(size_t id) const -> decltype(base.get_o(id))
		{
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

		inline __device__ auto get_o(size_t id) -> decltype(base.get_o(id))
		{
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
		__device__ __host__ inline auto get(size_t id) -> decltype(base.template get<p>(grid_key_dx<1>(0)))
		{
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
		template <typename encap_S, unsigned int ...args> void set_o(size_t i, const encap_S & obj)
		{
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
		__device__ void set(size_t id, const vector_gpu_ker<T_,layout_base> & v, size_t src)
		{
			base.set(id,v.base,src);
		}

		/*! \brief Set the element of the vector v from another element of another vector
		 *
		 * \param id element id
		 * \param v vector source
		 * \param src source element
		 *
		 */
		template<unsigned int ... prp> __device__ void set(size_t id, const vector_gpu_ker<T_,layout_base> & v, size_t src)
		{
			base.template set<prp...>(id,v.base,src);
		}

		/*! \brief operator= this operator absorb the pointers, consider that this object wrap device pointers
		 *
		 * \param object to copy
		 *
		 */
		vector_gpu_ker<T,layout_base> & operator=(const vector_gpu_ker<T,layout_base> & v)
		{
			v_size = v.v_size;
			base = v.base;

			return *this;
		}

		void * internal_get_size_pointer()	{return &v_size;}

		void print_size()
		{
			std::cout << "the size of: " << demangle(typeid(self_type).name()) << " is " << sizeof(self_type) << std::endl;
			std::cout << "    " << demangle(typeid(decltype(v_size)).name()) << ":" << sizeof(decltype(v_size)) << std::endl;
			std::cout << "    " << demangle(typeid(decltype(base)).name()) << ":" << sizeof(decltype(base)) << std::endl;
		}
	};

}

#endif /* MAP_VECTOR_CUDA_HPP_ */
