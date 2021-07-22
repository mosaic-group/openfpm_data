/*
 * map_vector_std_cuda_ker.cuh
 *
 *  Created on: Mar 9, 2019
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_STD_CUDA_KER_CUH_
#define MAP_VECTOR_STD_CUDA_KER_CUH_


/*! \brief grid interface available when on gpu
 *
 * \tparam n_buf number of template buffers
 *
 */

template<typename T, template <typename> class layout_base>
struct vector_custd_ker
{
	typedef vector_custd_ker<T,layout_base> self_type;

	typedef typename apply_transform<layout_base,T>::type T_;

	//! Actual size of the vector, warning: it is not the space allocated in grid
	//! grid size increase by a fixed amount every time we need a vector bigger than
	//! the actually allocated space
	unsigned int v_size;

	//! 1-D static grid
	vector_gpu_ker<T_,layout_base> base;

public:

	//! it define that it is a vector
	typedef int yes_i_am_vector;

	//! Type of the encapsulation memory parameter
	typedef typename memory_traits_inte<aggregate<T_>>::type layout_type;

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef typename grid_base<1,T_,CudaMemory,typename memory_traits_inte<T_>::type>::container container;

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
	__device__ inline auto at(size_t id) const -> decltype(base.template get<0>(0))
	{
		return base.template get<0>(id);
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
	__device__ inline auto at(size_t id) -> decltype(base.template get<0>(0))
	{
		return base.template get<0>(id);
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
	__device__ inline auto get(size_t id) const -> decltype(base.template get<0>(0))
	{
		return base.template get<0>(id);
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
	__device__ inline auto get(size_t id) -> decltype(base.template get<0>(0))
	{
		return base.template get<0>(id);
	}

	/*! \brief Return the pointer to the chunk of memory
	 *
	 * \return the pointer to the chunk of memory
	 *
	 */
	__device__ void * getPointer()
	{
		return &base.template get<0>(0);
	}

	/*! \brief Return the pointer to the chunk of memory
	 *
	 * \return the pointer to the chunk of memory
	 *
	 */
	__device__ const void * getPointer() const
	{
		return &base.template get<0>(0);
	}


	vector_custd_ker()
	{}

	vector_custd_ker(int v_size, const vector_gpu_ker<T_,layout_base> & cpy)
	:v_size(v_size),base(cpy)
	{}
};




#endif /* MAP_VECTOR_STD_CUDA_KER_CUH_ */
