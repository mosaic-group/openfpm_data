/*
 * memory_stride.hpp
 *
 *  Created on: Aug 17, 2014
 *      Author: Pietro Incardona
 */

#ifndef MEMORY_ARRAY_HPP_
#define MEMORY_ARRAY_HPP_

#include "memory/memory.hpp"
#include "Memleak_check.hpp"

#ifdef __NVCC__
#else
#define __host__
#define __device__
#endif

/*!
 *
 * \brief This class give a representation to a chunk or memory
 *
 *	This class give a representation to a chunk of memory as an array of T objects
 *
 * \param T This is the object type
 *
 */

template<typename T>
class memory_array
{

#if defined(__GNUG__) || defined(__clang__)

	//! Internal pointer
	T __attribute__((aligned(16))) * ptr;

#else

	//! Internal pointer
	T * ptr;

#endif

	//! number of elements
	size_t sz;

	/*! \brief Get the i element
	 *
	 * \param i element
	 *
	 * \return the element i
	 *
	 */
	T get(mem_id i)
	{
		return ptr[i];
	}


	public:

	/*! \brief Initialize the memory array
	 *
	 * \param ptr pointer
	 * \param sz number of elements in the array
	 * \param init indicate if you have to initialize the memory
	 *
	 * \return the element i
	 *
	 */
	void initialize(void * ptr, size_t sz, bool init)
	{
		this->ptr = static_cast<T *>(ptr);

#ifdef SE_CLASS2
		check_valid(ptr,sz);
#endif

		// Initialize the constructors

		if (init == false)
			new (ptr)T[sz];

		this->sz = sz;
	}

	//! Set the internal pointer to the indicated chunk of memory
	void set_pointer(void * ptr_)
	{
		ptr = static_cast<T *>(ptr_);
	}

	//! Return the pointer
	void * get_pointer()
	{
		return ptr;
	}

	/*! \brief Set from another memory array
	 *
	 * \param the other memory_array
	 *
	 */
	void bind_ref(const memory_array<T> & ref)
	{
		ptr = ref.ptr;
	}

	/*! \brief Access element an element of the array
	 *
	 * \param i element to access
	 *
	 * \return a teference to the object
	 *
	 */
	__host__ __device__ T & operator[](mem_id i)
	{
		return ptr[i];
	}

	/*! \brief Access element an element of the array
	 *
	 * \param i element to access
	 *
	 * \return a teference to the object
	 *
	 */
	__host__ __device__ const T & operator[](mem_id i) const
	{
		return ptr[i];
	}

	/*! \brief swap the two objects memory
	 *
	 * \param obj memory to swap with
	 *
	 */
	void swap(memory_array<T> & obj)
	{
		size_t sz_tmp = sz;
		sz = obj.sz;
		obj.sz = sz_tmp;

		T * ptr_tmp = ptr;
		ptr = obj.ptr;
		obj.ptr = ptr_tmp;
	}

	/*! \brief Deinitialize the memory
	 *
	 *
	 *
	 */
	void deinit()
	{
		// Call the destructor of every objects
		for (size_t i = 0 ; i < sz ; i++)
		{
			(&ptr[i])->~T();
		}
	}

	//! Default constructor
	memory_array()
	:ptr(NULL),sz(0)
	{};

	/*! \brief Memory array constructor
	 *
	 * \param ptr Memory pointer
	 * \param sz size
	 * \param init specify if the pointer is initialized
	 *
	 */
	memory_array(void * ptr, size_t sz, bool init)
	{
		initialize(ptr,sz,init);
	}

	/*! \brief Destructor
	 *
	 *
	 *
	 */
	~memory_array()
	{
	};
};


#endif
