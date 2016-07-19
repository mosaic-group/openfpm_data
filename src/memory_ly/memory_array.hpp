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

	T * ptr;

#endif

	//! return the i element
	T get(mem_id i)
	{
		return ptr[i];
	}

	public:

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

	//! operator[]
	T & operator[](mem_id i)
	{
		return ptr[i];
	}

	//! Default constructor
	memory_array()	{};

	/*! \brief Memory array constructor
	 *
	 * \param ptr Memory pointer
	 * \param sz size
	 * \param init specify if the pointer is initialized
	 *
	 */
	memory_array(void * ptr, size_t sz, bool init)
	: ptr(static_cast<T *>(ptr))
	{
#ifdef SE_CLASS2
		check_valid(ptr,sz);
#endif

		// Initialize the constructors

		if (init == false)
			new (ptr)T[sz];
	};
};


#endif
