/*
 * memory_stride.hpp
 *
 *  Created on: Aug 17, 2014
 *      Author: Pietro Incardona
 */

#ifndef MEMORY_ARRAY_HPP_
#define MEMORY_ARRAY_HPP_

#include <memory.hpp>

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
	//! Internal pointer
	T * ptr;

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

	//! Constructor
	memory_array(void * ptr, size_t sz)
	: ptr(static_cast<T *>(ptr))
	{
		// Initialize the constructors
//		new (ptr)T[sz];

	};
};


#endif
