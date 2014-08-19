/*
 * memory_stride.hpp
 *
 *  Created on: Aug 17, 2014
 *      Author: Pietro Incardona
 */

/*!
 *
 * \brief This class give a representation to a chunk or memory
 *
 *	This class give a representation to a chunk of memory as an array of T objects
 *
 * \param T This is the object type
 *
 */

#include <memory.hpp>

template<typename T>
class memory_array
{
	T * ptr;
	memory_array(void * ptr)
	: ptr(ptr)
	{};

	//! return the i element
	T get(mem_id i)
	{
		return ptr[i];
	}

	//! Set the internal pointer to the indicated chunk of memory
	void set_pointer(void * ptr_)
	{
		ptr = ptr_;
	}
};



