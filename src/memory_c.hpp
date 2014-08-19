/*
 * memory_c.hpp
 *
 *  Created on: Aug 17, 2014
 *      Author: Pietro Incardona
 */

/*!
 * \brief This class is a container for the memory interface
 *
 * This class is a container for the memory interface. It give the possibility
 * to have two specialization, one when the memory interface is full known
 * at compile time, and one when is not-known at compile time.
 * It internally store two objects
 *
 * mem is object used to allocate device/host/... memory
 * mem_r is a representation of this memory as an array of objects T
 *
 */

template<typename T, typename D = void *>
class memory_c
{
	//! compile time specialization object that allocate memory
	D mem;

	//! object that represent the memory as an array of objects T
	memory_array<T> mem_r;

	//! constructor
	memory_c(){}
};

template<typename T>
class memory_c<typename T,void *>
{
	//! Reference to an object to allocate memory
	memory & mem;

	//! object that represent the memory as an array of objects T
	memory_array<T> mem_r

	//! constructor
	memory_c(memory & mem)
	:mem(mem)
	{}
};


