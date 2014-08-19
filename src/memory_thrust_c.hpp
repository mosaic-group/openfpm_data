/*
 * memory_thrust_c.hpp
 *
 *  Created on: Aug 17, 2014
 *      Author: Pietro Incardona
 */

/*!
 * \brief This class is a container for the memory interface
 *
 * This class is a container for the memory interface with thrust.
 *
 * \see memory_c for further explanation
 */

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

template<typename T>
class memory_thrust_c
{
	//! compile time specialization object that allocate memory
	thrust::device_vector<T> mem;

	//! object that represent the memory as an array of objects T
	memory_array<T> mem_r;

	//! constructor
	memory_c(){}
};


