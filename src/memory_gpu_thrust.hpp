#ifndef MEMORY_GPU_THRUST_HPP
#define MEMORY_GPU_THRUST_HPP

#include <boost/shared_array.hpp>
#include <algorithm>
#include <cstdlib>
#include "reduce_type.hpp"

template<typename T>
class memory_gpu_thrust;

#include "memory_cpu_thrust.hpp"

typedef long int mem_id;

template<typename T>
class memory_gpu_thrust
{
    public:

	  void * dv;
	  void * pdv;
	  void * ppdv;

	public:

	  memory_gpu_thrust();

	  /**
	   *
	   * Create and load a memory block from memory_cpu_thrust
	   *
	   */
	  void load(memory_cpu_thrust<T> & h_vec);

	  /**
	   *
	   * Create and load a memory block from a boost::array
	   *
	   */
	  void load(boost::shared_array<T> & mem);

	  /**
	   *
	   * Sort all the member
	   *
	   */
	  void sort();

	  /**
	   *
	   * Resize the memory block
	   *
	   */
	  void resize(size_t sz);

	  /**
	   *
	   * Compute a reducion
	   *
	   */
	  typename reduce_type<T>::type reduce();

	  /**
	   *
	   * Get the internal memory pointer
	   *
	   */
	  T* getPointer();

	  /**
	   *
	   * Get the internal thrust object
	   *
	   */
	  void * getThrustObj();
//	  void foreach();
//	  void test();
};

#endif
