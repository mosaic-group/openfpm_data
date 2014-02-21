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

	public:

	  void load(memory_cpu_thrust<T> & h_vec);
	  void load(boost::shared_array<T> & mem);
	  void sort();
	  void resize(size_t sz);
	  typename reduce_type<T>::type reduce();
	  T* getPointer();
	  void * getThrustObj();
//	  void foreach();
//	  void test();
};

#endif
