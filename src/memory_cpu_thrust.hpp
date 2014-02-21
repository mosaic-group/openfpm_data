#ifndef MEMORY_CPU_THRUST_HPP
#define MEMORY_CPU_THRUST_HPP

#include <algorithm>
#include <cstdlib>

template<typename T>
class memory_cpu_thrust;

#include "memory_gpu_thrust.hpp"

typedef long int mem_id;

template<typename T>
class memory_cpu_thrust
{
	  void * dv;

	public:

	  void load(memory_gpu_thrust<T> & d_vec);
	  void load(boost::shared_array<T> & mem);
	  void sort();
	  void resize(size_t sz);
	  T reduce();
	  void * getThrustObj();
//	  void foreach();
//	  void test();
};

#endif
