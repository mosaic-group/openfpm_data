
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <algorithm>
#include <cstdlib>

#include "memory_cpu_thrust.hpp"

typedef long int mem_id;

#define STC_D(dv) (static_cast<thrust::device_vector<T> *>(dv))
#define STC_H(dv) (static_cast<thrust::host_vector<T> *>(dv))

template<class T> void memory_cpu_thrust<T>::load(memory_gpu_thrust<T> & d_vec)
{
	// transfer data to the device
	thrust::copy(STC_H(dv)->begin(), STC_H(dv)->end(), STC_D(d_vec.dv)->begin());
}

template<class T> void memory_cpu_thrust<T>::sort()
{
	// sort data into device
	thrust::sort(STC_H(dv)->begin(), STC_H(dv)->end());
}

template<class T> void memory_gpu_thrust<T>::resize(size_t sz)
{
	// sort data into device
	return STC_H(dv)->resize(sz);
}

template<class T> T memory_cpu_thrust<T>::reduce()
{
	// sort data into device
	return thrust::reduce(STC_H(dv)->begin(), STC_H(dv)->end());
}

template<class T>  void * memory_cpu_thrust<T>::getThrustObj()
{
	return static_cast<void *>(STC_H(dv));
}

// Compiled Forced specialization

template class memory_cpu_thrust<int>;
template class memory_cpu_thrust<float>;
template class memory_cpu_thrust<double>;


