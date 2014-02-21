#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <algorithm>
#include <cstdlib>
#include "memory_gpu_thrust.hpp"

typedef long int mem_id;

#define STC_D(dv) (static_cast<thrust::device_vector<T> *>(dv))
#define STC_H(dv) (static_cast<thrust::host_vector<T> *>(dv))

template<class T> void memory_gpu_thrust<T>::load(memory_cpu_thrust<T> & h_vec)
{
	// transfer data to the device
	thrust::copy(STC_D(dv)->begin(), STC_D(dv)->end(), STC_H(h_vec.getThrustObj())->begin());
}

template<class T> void memory_gpu_thrust<T>::load(boost::shared_array<T> & mem)
{
	// transfer data to the device
	thrust::copy(STC_D(dv)->begin(), STC_D(dv)->end(), mem.get());
}

template<class T> void memory_gpu_thrust<T>::sort()
{
	// sort data into device
	thrust::sort(STC_D(dv)->begin(), STC_D(dv)->end());
}

template<class T> void memory_gpu_thrust<T>::resize(size_t sz)
{
	if (dv == NULL)
	{dv = new thrust::device_vector<T>();}

	// resize data vector
	return STC_D(dv)->resize(sz);
}

template<class T> typename reduce_type<T>::type memory_gpu_thrust<T>::reduce()
{
	// sort data into device
	return thrust::reduce(STC_D(dv)->begin(), STC_D(dv)->end());
}

template<class T> T* memory_gpu_thrust<T>::getPointer()
{
	// return the pointer to data

	return thrust::raw_pointer_cast(&(*STC_D(dv))[0]);
}

template<class T>  void * memory_gpu_thrust<T>::getThrustObj()
{
	return static_cast<void *>(dv);
}

// Compiled Forced specialization

template class memory_gpu_thrust<int>;
template class memory_gpu_thrust<float>;
template class memory_gpu_thrust<double>;


