/*
 * cuda_grid_gpu_funcs.cuh
 *
 *  Created on: Aug 20, 2018
 *      Author: i-bird
 */

#ifndef CUDA_GRID_GPU_FUNCS_CUH_
#define CUDA_GRID_GPU_FUNCS_CUH_

#include "map_grid_cuda_ker.cuh"

#if defined(CUDA_GPU) && defined(__NVCC__)

template<unsigned int dim, typename grid_type>
struct copy_ndim_grid_impl
{
	static __device__ void copy(grid_type & src, grid_type & dst)
	{
		unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

		if (i >= src.getGrid().size())
		{return;}

		auto key_src = src.getGrid().InvLinId(i);

		dst.get_o(key_src) = src.get_o(key_src);
	}
};

template<typename grid_type>
struct copy_ndim_grid_impl<2,grid_type>
{
	static __device__ void copy(grid_type & src, grid_type & dst)
	{
		grid_key_dx<2> key_src;
		key_src.set_d(0,threadIdx.x + blockIdx.x * blockDim.x);
		key_src.set_d(1,threadIdx.y + blockIdx.y * blockDim.y);

		dst.get_o(key_src) = src.get_o(key_src);
	}
};

template<typename grid_type>
struct copy_ndim_grid_impl<3,grid_type>
{
	static __device__ void copy(grid_type & src, grid_type & dst)
	{
		grid_key_dx<3> key_src;
		key_src.set_d(0,threadIdx.x + blockIdx.x * blockDim.x);
		key_src.set_d(1,threadIdx.y + blockIdx.y * blockDim.y);
		key_src.set_d(2,threadIdx.y + blockIdx.y * blockDim.y);

		dst.get_o(key_src) = src.get_o(key_src);
	}
};

template<unsigned int dim, typename grid_type>
__global__ void copy_ndim_grid_device(grid_type src, grid_type dst)
{
	copy_ndim_grid_impl<dim,grid_type>::copy(src,dst);
}


#endif


template<bool inte_or_lin,unsigned int dim, typename T>
struct grid_toKernelImpl
{
	template<typename grid_type> static grid_gpu_ker<dim,T,memory_traits_lin> toKernel(grid_type & gc)
	{
		grid_gpu_ker<dim,T,memory_traits_lin> g(gc.getGrid());

		g.data_.mem = gc.get_internal_data_().mem;
		// Increment the reference of mem
		g.data_.mem->incRef();
		g.data_.mem_r.bind_ref(gc.get_internal_data_().mem_r);
		g.data_.switchToDevicePtrNoCopy();

		return g;
	}
};

template<unsigned int dim, typename T>
struct grid_toKernelImpl<true,dim,T>
{
	template<typename grid_type> static grid_gpu_ker<dim,T,memory_traits_inte> toKernel(grid_type & gc)
	{
		grid_gpu_ker<dim,T,memory_traits_inte> g(gc.getGrid());
		copy_switch_memory_c_no_cpy<T> cp_mc(gc.get_internal_data_(),g.data_);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp_mc);

		return g;
	}
};

#endif /* CUDA_GRID_GPU_FUNCS_CUH_ */
