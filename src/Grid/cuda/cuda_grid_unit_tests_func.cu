#include "config.h"
#include <Grid/map_grid.hpp>
#include "Point_test.hpp"
#include <stdio.h>
#include "cuda_grid_unit_tests_func.cuh"

template<typename grid_type1, typename grid_type2>
__global__ void grid_gradient_vector(grid_type1 g1, grid_type2 g2, ite_gpu<3> ite_gpu)
{
	GRID_ID_3(ite_gpu);

	g2.template get<4>(key)[0] = (g1.template get<0>(key.move(0,1)) - g1.template get<0>(key.move(0,-1))) / 2.0;
	g2.template get<4>(key)[1] = (g1.template get<0>(key.move(1,1)) - g1.template get<0>(key.move(1,-1))) / 2.0;
	g2.template get<4>(key)[2] = (g1.template get<0>(key.move(2,1)) - g1.template get<0>(key.move(2,-1))) / 2.0;
}

template<typename grid_type>
__global__ void grid_fill_vector(grid_type g1,  ite_gpu<3> ite_gpu)
{
	GRID_ID_3(ite_gpu);

	g1.template get<4>(key)[0] = 1.0;
	g1.template get<4>(key)[1] = 2.0;
	g1.template get<4>(key)[2] = 3.0;
}

template<typename grid_type>
__global__ void grid_fill_vector2(grid_type g1,  ite_gpu<3> ite_gpu)
{
	GRID_ID_3(ite_gpu);

	g1.template get<4>(key)[0] = 1001.0;
	g1.template get<4>(key)[1] = 1002.0;
	g1.template get<4>(key)[2] = 1003.0;
}

template<typename grid_type>
__global__ void compute_stencil_grid(grid_type g1, grid_type g2, ite_gpu<3> ite_gpu)
{
	GRID_ID_3(ite_gpu);

	g2.template get<0>(key) = g1.template get<0>(key.move(0,1)) + g1.template get<0>(key.move(0,-1)) +
			                  g1.template get<0>(key.move(1,1)) + g1.template get<0>(key.move(1,-1)) +
							  g1.template get<0>(key.move(2,1)) + g1.template get<0>(key.move(2,-1)) -
							  6.0*g1.template get<0>(key);
}

__global__ void compute_stencil(float * prp_0, float * prp_1, int sz, grid_key_dx<3> start, grid_key_dx<3> stop)
{
	GRID_ID_3_TRAW(start,stop);

	prp_1[tz*sz*sz + ty*sz + tx] = prp_0[tz*sz*sz + ty*sz + tx + 1] + prp_0[tz*sz*sz + ty*sz + tx - 1] +
									   prp_0[tz*sz*sz + (ty + 1)*sz + tx] + prp_0[tz*sz*sz + (ty - 1)*sz + tx] +
									   prp_0[(tz + 1)*sz*sz + ty*sz + tx + 1] + prp_0[(tz - 1)*sz*sz + ty*sz + tx - 1] -
									   6.0*prp_0[tz*sz*sz + ty*sz + tx];
}

__global__ void fill_one(float * prp_0,int sz)
{
    // Thread index
    int tx = threadIdx.x + blockIdx.x * blockDim.x;
    int ty = threadIdx.y + blockIdx.y * blockDim.y;
    int tz = threadIdx.z + blockIdx.z * blockDim.z;

	prp_0[tz*sz*sz + ty*sz + tx] = 1.0f;
}

__global__ void fill_count(float * prp_0,int sz)
{
    // Thread index
    int tx = threadIdx.x + blockIdx.x * blockDim.x;
    int ty = threadIdx.y + blockIdx.y * blockDim.y;
    int tz = threadIdx.z + blockIdx.z * blockDim.z;

	prp_0[tz*sz*sz + ty*sz + tx] = tz*sz*sz + ty*sz + tx;
}

// call compute

void gpu_grid_3D_one(grid_gpu<3,Point_aggr_test> & g)
{
    // Setup execution parameters
    dim3 threads(8,8,8);
    dim3 grid(8,8,8);

    float * prp_0 = (float *)g.getDeviceBuffer<0>();

	CUDA_LAUNCH_DIM3((fill_one), grid, threads ,prp_0,64);
}

// call compute

void gpu_grid_3D_compute(grid_gpu<3,Point_aggr_test> & g)
{
    // Setup execution parameters
    dim3 threads(8,8,8);
    dim3 grid(8,8,8);

    float * prp_0 = (float *)g.getDeviceBuffer<0>();

	CUDA_LAUNCH_DIM3((fill_count), grid, threads, prp_0,64);
}

void gpu_grid_3D_compute_stencil(grid_gpu<3,Point_aggr_test> & g1, grid_gpu<3,Point_aggr_test> & g2,
								 grid_key_dx<3> & start, grid_key_dx<3> & stop)
{
    // Setup execution parameters

    float * prp_0 = (float *)g1.getDeviceBuffer<0>();
    float * prp_1 = (float *)g2.getDeviceBuffer<0>();

    auto gpu_it = g2.getGPUIterator(start,stop);

    CUDA_LAUNCH_DIM3(compute_stencil, gpu_it.wthr, gpu_it.thr,prp_0,prp_1,64,start,stop);
}

void gpu_grid_3D_compute_grid_stencil(grid_gpu<3,Point_aggr_test> & g1, grid_gpu<3,Point_aggr_test> & g2,
		 	 	 	 	 	 	 	 grid_key_dx<3> & start, grid_key_dx<3> & stop)
{
	auto gpu_it = g2.getGPUIterator(start,stop);

	auto g1k = g1.toKernel();
	auto g2k = g2.toKernel();

	CUDA_LAUNCH_DIM3(compute_stencil_grid, gpu_it.wthr, gpu_it.thr,g1k,g2k,gpu_it);
}

void gpu_grid_fill_vector(grid_gpu<3,Point_aggr_test> & g1, grid_key_dx<3> & start, grid_key_dx<3> & stop)
{
	auto gpu_it = g1.getGPUIterator(start,stop);

	CUDA_LAUNCH_DIM3(grid_fill_vector, gpu_it.wthr, gpu_it.thr ,g1.toKernel(),gpu_it);
}

void gpu_grid_fill_vector2(grid_gpu<3,Point_aggr_test> & g1, grid_key_dx<3> & start, grid_key_dx<3> & stop)
{
	auto gpu_it = g1.getGPUIterator(start,stop);

	CUDA_LAUNCH_DIM3(grid_fill_vector2, gpu_it.wthr, gpu_it.thr ,g1.toKernel(),gpu_it);
}

void gpu_grid_gradient_vector(grid_gpu<3,Point_aggr_test> & g1, grid_gpu<3,Point_aggr_test> & g2, grid_key_dx<3> & start, grid_key_dx<3> & stop)
{
	auto gpu_it = g1.getGPUIterator(start,stop);

	CUDA_LAUNCH_DIM3(grid_gradient_vector, gpu_it.wthr, gpu_it.thr ,g1.toKernel(),g2.toKernel(),gpu_it);
}

