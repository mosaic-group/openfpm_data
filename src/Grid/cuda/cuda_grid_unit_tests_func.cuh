/*
 * cuda_gpu_compute.cuh
 *
 *  Created on: Sep 29, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_CUDA_GPU_COMPUTE_CUH_
#define OPENFPM_DATA_SRC_GRID_CUDA_GPU_COMPUTE_CUH_

typedef aggregate<float, float, float, float, float [3], float [3][3]> Point_aggr_test;

void gpu_grid_3D_compute(grid_gpu<3,Point_aggr_test> & g);
void gpu_grid_3D_compute_stencil(grid_gpu<3,Point_aggr_test> & g1, grid_gpu<3,Point_aggr_test> & g2,
		 	 	 	 	 	 	 grid_key_dx<3> & key1, grid_key_dx<3> & key2);
void gpu_grid_3D_one(grid_gpu<3,Point_aggr_test> & g);
void gpu_grid_3D_compute_grid_stencil(grid_gpu<3,Point_aggr_test> & g1, grid_gpu<3,Point_aggr_test> & g2,
		 	 	 	 	 	 	 	 grid_key_dx<3> & start, grid_key_dx<3> & stop);

void gpu_grid_fill_vector(grid_gpu<3,Point_aggr_test> & g1, grid_key_dx<3> & start, grid_key_dx<3> & stop);

void gpu_grid_fill_vector2(grid_gpu<3,Point_aggr_test> & g1, grid_key_dx<3> & start, grid_key_dx<3> & stop);

void gpu_grid_gradient_vector(grid_gpu<3,Point_aggr_test> & g1, grid_gpu<3,Point_aggr_test> & g2, grid_key_dx<3> & start, grid_key_dx<3> & stop);


#endif /* OPENFPM_DATA_SRC_GRID_CUDA_GPU_COMPUTE_CUH_ */
