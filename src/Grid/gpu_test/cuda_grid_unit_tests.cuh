/*
 * cuda_gpu_compute.cuh
 *
 *  Created on: Sep 29, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_CUDA_GPU_COMPUTE_CUH_
#define OPENFPM_DATA_SRC_GRID_CUDA_GPU_COMPUTE_CUH_


void gpu_grid_3D_compute(grid_gpu<3,Point_test<float>> & g);
void gpu_grid_3D_compute_stencil(grid_gpu<3,Point_test<float>> & g1, grid_gpu<3,Point_test<float>> & g2,
		 	 	 	 	 	 	 grid_key_dx<3> & key1, grid_key_dx<3> & key2);
void gpu_grid_3D_one(grid_gpu<3,Point_test<float>> & g);
void gpu_grid_3D_compute_grid_stencil(grid_gpu<3,Point_test<float>> & g1, grid_gpu<3,Point_test<float>> & g2,
		 	 	 	 	 	 	 	 grid_key_dx<3> & start, grid_key_dx<3> & stop);

#endif /* OPENFPM_DATA_SRC_GRID_CUDA_GPU_COMPUTE_CUH_ */
