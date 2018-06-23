/*
 * cuda_util.hpp
 *
 *  Created on: Jun 13, 2018
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_CUDA_UTIL_HPP_
#define OPENFPM_DATA_SRC_UTIL_CUDA_UTIL_HPP_

#ifdef CUDA_GPU


#define CUDA_SAFE(cuda_call) \
cuda_call; \
{\
	cudaError_t e = cudaPeekAtLastError();\
    if (e != cudaSuccess)\
    {\
        std::string error = cudaGetErrorString(e);\
        std::cout << "Cuda Error in: " << __FILE__ << ":" << __LINE__ << " " << error << std::endl;\
    }\
}

#else
#define __host__
#define __device__
#endif

#endif /* OPENFPM_DATA_SRC_UTIL_CUDA_UTIL_HPP_ */
