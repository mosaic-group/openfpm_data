/*
 * cuda_launch.hpp
 *
 *  Created on: Jan 14, 2019
 *      Author: i-bird
 */

#ifndef CUDA_LAUNCH_HPP_
#define CUDA_LAUNCH_HPP_


#ifdef CUDA_GPU

#if defined(SE_CLASS1) || defined(CUDA_CHECK_LAUNCH)

#include "cuda_kernel_error_checker.hpp"

#define CUDA_LAUNCH(cuda_call,ite, ...) \
        {\
	    CHECK_SE_CLASS1_PRE\
		cuda_call<<<ite.wthr,ite.thr>>>(__VA_ARGS__); \
		cudaDeviceSynchronize(); \
		{\
			cudaError_t e = cudaGetLastError();\
			if (e != cudaSuccess)\
			{\
				std::string error = cudaGetErrorString(e);\
				std::cout << "Cuda Error in: " << __FILE__ << ":" << __LINE__ << " " << error << std::endl;\
			}\
			CHECK_SE_CLASS1_POST(#cuda_call,__VA_ARGS__)\
		}\
        }

#else

#define CUDA_LAUNCH(cuda_call,ite, ...) \
		cuda_call<<<ite.wthr,ite.thr>>>(__VA_ARGS__);


#endif

#endif

#endif /* CUDA_LAUNCH_HPP_ */
