/*
 * cuda_util.hpp
 *
 *  Created on: Jun 13, 2018
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_CUDA_UTIL_HPP_
#define OPENFPM_DATA_SRC_UTIL_CUDA_UTIL_HPP_

#include "config.h"
#ifdef CUDA_GPU
#include <cuda_runtime.h>
#endif


#ifdef CUDA_GPU

#include "util/cuda/ofp_context.hxx"

	#ifndef __NVCC__


	#else

		#ifndef __host__
		#define __host__
		#define __device__
		#endif

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

	#endif
#else

#ifndef __host__
#define __host__
#define __device__
#endif

namespace mgpu
{
	// Stub class for modern gpu

	struct ofp_context_t
	{
		ofp_context_t(bool print_prop = true, int dev_num = 0)
		{}
	};
}

#endif


#endif /* OPENFPM_DATA_SRC_UTIL_CUDA_UTIL_HPP_ */
