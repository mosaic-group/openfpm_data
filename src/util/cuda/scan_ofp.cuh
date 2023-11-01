/*
 * scan_ofp.hpp
 *
 *  Created on: May 15, 2019
 *      Author: i-bird
 */

#ifndef SCAN_OFP_HPP_
#define SCAN_OFP_HPP_

#ifdef __NVCC__

#include "util/cuda_launch.hpp"
#include "util/ofp_context.hpp"

#if CUDART_VERSION >= 11000
	// Here we have for sure CUDA >= 11
	#ifndef CUDA_ON_CPU
		#ifdef __HIP__
			#include "hipcub/hipcub.hpp"
		#else
			#include "cub/cub.cuh"
		#endif
	#endif
#else
	#include "cub_old/cub.cuh"
#endif


namespace openfpm
{
	template<typename input_it, typename output_it>
			 void scan(input_it input, int count, output_it output, gpu::ofp_context_t& gpuContext)
	{
#ifdef CUDA_ON_CPU

	if (count == 0)	{return;}

	auto prec = input[0];
	output[0] = 0;
	for (int i = 1 ; i < count ; i++)
	{
		auto next = prec + output[i-1];
		prec = input[i];
		output[i] = next;
	}

#else
	if (count == 0)	return;

	#ifdef __HIP__

		size_t temp_storage_bytes = 0;
		hipcub::DeviceScan::ExclusiveSum(NULL,
			temp_storage_bytes,input, output, count);

		auto & temporal = gpuContext.getTemporalCUB();
		temporal.resize(temp_storage_bytes);

		hipcub::DeviceScan::ExclusiveSum(temporal.template getDeviceBuffer<0>(),
			temp_storage_bytes, input, output, count);

	#else

		size_t temp_storage_bytes = 0;
		cub::DeviceScan::ExclusiveSum(NULL,
			temp_storage_bytes, input, output, count);

		auto & temporal = gpuContext.getTemporalCUB();
		temporal.resize(temp_storage_bytes);

		cub::DeviceScan::ExclusiveSum(temporal.template getDeviceBuffer<0>(),
			temp_storage_bytes,input, output, count);

	#endif
#endif
	}
}

#endif /* __NVCC__ */

#endif /* SCAN_OFP_HPP_ */
