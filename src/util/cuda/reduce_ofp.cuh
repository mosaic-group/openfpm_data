/*
 * reduce_ofp.hpp
 *
 *  Created on: May 15, 2019
 *      Author: i-bird
 */

#ifndef REDUCE_OFP_HPP_
#define REDUCE_OFP_HPP_

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
	template<typename input_it, typename output_it, typename reduce_op>
			void reduce(input_it input, int count, output_it output, reduce_op op, gpu::ofp_context_t& context)
	{
#ifdef CUDA_ON_CPU

	output[0] = 0;
	for (int i = 0 ; i < count ; i++)
	{
		output[0] = op(output[0],input[i]);
	}

#else

	#ifdef __HIP__

		size_t temp_storage_bytes = 0;
		hipcub::DeviceReduce::Reduce(NULL,
			temp_storage_bytes,input, output, count, op, false);

		auto & temporal = context.getTemporalCUB();
		temporal.resize(temp_storage_bytes);

		hipcub::DeviceReduce::Reduce(temporal.template getDeviceBuffer<0>(),
			temp_storage_bytes,input, output, count, op, false);
	#else

		size_t temp_storage_bytes = 0;
		cub::DeviceReduce::Reduce(NULL,
			temp_storage_bytes, input, output, count, op, false);

		auto & temporal = context.getTemporalCUB();
		temporal.resize(temp_storage_bytes);

		cub::DeviceReduce::Reduce(temporal.template getDeviceBuffer<0>(),
			temp_storage_bytes, input, output, count, op, false);

	#endif
#endif
	}
}

#endif

#endif /* REDUCE_OFP_HPP_ */
