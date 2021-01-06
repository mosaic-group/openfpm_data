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

#if CUDART_VERSION >= 11000
	#ifndef CUDA_ON_CPU 
	// Here we have for sure CUDA >= 11
	#include "cub/cub.cuh"
	#ifndef REDUCE_WITH_CUB
		#define REDUCE_WITH_CUB
	#endif
	#endif
#else
	// Here we have old CUDA
	#include "cub_old/cub.cuh"
	#include "util/cuda/moderngpu/kernel_reduce.hxx"
#endif

#include "util/cuda/ofp_context.hxx"

namespace openfpm
{
	template<typename input_it, typename output_it, typename reduce_op>
			void reduce(input_it input, int count, output_it output, reduce_op op, mgpu::ofp_context_t& context)
	{
#ifdef CUDA_ON_CPU

	output[0] = 0;
	for (int i = 0 ; i < count ; i++)
	{
		output[0] = op(output[0],input[i]);
	}

#else
	#ifdef REDUCE_WITH_CUB

			void *d_temp_storage = NULL;
			size_t temp_storage_bytes = 0;
			cub::DeviceReduce::Reduce(d_temp_storage, temp_storage_bytes,input,
																		output,
																		count,
																		op,
																		false);

			auto & temporal = context.getTemporalCUB();
			temporal.resize(temp_storage_bytes);

			// Run
			cub::DeviceReduce::Reduce(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,input,
					output,
					count,
					op,
					false);

	#else
			mgpu::reduce(input,count,output,op,context);
	#endif
#endif
	}
}

#endif

#endif /* REDUCE_OFP_HPP_ */