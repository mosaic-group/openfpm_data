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

#if CUDART_VERSION >= 11000
	#ifndef CUDA_ON_CPU 
	// Here we have for sure CUDA >= 11
	#ifdef __HIP__
		#include "hipcub/hipcub.hpp"
	#else
		#include "cub/cub.cuh"
	#endif
	#ifndef SCAN_WITH_CUB
		#define SCAN_WITH_CUB
	#endif
	#endif
#else
	// Here we have old CUDA
	#include "cub_old/cub.cuh"
	#include "util/cuda/moderngpu/kernel_scan.hxx"
#endif

#include "util/cuda/ofp_context.hxx"

namespace openfpm
{
	template<typename input_it, typename output_it>
			 void scan(input_it input, int count, output_it output, mgpu::ofp_context_t& context)
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
	#ifdef SCAN_WITH_CUB

			#ifdef __HIP__

				if (count == 0)	{return;}

				void *d_temp_storage = NULL;
				size_t temp_storage_bytes = 0;
				hipcub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes,input,
																			output,
																			count);

				auto & temporal = context.getTemporalCUB();
				temporal.resize(temp_storage_bytes);

				// Run
				hipcub::DeviceScan::ExclusiveSum(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,input,
						output,
						count);

			#else

				void *d_temp_storage = NULL;
				size_t temp_storage_bytes = 0;
				cub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes,input,
																			output,
																			count);

				auto & temporal = context.getTemporalCUB();
				temporal.resize(temp_storage_bytes);

				// Run
				cub::DeviceScan::ExclusiveSum(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,input,
						output,
						count);

			#endif

	#else
			mgpu::scan(input,count,output,context);
	#endif
#endif
	}
}

#endif /* __NVCC__ */

#endif /* SCAN_OFP_HPP_ */
