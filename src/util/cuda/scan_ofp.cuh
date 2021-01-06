/*
 * scan_ofp.hpp
 *
 *  Created on: May 15, 2019
 *      Author: i-bird
 */

#ifndef SCAN_OFP_HPP_
#define SCAN_OFP_HPP_

#ifdef __NVCC__

#if defined(CUDA_ON_CPU)
#include "util/cuda_launch.hpp"
#endif

#if CUDART_VERSION < 11000
#include "cub_old/cub.cuh"
#include "util/cuda/moderngpu/kernel_scan.hxx"
#else
#include "cub/cub.cuh"
#ifndef SCAN_WITH_CUB
#define SCAN_WITH_CUB
#endif
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

	#else
			mgpu::scan(input,count,output,context);
	#endif
#endif
	}
}

#endif

#endif /* SCAN_OFP_HPP_ */
