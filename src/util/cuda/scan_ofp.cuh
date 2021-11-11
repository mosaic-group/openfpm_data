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
	#include "cub/cub.cuh"
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

    typename std::remove_reference<decltype(input[0])>::type scan_a = 0;
	// No particular speed-up from omp amd simd
//    #pragma omp parallel for simd reduction(inscan, +:scan_a)
    for(int i = 0; i < count; i++)
	{
        output[i] = scan_a;
//        #pragma omp scan exclusive(scan_a)
        scan_a += input[i];
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

#endif /* __NVCC__ */

#endif /* SCAN_OFP_HPP_ */
