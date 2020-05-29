/*
 * scan_ofp.hpp
 *
 *  Created on: May 15, 2019
 *      Author: i-bird
 */

#ifndef SCAN_OFP_HPP_
#define SCAN_OFP_HPP_

#if defined(__HIPCC__) || defined(__HIPIFY__)
// On hip only CUB is activated
#define SCAN_WITH_CUB
#endif

#if defined(__NVCC__) || defined(__HIPCC__)

#if defined(__HIPCC__)
#include "hipcub/hipcub.hpp"
#else
#include "cub/cub.cuh"
#endif

#ifdef __NVCC__
#include "util/cuda/moderngpu/kernel_scan.hxx"
#endif
#include "ofp_context.cuh"

namespace openfpm
{
#ifdef __HIPCC__
        template<typename input_it, typename output_it>
        void scan(input_it input, int count, output_it output, mgpu::ofp_context_t& context)
#else
        template<mgpu::scan_type_t scan_type = mgpu::scan_type_exc,
                        typename launch_arg_t = mgpu::empty_t,
                        typename input_it, typename output_it>
                        void scan(input_it input, int count, output_it output, mgpu::ofp_context_t& context)
#endif
	{
#ifdef SCAN_WITH_CUB

	    void *d_temp_storage = NULL;
	    size_t temp_storage_bytes = 0;
	    #ifdef __NVCC__
	    cub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes,input,
	    		                                                    output,
	    		                                                    count);
	    #else
            hipcub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes,input,
                                                                            output,
                                                                            count);
	    #endif

	    auto & temporal = context.getTemporalCUB();
	    temporal.resize(temp_storage_bytes);

	    // Run
	    #ifdef __NVCC__
	    cub::DeviceScan::ExclusiveSum(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,input,
	            output,
	            count);
	    #else
            hipcub::DeviceScan::ExclusiveSum(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,input,
                    output,
                    count);
	    #endif

#else
		mgpu::scan(input,count,output,context);
#endif
	}
}

#endif

#endif /* SCAN_OFP_HPP_ */
