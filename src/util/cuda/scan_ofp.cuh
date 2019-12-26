/*
 * scan_ofp.hpp
 *
 *  Created on: May 15, 2019
 *      Author: i-bird
 */

#ifndef SCAN_OFP_HPP_
#define SCAN_OFP_HPP_

#if defined(__NVCC__) || defined(__HIPCC__)

#include "cub/cub.cuh"
#include "util/cuda/moderngpu/kernel_scan.hxx"
#include "util/cuda/ofp_context.hxx"

namespace openfpm
{
	template<mgpu::scan_type_t scan_type = mgpu::scan_type_exc,
			typename launch_arg_t = mgpu::empty_t,
			typename input_it, typename output_it>
			void scan(input_it input, int count, output_it output, mgpu::ofp_context_t& context)
	{
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
	}
}

#endif

#endif /* SCAN_OFP_HPP_ */
