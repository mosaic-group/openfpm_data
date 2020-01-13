/*
 * sort_ofp.cuh
 *
 *  Created on: Aug 23, 2019
 *      Author: i-bird
 */

#ifndef SORT_OFP_CUH_
#define SORT_OFP_CUH_

#ifdef __HIPCC__
// On hip only CUB is activated
#define SORT_WITH_CUB
#endif

#if defined(__NVCC__) || defined(__HIPCC__)

#include "hipcub/hipcub.hpp"
#ifdef __NVCC__
#include "util/cuda/moderngpu/kernel_scan.hxx"
#endif
#ifdef __HIPCC__
#include "util/cuda/modern_gpu_hipcc.hpp"
#endif
#include "ofp_context.cuh"

namespace openfpm
{
	template<typename key_t, typename val_t,
	  typename comp_t>
	void sort(key_t* keys_input, val_t* vals_input, int count,
	  comp_t comp, mgpu::ofp_context_t& context)
	{
#ifdef SORT_WITH_CUB

		if (std::is_same<mgpu::template less_t<key_t>,comp_t>::value == true)
		{
			void *d_temp_storage = NULL;
			size_t temp_storage_bytes = 0;
			#ifdef  __NVCC__
			cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
			#else
			hipcub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
			#endif

			auto & temporal = context.getTemporalCUB();
			temporal.resize(temp_storage_bytes);

			// Run
			#ifdef  __NVCC__
			cub::DeviceRadixSort::SortPairs(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
			#else
			hipcub::DeviceRadixSort::SortPairs(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
			#endif
		}
		else if (std::is_same<mgpu::template greater_t<key_t>,comp_t>::value == true)
		{
			void *d_temp_storage = NULL;
			size_t temp_storage_bytes = 0;
			#ifdef  __NVCC__
			cub::DeviceRadixSort::SortPairsDescending(d_temp_storage, temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
			#else
			hipcub::DeviceRadixSort::SortPairsDescending(d_temp_storage, temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
			#endif

			auto & temporal = context.getTemporalCUB();
			temporal.resize(temp_storage_bytes);

			// Run
			#ifdef  __NVCC__
			cub::DeviceRadixSort::SortPairsDescending(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
			#else
			hipcub::DeviceRadixSort::SortPairsDescending(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
			#endif
		}

#elif !defined(__HIPIFY__)
		mgpu::mergesort(keys_input,vals_input,count,comp,context);
#endif
	}
}

#endif


#endif /* SORT_OFP_CUH_ */
