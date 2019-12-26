/*
 * sort_ofp.cuh
 *
 *  Created on: Aug 23, 2019
 *      Author: i-bird
 */

#ifndef SORT_OFP_CUH_
#define SORT_OFP_CUH_


#if defined(__NVCC__) || defined(__HIPCC__)

#include "cub/cub.cuh"
#include "util/cuda/moderngpu/kernel_mergesort.hxx"
#include "util/cuda/ofp_context.hxx"

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
			cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);

			auto & temporal = context.getTemporalCUB();
			temporal.resize(temp_storage_bytes);

			// Run
			cub::DeviceRadixSort::SortPairs(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
		}
		else if (std::is_same<mgpu::template greater_t<key_t>,comp_t>::value == true)
		{
			void *d_temp_storage = NULL;
			size_t temp_storage_bytes = 0;
			cub::DeviceRadixSort::SortPairsDescending(d_temp_storage, temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);

			auto & temporal = context.getTemporalCUB();
			temporal.resize(temp_storage_bytes);

			// Run
			cub::DeviceRadixSort::SortPairsDescending(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,keys_input,keys_input,vals_input,vals_input,count);
		}

#else
		mgpu::mergesort(keys_input,vals_input,count,comp,context);
#endif
	}
}

#endif


#endif /* SORT_OFP_CUH_ */
