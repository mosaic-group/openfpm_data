/*
 * segreduce_ofp.hpp
 *
 *  Created on: May 15, 2019
 *      Author: i-bird
 */

 #ifndef SEGREDUCE_OFP_HPP_
 #define SEGREDUCE_OFP_HPP_
 
 #ifdef __NVCC__
 
 #include "util/cuda_util.hpp"
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
    template<typename input_it,
             typename segments_it, typename output_it, typename op_t, typename type_t>
    void segreduce(input_it input, int count, segments_it segments,
                    int num_segments, output_it output, op_t op, type_t init,
                    gpu::ofp_context_t& gpuContext)
     {
 #ifdef CUDA_ON_CPU
 
        int i = 0;
        for ( ; i < num_segments - 1; i++)
        {
            int j = segments[i];
            output[i] = init;
            if (j == segments[i+1]) {continue;}
            output[i] = input[j];
            ++j;
            for ( ; j < segments[i+1] ; j++)
            {
                output[i] = op(output[i],input[j]);
            }
        }

        // Last segment
        int j = segments[i];
        if (j != count)
        {
            output[i] = input[j];
            ++j;
            for ( ; j < count ; j++)
            {
                output[i] = op(output[i],input[j]);
            }
        }
 
 #else
        #ifdef __HIP__

            size_t temp_storage_bytes = 0;

            hipcub::DeviceSegmentedReduce::Reduce(NULL, temp_storage_bytes, input, output,
                num_segments, segments, segments + 1, op, init);

            auto & temporal = gpuContext.getTemporalCUB();
            temporal.resize(temp_storage_bytes);

            hipcub::DeviceSegmentedReduce::Reduce(temporal.getDeviceBuffer<0>(), temp_storage_bytes, input, output,
                num_segments, segments, segments + 1, op, init);

        #else

            size_t temp_storage_bytes = 0;

            cub::DeviceSegmentedReduce::Reduce(NULL, temp_storage_bytes, input, output,
                num_segments, segments, segments + 1, op, init);

            auto & temporal = gpuContext.getTemporalCUB();
            temporal.resize(temp_storage_bytes);

            cub::DeviceSegmentedReduce::Reduce(temporal.template getDeviceBuffer<0>(), temp_storage_bytes, input, output,
                num_segments, segments, segments + 1, op, init);

        #endif
 #endif
     }
 }
 
 #endif /* __NVCC__ */
 
 #endif /* SCAN_OFP_HPP_ */
 
