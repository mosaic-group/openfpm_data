/*
 * segreduce_ofp.hpp
 *
 *  Created on: May 15, 2019
 *      Author: i-bird
 */

 #ifndef SEGREDUCE_OFP_HPP_
 #define SEGREDUCE_OFP_HPP_
 
 #ifdef __NVCC__
 
 #include "Vector/map_vector.hpp"
 #include "util/cuda_launch.hpp"
 #include "util/cuda/segreduce_ofp.cuh"
 
 #if CUDART_VERSION >= 11000
     #ifndef CUDA_ON_CPU 
     // Here we have for sure CUDA >= 11
     #ifdef __HIP__
        #undef __CUDACC__
        #undef __CUDA__
        #include <thrust/reduce.h>
        #define __CUDACC__
        #define __CUDA__
     #else
        #include "util/cuda/moderngpu/kernel_segreduce.hxx"
     #endif
     #endif
 #else
    #include "util/cuda/moderngpu/kernel_segreduce.hxx"
 #endif
 #include "util/cuda/ofp_context.hxx"
 
template<typename segments_it, typename keys_type, typename output_it, typename seg_type, typename type_t>
__global__ void seg_to_keys(segments_it segs, keys_type keys, seg_type seg_out ,output_it output, int n_count, int num_segments,type_t init)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= num_segments)    {return;}

    int s = segs[tid];
    int s_p1 = (tid == num_segments -1)?n_count:segs[tid+1];

    int n_ele = s_p1 - s;

    seg_out.template get<1>(tid) = (s != s_p1);
    output[tid] = init;

    for (int j = 0 ; j < n_ele ; j++)
    {
        keys.template get<0>(s + j) = tid;
    }
}

template<typename output_it, typename out_tmp_type ,typename segs_type>
__global__ void realign_output(output_it out, out_tmp_type out_tmp, segs_type segs, int num_segments)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= num_segments)    {return;}

    int t = segs.template get<2>(tid);
    int to_copy = segs.template get<1>(tid);

    auto op = out_tmp.template get<0>(t);

    if (to_copy == 1)
    {out[tid] = op;}
}

 namespace openfpm
 {
    template<typename input_it,
             typename segments_it, typename output_it, typename op_t, typename type_t>
    void segreduce(input_it input, int count, segments_it segments,
                    int num_segments, output_it output, op_t op, type_t init,
                    mgpu::ofp_context_t & context)
     {
 #ifdef CUDA_ON_CPU
 
        int i = 0;
        for ( ; i < num_segments - 1; i++)
        {
            int j = segments[i];
            output[i] = input[j];
            ++j;
            for ( ; j < segments[i+1] ; j++)
            {
                output[i] = op(output[i],input[j]);
            }
        }

        // Last segment
        int j = segments[i];
        output[i] = input[j];
        ++j;
        for ( ; j < count ; j++)
        {
            output[i] = op(output[i],input[j]);
        }
 
 #else

        #ifdef __HIP__

            typedef typename std::remove_pointer<segments_it>::type index_type;
            typedef typename std::remove_pointer<output_it>::type out_type;

            openfpm::vector_gpu<aggregate<index_type>> keys;
            keys.resize(count);

            openfpm::vector_gpu<aggregate<index_type,index_type,index_type>> segs_out;
            segs_out.resize(num_segments);

            openfpm::vector_gpu<aggregate<out_type>> out_tmp;
            out_tmp.resize(num_segments);

            grid_sm<1,void> g(num_segments);

            auto it = g.getGPUIterator();

            CUDA_LAUNCH(seg_to_keys,it,segments,keys.toKernel(),segs_out.toKernel(),output,count,num_segments,init);

            openfpm::scan((index_type *)segs_out.template getDeviceBuffer<1>(),num_segments,(index_type *)segs_out.template getDeviceBuffer<2>(),context);

            thrust::pair<index_type *,out_type *> new_end;
            new_end = thrust::reduce_by_key(thrust::device, (segments_it)keys.template getDeviceBuffer<0>(),((segments_it)keys.template getDeviceBuffer<0>()) + count, 
                                            input, 
                                            (segments_it)segs_out.template getDeviceBuffer<0>(), 
                                            (output_it)out_tmp.template getDeviceBuffer<0>(),
                                            thrust::equal_to<int>(),
                                            op);

            // .. Not so easy to emulate a segmented reduce we have to track the zeros segments and realign the output

            CUDA_LAUNCH(realign_output,it,output,out_tmp.toKernel(),segs_out.toKernel(),num_segments);

        #else

            mgpu::segreduce(input,count,segments,num_segments,output,op,init,context);

        #endif

 #endif
     }
 }
 
 #endif /* __NVCC__ */
 
 #endif /* SCAN_OFP_HPP_ */
 