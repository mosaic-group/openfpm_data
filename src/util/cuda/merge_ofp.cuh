/*
 * segreduce_ofp.hpp
 *
 *  Created on: May 15, 2019
 *      Author: i-bird
 */

 #ifndef MERGE_OFP_HPP_
 #define MERGE_OFP_HPP_
 
 #ifdef __NVCC__
 
 #include "Vector/map_vector.hpp"
 #include "util/cuda_util.hpp"
 
 #ifndef CUDA_ON_CPU
     // Here we have for sure CUDA >= 11
     #ifdef __HIP__
        #undef __CUDACC__
        #undef __CUDA__
        #include <thrust/merge.h>
        #include <thrust/execution_policy.h>
        #define __CUDACC__
        #define __CUDA__
     #else
        #include <thrust/merge.h>
        #include <thrust/execution_policy.h>
     #endif
 #endif
 

 namespace openfpm
 {
    template<typename a_keys_it, typename a_vals_it,
             typename b_keys_it, typename b_vals_it,
             typename c_keys_it, typename c_vals_it,
             typename comp_t, typename context_t>
    void merge(a_keys_it a_keys, a_vals_it a_vals, int a_count,
               b_keys_it b_keys, b_vals_it b_vals, int b_count,
            c_keys_it c_keys, c_vals_it c_vals, comp_t comp, context_t& gpuContext)
    {
 #ifdef CUDA_ON_CPU
 
        int a_it = 0;
        int b_it = 0;
        int c_it = 0;

        while (a_it < a_count || b_it < b_count)
        {
            if (a_it < a_count)
            {
                if (b_it < b_count)
                {
                    if (comp(b_keys[b_it],a_keys[a_it]))
                    {
                        c_keys[c_it] = b_keys[b_it];
                        c_vals[c_it] = b_vals[b_it];
                        c_it++;
                        b_it++;
                    }
                    else
                    {
                        c_keys[c_it] = a_keys[a_it];
                        c_vals[c_it] = a_vals[a_it];
                        c_it++;
                        a_it++;
                    }
                }
                else
                {
                    c_keys[c_it] = a_keys[a_it];
                    c_vals[c_it] = a_vals[a_it];
                    c_it++;
                    a_it++;
                }
            }
            else
            {
                c_keys[c_it] = b_keys[b_it];
                c_vals[c_it] = b_vals[b_it];
                c_it++;
                b_it++;
            }
        }
 
 #else

        #ifdef __HIP__

            thrust::merge_by_key(thrust::device, a_keys,a_keys + a_count, 
                                                 b_keys,b_keys + b_count, 
                                                 a_vals,b_vals,
                                                 c_keys,c_vals,comp);

        #else

            thrust::merge_by_key(thrust::device, a_keys,a_keys + a_count, 
                                                 b_keys,b_keys + b_count, 
                                                 a_vals,b_vals,
                                                 c_keys,c_vals,comp);

        #endif

 #endif
    }
 }
 
 #endif /* __NVCC__ */
 
 #endif /* SCAN_OFP_HPP_ */
