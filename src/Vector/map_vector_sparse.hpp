/*
 * map_vector_sparse.hpp
 *
 *  Created on: Jan 22, 2019
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_SPARSE_HPP_
#define MAP_VECTOR_SPARSE_HPP_

#include "Vector/map_vector.hpp"
#include "Vector/cuda/map_vector_sparse_cuda_ker.cuh"
#include "Vector/cuda/map_vector_sparse_cuda_kernels.cuh"
#include "util/cuda/ofp_context.hxx"
#include <iostream>

#ifdef __NVCC__
#include "util/cuda/moderngpu/kernel_scan.hxx"
#include "util/cuda/moderngpu/kernel_mergesort.hxx"
#include "util/cuda/moderngpu/kernel_segreduce.hxx"
#include "util/cuda/moderngpu/kernel_merge.hxx"
#include "util/cuda/kernels.cuh"
#endif

enum flush_type
{
	FLUSH_ON_HOST = 0,
	FLUSH_ON_DEVICE = 1
};

template<typename OfpmVectorT>
using ValueTypeOf = typename std::remove_reference<OfpmVectorT>::type::value_type;

namespace openfpm
{
    constexpr int VECTOR_SPARSE_STANDARD = 1;
    constexpr int VECTOR_SPARSE_BLOCK = 2;

    template<unsigned int impl, typename block_functor>
    struct scalar_block_implementation_switch // Case for scalar implementations
    {
        template <unsigned int p, typename vector_index_type>
        static void extendSegments(vector_index_type & segments, size_t dataSize)
        {
#ifdef __NVCC__
            // Pass as there is nothing to append for mgpu
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif // __NVCC__
        }

        template <typename vector_index_type>
        static void trimSegments(vector_index_type & segments)
        {
#ifdef __NVCC__
            // Pass as there is nothing to append for mgpu
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif // __NVCC__
        }

        template <unsigned int pSegment, typename vector_reduction, typename T, typename vector_data_type, typename vector_index_type>
        static void segreduce(vector_data_type & vector_data,
                vector_index_type & segment_offset,
                vector_data_type & vector_data_red,
                mgpu::ofp_context_t & context)
        {
#ifdef __NVCC__
            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<typename vector_data_type::value_type::type,typename reduction_type::prop>::type red_type;
            typedef typename reduction_type::op_red<red_type> red_op;
            red_type init;
            init = 0;

            mgpu::segreduce(
                    (red_type *)vector_data.template getDeviceBuffer<reduction_type::prop::value>(), vector_data.size(),
                    (int *)segment_offset.template getDeviceBuffer<1>(), segment_offset.size(),
                    (red_type *)vector_data_red.template getDeviceBuffer<reduction_type::prop::value>(),
                    red_op(), init, context);
#else // __NVCC__
    std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif // __NVCC__
        }

        template <
                typename vector_data_type,
                typename vector_index_type,
                typename vector_index_dtmp_type,
                typename dim3T,
                typename Ti,
                typename ... v_reduce>
        static void solveConflicts(
                vector_index_type & vct_index,
                vector_index_type & vct_index_tmp,
                vector_index_type & vct_index_tmp2,
                vector_index_type & vct_index_tmp3,
                vector_index_dtmp_type & vct_index_dtmp,
                vector_data_type & vct_data,
                vector_data_type & vct_add_data,
                vector_data_type & vct_add_data_unique,
                vector_data_type & vct_data_tmp,
                dim3T wthr,
                dim3T thr,
                mgpu::ofp_context_t & context
                )
        {
#ifdef __NVCC__
            solve_conflicts<
                        decltype(vct_index_tmp.toKernel()),
                        decltype(vct_data.toKernel()),
                        decltype(vct_index_dtmp.toKernel()),
                        128,
                        v_reduce ...
                        >
                <<<wthr,thr>>>
                                            (vct_index_tmp.toKernel(),vct_data.toKernel(),
                                              vct_index_tmp2.toKernel(),vct_add_data_unique.toKernel(),
                                              vct_index_tmp3.toKernel(),vct_data_tmp.toKernel(),
                                              vct_index_dtmp.toKernel(),
                                              vct_index.size());

                // we scan tmp3
                mgpu::scan(
                        (Ti*)vct_index_dtmp.template getDeviceBuffer<0>(),
                        vct_index_dtmp.size(),
                        (Ti *)vct_index_dtmp.template getDeviceBuffer<1>(),
                        context);

                // get the size to resize vct_index and vct_data
                vct_index_dtmp.template deviceToHost<0,1>(vct_index_dtmp.size()-1,vct_index_dtmp.size()-1);
                int size = vct_index_dtmp.template get<1>(vct_index_dtmp.size()-1) + vct_index_dtmp.template get<0>(vct_index_dtmp.size()-1);

                vct_index.resize(size);
                vct_data.resize(size);

                realign<<<wthr,thr>>>(vct_index_tmp3.toKernel(),vct_data_tmp.toKernel(),
                                      vct_index.toKernel(), vct_data.toKernel(),
                                      vct_index_dtmp.toKernel());
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif // __NVCC__
        }
    };

    template<typename block_functor>
    struct scalar_block_implementation_switch<2, block_functor> // Case for blocked implementations
    {
        template <unsigned int p, typename vector_index_type>
        static void extendSegments(vector_index_type & segments, size_t dataSize)
        {
#ifdef __NVCC__
            // Append trailing element to segment (marks end of last segment)
            segments.resize(segments.size()+1);
            segments.template get<p>(segments.size() - 1) = dataSize;
            segments.template hostToDevice<p>(segments.size() - 1, segments.size() - 1);
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif // __NVCC__
        }

        template <typename vector_index_type>
        static void trimSegments(vector_index_type & segments)
        {
#ifdef __NVCC__
            // Append trailing element to segment (marks end of last segment)
            segments.resize(segments.size()-1);
            ////////////
#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif // __NVCC__
        }

        template <unsigned int pSegment, typename vector_reduction, typename T, typename vector_data_type, typename vector_index_type>
        static void segreduce(vector_data_type & vector_data,
                       vector_index_type & segment_offset,
                       vector_data_type & vector_data_red,
                       mgpu::ofp_context_t & context)
        {
#ifdef __NVCC__
            block_functor::template seg_reduce<pSegment, vector_reduction, T>(segment_offset, vector_data, vector_data_red);
#else // __NVCC__
    std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif // __NVCC__
        }

        template <
                typename vector_data_type,
                typename vector_index_type,
                typename vector_index_dtmp_type,
                typename dim3T,
                typename Ti,
                typename ... v_reduce>
        static void solveConflicts(
                vector_index_type & vct_index,
                vector_index_type & vct_index_tmp,
                vector_index_type & vct_index_tmp2,
                vector_index_type & vct_index_tmp3,
                vector_index_dtmp_type & vct_index_dtmp,
                vector_data_type & vct_data,
                vector_data_type & vct_add_data,
                vector_data_type & vct_add_data_unique,
                vector_data_type & vct_data_tmp,
                dim3T wthr,
                dim3T thr,
                mgpu::ofp_context_t & context
        )
        {
#ifdef __NVCC__
            block_functor::template solve_conflicts<
			            decltype(vct_index_tmp),
			            decltype(vct_data),
			            v_reduce ...>
			            (vct_index_tmp, vct_index_tmp2,
			            vct_data, vct_add_data_unique,
			            vct_index_tmp3, vct_add_data,
			            vct_index, vct_data_tmp,
			            context);
                vct_data_tmp.swap(vct_data);

#else // __NVCC__
            std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif // __NVCC__
        }
    };

	template<typename Ti>
	struct reorder
	{
		Ti id;
		Ti id2;

		bool operator<(const reorder & t) const
		{
			return id < t.id;
		}
	};

	/*! \brief this class is a functor for "for_each" algorithm
	 *
	 * This class is a functor for "for_each" algorithm. For each
	 * element of the boost::vector the operator() is called.
	 * Is mainly used to copy one encap into another encap object
	 *
	 * \tparam encap source
	 * \tparam encap dst
	 *
	 */
	template<typename vector_data_type,
	        typename vector_index_type,
	        typename vector_reduction,
	        typename block_functor,
	        unsigned int impl2, unsigned int pSegment=1>
	struct sparse_vector_reduction
	{
		//! Vector in which to the reduction
		vector_data_type & vector_data_red;

		//! Vector in which to the reduction
		vector_data_type & vector_data;

		//! segment of offsets
		vector_index_type & segment_offset;

		//! gpu context
		mgpu::ofp_context_t & context;

		/*! \brief constructor
		 *
		 * \param src source encapsulated object
		 * \param dst source encapsulated object
		 *
		 */
		inline sparse_vector_reduction(vector_data_type & vector_data_red,
									   vector_data_type & vector_data,
									   vector_index_type & segment_offset,
									   mgpu::ofp_context_t & context)
		:vector_data_red(vector_data_red),vector_data(vector_data),segment_offset(segment_offset),context(context)
		{};

		//! It call the copy function for each property
		template<typename T>
		inline void operator()(T& t) const
		{
#ifdef __NVCC__

            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<typename ValueTypeOf<vector_data_type>::type,typename reduction_type::prop>::type red_type;
            if (reduction_type::is_special() == false)
			{
			    scalar_block_implementation_switch<impl2, block_functor>::template segreduce<pSegment, vector_reduction, T>(
			            vector_data,
			            segment_offset,
			            vector_data_red,
			            context);
			}
#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif
		}
	};

	struct stub_block_functor
	{
		template<typename vector_index_type, typename vector_data_type>
		static bool compact(vector_index_type & starts, int poolSize, vector_data_type & src, vector_data_type & dst)
		{
			return true;
		}

		template<typename vector_index_type, typename vector_data_type>
		static bool reorder(vector_index_type & src_id, vector_data_type & data, vector_data_type & data_reord)
		{
			return true;
		}

        template<unsigned int pSegment, typename vector_reduction, typename T, typename vector_index_type, typename vector_data_type>
		static bool seg_reduce(vector_index_type & segments, vector_data_type & src, vector_data_type &  dst)
		{
			return true;
		}

		template<typename vector_index_type, typename vector_data_type, typename ... v_reduce>
		static bool solve_conflicts(vector_index_type &keys, vector_index_type &merge_indices,
                                    vector_data_type &data1, vector_data_type &data2,
                                    vector_index_type &indices_tmp, vector_data_type &data_tmp,
                                    vector_index_type &keysOut, vector_data_type &dataOut,
                                    mgpu::ofp_context_t & context)
		{
			return true;
		}
	};

	/*! \brief this class is a functor for "for_each" algorithm
	 *
	 * This class is a functor for "for_each" algorithm. For each
	 * element of the boost::vector the operator() is called.
	 * Is mainly used to copy one encap into another encap object
	 *
	 * \tparam encap source
	 * \tparam encap dst
	 *
	 */
	template<typename vector_data_type, typename vector_index_type, typename vector_reduction>
	struct sparse_vector_special
	{
		//! Vector in which to the reduction
		vector_data_type & vector_data_red;

		//! Vector in which to the reduction
		vector_data_type & vector_data;

		//! segment of offsets
		vector_index_type & segment_offset;

		//! gpu context
		mgpu::ofp_context_t & context;

		/*! \brief constructor
		 *
		 * \param src source encapsulated object
		 * \param dst source encapsulated object
		 *
		 */
		inline sparse_vector_special(vector_data_type & vector_data_red,
									   vector_data_type & vector_data,
									   vector_index_type & segment_offset,
									   mgpu::ofp_context_t & context)
		:vector_data_red(vector_data_red),vector_data(vector_data),segment_offset(segment_offset),context(context)
		{};

		//! It call the copy function for each property
		template<typename T>
		inline void operator()(T& t) const
		{
#ifdef __NVCC__

			typedef typename boost::mpl::at<vector_reduction,T>::type reduction_type;

			// reduction type
			typedef typename boost::mpl::at<typename vector_data_type::value_type::type,typename reduction_type::prop>::type red_type;

			if (reduction_type::is_special() == true)
			{
				auto ite = segment_offset.getGPUIterator();

				reduce_from_offset<decltype(segment_offset.toKernel()),decltype(vector_data_red.toKernel()),reduction_type><<<ite.wthr,ite.thr>>>(segment_offset.toKernel(),vector_data_red.toKernel(),vector_data.size());
			}

#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: this file si supposed to be compiled with nvcc" << std::endl;
#endif
		}
	};

	template<typename T,
			 typename Ti = int,
			 typename Memory=HeapMemory,
			 typename layout=typename memory_traits_lin<T>::type,
			 template<typename> class layout_base=memory_traits_lin ,
			 typename grow_p=grow_policy_double,
			 unsigned int impl=vect_isel<T>::value,
			 unsigned int impl2 = VECTOR_SPARSE_STANDARD,
			 typename block_functor = stub_block_functor>
	class vector_sparse
	{
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_index;
		vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p,impl> vct_data;
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_m_index;

		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_add_index;
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_rem_index;
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_nadd_index;
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_nrem_index;
		vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p> vct_add_data;
		vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p> vct_add_data_reord;

		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_add_index_cont_0;
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_add_index_cont_1;
		vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p> vct_add_data_cont;
		vector<aggregate<Ti,Ti>,Memory,typename layout_base<aggregate<Ti,Ti>>::type,layout_base,grow_p> vct_add_index_unique;

		vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p,impl> vct_add_data_unique;

		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> starts;
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_index_tmp;
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_index_tmp2;
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_index_tmp3;
		vector<aggregate<Ti,Ti,Ti>,Memory,typename layout_base<aggregate<Ti,Ti,Ti>>::type,layout_base,grow_p> vct_index_dtmp;
		vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p,impl> vct_data_tmp;

		CudaMemory mem;

		openfpm::vector<reorder<Ti>> reorder_add_index_cpu;

		T bck;

		size_t max_ele;

		int n_gpu_add_block_slot;
		int n_gpu_rem_block_slot;

		/*! \brief get the element i
		 *
		 * search the element x
		 *
		 * \param i element i
		 */
		template<bool prefetch>
		Ti _branchfree_search(Ti x, Ti & id) const
		{
			if (vct_index.size() == 0)	{return -1;}
			const Ti *base = &vct_index.template get<0>(0);
			Ti n = vct_data.size();
			while (n > 1)
			{
				Ti half = n / 2;
				if (prefetch)
				{
					__builtin_prefetch(base + half/2, 0, 0);
					__builtin_prefetch(base + half + half/2, 0, 0);
				}
				base = (base[half] < x) ? base+half : base;
				n -= half;
			}

			int off = (*base < x);
			id = base - &vct_index.template get<0>(0) + off;
			return *(base + off);
		}

		template<typename ... v_reduce>
		void flush_on_gpu_insert(vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> & vct_add_index_cont_0,
				  vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> & vct_add_index_cont_1,
				  vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p> & vct_add_data_reord,
				  mgpu::ofp_context_t & context,
				  int i)
		{
#ifdef __NVCC__
			// Add 0 to the last element to vct_nadd_index
			vct_nadd_index.resize(vct_nadd_index.size()+1);
			vct_nadd_index.template get<0>(vct_nadd_index.size()-1) = 0;
			vct_nadd_index.template hostToDevice<0>(vct_nadd_index.size()-1,vct_nadd_index.size()-1);

			// Merge the list of inserted points for each block
			starts.resize(vct_nadd_index.size());

			mgpu::scan((Ti *)vct_nadd_index.template getDeviceBuffer<0>(),
			            vct_nadd_index.size(),
			            (Ti *)starts.template getDeviceBuffer<0>() ,
                        context);

			starts.template deviceToHost<0>(starts.size()-1,starts.size()-1);
			size_t n_ele = starts.template get<0>(starts.size()-1);

			// we reuse vct_nadd_index
			vct_add_index_cont_0.resize(n_ele);
			vct_add_index_cont_1.resize(n_ele);
			vct_add_data_cont.resize(n_ele);

			dim3 wthr;
			dim3 thr;
			wthr.x = vct_nadd_index.size()-1;
			wthr.y = 1;
			wthr.z = 1;
			thr.x = 128;
			thr.y = 1;
			thr.z = 1;

			if (impl2 == VECTOR_SPARSE_STANDARD)
			{
				construct_insert_list<<<wthr,thr>>>(vct_add_index.toKernel(),
										vct_nadd_index.toKernel(),
										starts.toKernel(),
										vct_add_index_cont_0.toKernel(),
										vct_add_index_cont_1.toKernel(),
										vct_add_data.toKernel(),
										vct_add_data_cont.toKernel(),
										n_gpu_add_block_slot);
			}
			else
			{
				construct_insert_list_key_only<<<wthr,thr>>>(vct_add_index.toKernel(),
										vct_nadd_index.toKernel(),
										starts.toKernel(),
										vct_add_index_cont_0.toKernel(),
										vct_add_index_cont_1.toKernel(),
										n_gpu_add_block_slot);

				block_functor::compact(starts, n_gpu_add_block_slot, vct_add_data, vct_add_data_cont);
			}

            // At this point we can check whether we have not inserted anything actually,
            // in this case, return without further ado...
			if (vct_add_data_cont.size() == 0)
            {
                return;
            }

			// now we sort
			mergesort(
			        (Ti *)vct_add_index_cont_0.template getDeviceBuffer<0>(),
                    (Ti *)vct_add_index_cont_1.template getDeviceBuffer<0>(),
					vct_add_index_cont_0.size(),
					mgpu::template less_t<Ti>(),
                    context);

			auto ite = vct_add_index_cont_0.getGPUIterator();

			vct_add_data_reord.resize(n_ele);
			// Now we reorder the data vector accordingly to the indexes

			if (impl2 == VECTOR_SPARSE_STANDARD)
			{
				reorder_vector_data<<<ite.wthr,ite.thr>>>(vct_add_index_cont_1.toKernel(),vct_add_data_cont.toKernel(),vct_add_data_reord.toKernel());
			}
			else
			{
				block_functor::reorder(vct_add_index_cont_1,vct_add_data_cont,vct_add_data_reord);
			}

			mem.allocate(sizeof(int));
			mem.fill(0);
			vct_add_index_unique.resize(vct_add_index_cont_0.size()+1);

			ite = vct_add_index_cont_0.getGPUIterator();

			// produce unique index list
			// Find the buffer bases
			CUDA_LAUNCH(
			        (
                        find_buffer_offsets_zero
                                <0,
                                decltype(vct_add_index_cont_0.toKernel()),
                                decltype(vct_add_index_unique.toKernel())
                                >
                    ),
                    ite,
                    vct_add_index_cont_0.toKernel(),
                    (int *)mem.getDevicePointer(),
                    vct_add_index_unique.toKernel());

			mem.deviceToHost();
			int n_ele_unique = *(int *)mem.getPointer();

			vct_add_index_unique.resize(n_ele_unique);
			vct_add_data_unique.resize(n_ele_unique);

			mgpu::mergesort(
			        (Ti *)vct_add_index_unique.template getDeviceBuffer<1>(),
                    (Ti *)vct_add_index_unique.template getDeviceBuffer<0>(),
                    vct_add_index_unique.size(),
                    mgpu::template less_t<Ti>(),
                    context);

			typedef boost::mpl::vector<v_reduce...> vv_reduce;

			// Now we can do a segmented reduction
			scalar_block_implementation_switch<impl2, block_functor>
			        ::template extendSegments<1>(vct_add_index_unique, vct_add_data_reord.size());

			sparse_vector_reduction<decltype(vct_add_data),decltype(vct_add_index_unique),vv_reduce,block_functor,impl2>
			        svr(
			                vct_add_data_unique,
			                vct_add_data_reord,
			                vct_add_index_unique,
			                context);
			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(svr);

			sparse_vector_special<decltype(vct_add_data),decltype(vct_add_index_unique),vv_reduce> svr2(vct_add_data_unique,vct_add_data_reord,vct_add_index_unique,context);
			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(svr2);

			scalar_block_implementation_switch<impl2, block_functor>
			        ::trimSegments(vct_add_index_unique);

			// Then we merge the two list vct_index and vct_add_index_unique

			// index to get merge index
			vct_m_index.resize(vct_index.size() + vct_add_index_unique.size());

			ite = vct_m_index.getGPUIterator();
			set_indexes<0><<<ite.wthr,ite.thr>>>(vct_m_index.toKernel(),0);

			ite = vct_add_index_unique.getGPUIterator();
			set_indexes<1><<<ite.wthr,ite.thr>>>(vct_add_index_unique.toKernel(),vct_index.size());

			// after merge we solve the last conflicts, running across the vector again and spitting 1 when there is something to merge
			// we reorder the data array also

			vct_index_tmp.resize(vct_index.size() + vct_add_index_unique.size());
			vct_index_tmp2.resize(vct_index.size() + vct_add_index_unique.size());
			vct_index_tmp3.resize(vct_index.size() + vct_add_index_unique.size());
			vct_data_tmp.resize(vct_index.size() + vct_add_index_unique.size());

			wthr.x = vct_index_tmp.size() / 128 + (vct_index_tmp.size() % 128 != 0);
			wthr.y = 1;
			wthr.z = 1;
			thr.x = 128;
			thr.y = 1;
			thr.z = 1;

			vct_index_dtmp.resize(wthr.x);

			// we merge with vct_index with vct_add_index_unique in vct_index_tmp, vct_intex_tmp2 contain the merging index
			//
			mgpu::merge((Ti *)vct_index.template getDeviceBuffer<0>(),(Ti *)vct_m_index.template getDeviceBuffer<0>(),vct_index.size(),
						(Ti *)vct_add_index_unique.template getDeviceBuffer<0>(),(Ti *)vct_add_index_unique.template getDeviceBuffer<1>(),vct_add_index_unique.size(),
						(Ti *)vct_index_tmp.template getDeviceBuffer<0>(),(Ti *)vct_index_tmp2.template getDeviceBuffer<0>(),mgpu::less_t<Ti>(),context);

            // Now perform the right solve_conflicts according to impl2
            scalar_block_implementation_switch<impl2, block_functor>::template solveConflicts<
                    decltype(vct_data),
                    decltype(vct_index),
                    decltype(vct_index_dtmp),
                    dim3,
                    Ti,
                    v_reduce ...
                    >
                    (
                        vct_index,
                        vct_index_tmp,
                        vct_index_tmp2,
                        vct_index_tmp3,
                        vct_index_dtmp,
                        vct_data,
                        vct_add_data,
                        vct_add_data_unique,
                        vct_data_tmp,
                        wthr,
                        thr,
                        context
                    );
#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif
		}


		void flush_on_gpu_remove(
				  mgpu::ofp_context_t & context)
		{
#ifdef __NVCC__

			// Add 0 to the last element to vct_nadd_index
			vct_nrem_index.resize(vct_nrem_index.size()+1);
			vct_nrem_index.template get<0>(vct_nrem_index.size()-1) = 0;
			vct_nrem_index.template hostToDevice<0>(vct_nrem_index.size()-1,vct_nrem_index.size()-1);

			// Merge the list of inserted points for each block
			starts.resize(vct_nrem_index.size());

			mgpu::scan((Ti *)vct_nrem_index.template getDeviceBuffer<0>(), vct_nrem_index.size(), (Ti *)starts.template getDeviceBuffer<0>() , context);

			starts.template deviceToHost<0>(starts.size()-1,starts.size()-1);
			size_t n_ele = starts.template get<0>(starts.size()-1);

			// we reuse vct_nadd_index
			vct_add_index_cont_0.resize(n_ele);
			vct_add_index_cont_1.resize(n_ele);

			dim3 wthr;
			dim3 thr;
			wthr.x = vct_nrem_index.size()-1;
			wthr.y = 1;
			wthr.z = 1;
			thr.x = 128;
			thr.y = 1;
			thr.z = 1;

			construct_remove_list<<<wthr,thr>>>(vct_rem_index.toKernel(),
										vct_nrem_index.toKernel(),
										starts.toKernel(),
										vct_add_index_cont_0.toKernel(),
										vct_add_index_cont_1.toKernel(),
										n_gpu_rem_block_slot);

			// now we sort
			mergesort((Ti *)vct_add_index_cont_0.template getDeviceBuffer<0>(),(Ti *)vct_add_index_cont_1.template getDeviceBuffer<0>(),
					vct_add_index_cont_0.size(), mgpu::template less_t<Ti>(), context);

			auto ite = vct_add_index_cont_0.getGPUIterator();

			mem.allocate(sizeof(int));
			mem.fill(0);
			vct_add_index_unique.resize(vct_add_index_cont_0.size()+1);

			ite = vct_add_index_cont_0.getGPUIterator();

			// produce unique index list
			// Find the buffer bases
			CUDA_LAUNCH((find_buffer_offsets_zero<0,decltype(vct_add_index_cont_0.toKernel()),decltype(vct_add_index_unique.toKernel())>),
					    ite,
					    vct_add_index_cont_0.toKernel(),(int *)mem.getDevicePointer(),vct_add_index_unique.toKernel());

			mem.deviceToHost();
			int n_ele_unique = *(int *)mem.getPointer();

			vct_add_index_unique.resize(n_ele_unique);

			mgpu::mergesort((Ti *)vct_add_index_unique.template getDeviceBuffer<1>(),(Ti *)vct_add_index_unique.template getDeviceBuffer<0>(),
							vct_add_index_unique.size(),mgpu::template less_t<Ti>(),context);

			// Then we merge the two list vct_index and vct_add_index_unique

			// index to get merge index
			vct_m_index.resize(vct_index.size() + vct_add_index_unique.size());

			ite = vct_m_index.getGPUIterator();
			set_indexes<0><<<ite.wthr,ite.thr>>>(vct_m_index.toKernel(),0);

			ite = vct_add_index_unique.getGPUIterator();
			set_indexes<1><<<ite.wthr,ite.thr>>>(vct_add_index_unique.toKernel(),vct_index.size());

			// after merge we solve the last conflicts, running across the vector again and spitting 1 when there is something to merge
			// we reorder the data array also

			vct_index_tmp.resize(vct_index.size() + vct_add_index_unique.size());
			vct_index_tmp2.resize(vct_index.size() + vct_add_index_unique.size());

			wthr.x = vct_index_tmp.size() / 128 + (vct_index_tmp.size() % 128 != 0);
			wthr.y = 1;
			wthr.z = 1;
			thr.x = 128;
			thr.y = 1;
			thr.z = 1;

			vct_index_dtmp.resize(wthr.x);

			// we merge with vct_index with vct_add_index_unique in vct_index_tmp, vct_intex_tmp2 contain the merging index
			//
			mgpu::merge((Ti *)vct_index.template getDeviceBuffer<0>(),(Ti *)vct_m_index.template getDeviceBuffer<0>(),vct_index.size(),
						(Ti *)vct_add_index_unique.template getDeviceBuffer<0>(),(Ti *)vct_add_index_unique.template getDeviceBuffer<1>(),vct_add_index_unique.size(),
						(Ti *)vct_index_tmp.template getDeviceBuffer<0>(),(Ti *)vct_index_tmp2.template getDeviceBuffer<0>(),mgpu::less_t<Ti>(),context);

			vct_index_tmp3.resize(128*wthr.x);

			solve_conflicts_remove<decltype(vct_index_tmp.toKernel()),decltype(vct_index_dtmp.toKernel()),128>
			<<<wthr,thr>>>
										(vct_index_tmp.toKernel(),
										 vct_index_tmp2.toKernel(),
										 vct_index_tmp3.toKernel(),
										 vct_m_index.toKernel(),
										  vct_index_dtmp.toKernel(),
										  vct_index.size());

			// we scan tmp3
			mgpu::scan((Ti*)vct_index_dtmp.template getDeviceBuffer<0>(),vct_index_dtmp.size(),(Ti *)vct_index_dtmp.template getDeviceBuffer<1>(),context);

			// get the size to resize vct_index and vct_data
			vct_index_dtmp.template deviceToHost<0,1>(vct_index_dtmp.size()-1,vct_index_dtmp.size()-1);
			int size = vct_index_dtmp.template get<1>(vct_index_dtmp.size()-1) + vct_index_dtmp.template get<0>(vct_index_dtmp.size()-1);

			vct_data_tmp.resize(size);
			vct_index.resize(size);

			realign_remove<<<wthr,thr>>>(vct_index_tmp3.toKernel(),vct_m_index.toKernel(),vct_data.toKernel(),
								  vct_index.toKernel(),vct_data_tmp.toKernel(),
								  vct_index_dtmp.toKernel());

			vct_data.swap(vct_data_tmp);

#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: you are suppose to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif
		}

		template<typename ... v_reduce>
		void flush_on_gpu(vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> & vct_add_index_cont_0,
						  vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> & vct_add_index_cont_1,
						  vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p> & vct_add_data_reord,
						  mgpu::ofp_context_t & context,
						  int i)
		{
			flush_on_gpu_insert<v_reduce ... >(vct_add_index_cont_0,vct_add_index_cont_1,vct_add_data_reord,context,i);
		}

		void flush_on_cpu()
		{
			// First copy the added index to reorder
			reorder_add_index_cpu.resize(vct_add_index.size());

			for (size_t i = 0 ; i < reorder_add_index_cpu.size() ; i++)
			{
				reorder_add_index_cpu.get(i).id = vct_add_index.template get<0>(i);
				reorder_add_index_cpu.get(i).id2 = i;
			}

			reorder_add_index_cpu.sort();

			// merge the the data

			vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p,impl> vct_data_tmp;
			vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> vct_index_tmp;

			vct_data_tmp.resize(vct_data.size() + vct_add_data.size());
			vct_index_tmp.resize(vct_index.size() + vct_add_index.size());

			Ti di = 0;
			Ti ai = 0;
			size_t i = 0;

			for ( ; i < vct_data_tmp.size() && ai < vct_add_index.size() && di < vct_index.size() ; i++)
			{
				Ti id_a = reorder_add_index_cpu.get(ai).id;
				Ti id_d = vct_index.template get<0>(di);

				if (  id_a <= id_d )
				{
					Ti pos = reorder_add_index_cpu.get(ai).id2;

					vct_index_tmp.template get<0>(i) = reorder_add_index_cpu.get(ai).id;
					vct_data_tmp.get(i) = vct_add_data.get(pos);

					if (id_a == id_d)
					{di++;}

					// duplicate, set again
					if (ai+1 < reorder_add_index_cpu.size() && reorder_add_index_cpu.get(ai+1).id == id_a)
					{i--;}

					ai++;
				}
				else
				{
					vct_index_tmp.template get<0>(i) = vct_index.template get<0>(di);
					vct_data_tmp.get(i) = vct_data.get(di);
					di++;
				}
			}

			for ( ; ai < vct_add_index.size() ; ai++, i++)
			{
				Ti pos = reorder_add_index_cpu.get(ai).id2;

				vct_index_tmp.template get<0>(i) = reorder_add_index_cpu.get(ai).id;
				vct_data_tmp.get(i) = vct_add_data.get(pos);
			}

			for ( ; di < vct_index.size() ; di++ , i++)
			{
				vct_index_tmp.template get<0>(i) = vct_index.template get<0>(di);
				vct_data_tmp.get(i) = vct_data.get(di);
			}

			// if there are duplicate vct_index_tmp can store less entry
			vct_index_tmp.resize(i);
			vct_data_tmp.resize(i);

			vct_index.swap(vct_index_tmp);
			vct_data.swap(vct_data_tmp);

			vct_add_data.clear();
			vct_add_index.clear();
		}

	public:

		vector_sparse()
		:max_ele(0)
		{
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \tparam p Property to get
		 * \param id Element to get
		 *
		 * \return the element value requested
		 *
		 */
		template <unsigned int p>
		inline auto get(Ti id) const -> decltype(vct_data.template get<p>(id))
		{
			Ti di;
			Ti v = this->_branchfree_search<false>(id,di);
			return (v == id)?vct_data.template get<p>(di):bck.template get<p>();
		}

		/*! \brief resize to n elements
		 *
		 * \param n elements
		 *
		 */
		void resize(size_t n)
		{
			max_ele = n;
		}

		/*! \brief
		 *
		 * \warning After using this function to move out the vector of the indexes, this object become useless and
		 *          must be destroyed
		 *
		 * \param iv
		 *
		 */
		void swapIndexVector(vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> & iv)
		{
			vct_index.swap(iv);
		}

		/*! \brief Set the background to bck (which value get must return when the value is not find)
		 *
		 * \param bck
		 *
		 */
		template <unsigned int p>
		auto getBackground() -> decltype(bck.template get<p>())
		{
			return bck.template get<p>();
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		template <unsigned int p>
		auto insert(Ti ele) -> decltype(vct_data.template get<p>(0))
		{
			vct_add_index.add();
			vct_add_index.template get<0>(vct_add_index.size()-1) = ele;
			vct_add_data.add();
			return vct_add_data.template get<p>(vct_add_data.size()-1);
		}

		/*! \brief merge the added element to the main data array but save the insert buffer in v
		 *
		 * \param v insert buffer
		 *
		 * \param opt options
		 *
		 */
		template<typename ... v_reduce>
		void flush_v(vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> & vct_add_index_cont_0,
				     mgpu::ofp_context_t & context,
				     flush_type opt = FLUSH_ON_HOST,
				     int i = 0)
		{
			if (opt & flush_type::FLUSH_ON_DEVICE)
			{this->flush_on_gpu<v_reduce ... >(vct_add_index_cont_0,vct_add_index_cont_1,vct_add_data_reord,context,i);}
			else
			{this->flush_on_cpu();}
		}

		/*! \brief merge the added element to the main data array but save the insert buffer in v
		 *
		 * \param v insert buffer
		 *
		 * \param opt options
		 *
		 */
		template<typename ... v_reduce>
		void flush_vd(vector<T,Memory,typename layout_base<T>::type,layout_base,grow_p> & vct_add_data_reord,
				     mgpu::ofp_context_t & context,
				     flush_type opt = FLUSH_ON_HOST,
				     int i = 0)
		{
			if (opt & flush_type::FLUSH_ON_DEVICE)
			{this->flush_on_gpu<v_reduce ... >(vct_add_index_cont_0,vct_add_index_cont_1,vct_add_data_reord,context,i);}
			else
			{this->flush_on_cpu();}
		}

		/*! \brief merge the added element to the main data array
		 *
		 * \param opt options
		 *
		 */
		template<typename ... v_reduce>
		void flush(mgpu::ofp_context_t & context, flush_type opt = FLUSH_ON_HOST, int i = 0)
		{
			if (opt & flush_type::FLUSH_ON_DEVICE)
			{this->flush_on_gpu<v_reduce ... >(vct_add_index_cont_0,vct_add_index_cont_1,vct_add_data_reord,context,i);}
			else
			{this->flush_on_cpu();}
		}

		/*! \brief merge the added element to the main data array
		 *
		 * \param opt options
		 *
		 */
		void flush_remove(mgpu::ofp_context_t & context, flush_type opt = FLUSH_ON_HOST)
		{
			if (opt & flush_type::FLUSH_ON_DEVICE)
			{this->flush_on_gpu_remove(context);}
			else
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " error, flush_remove on CPU has not implemented yet";
			}
		}

		/*! \brief Return how many element you have in this map
		 *
		 * \return the number of elements
		 */
		size_t size()
		{
			return vct_index.size();
		}

		/*! \brief Return the sorted vector of the indexes
		 *
		 * \return return the sorted vector of the indexes
		 */
		vector<aggregate<Ti>,Memory,typename layout_base<aggregate<Ti>>::type,layout_base,grow_p> &
		private_get_vct_index()
		{
			return vct_index;
		}

		/*! \brief Transfer from device to host
		 *
		 * \tparam set of parameters to transfer to host
		 *
		 */
		template<unsigned int ... prp>
		void deviceToHost()
		{
			vct_index.template deviceToHost<0>();
			vct_data.template deviceToHost<prp...>();
		}

		/*! \brief toKernel function transform this structure into one that can be used on GPU
		 *
		 * \return structure that can be used on GPU
		 *
		 */
		vector_sparse_gpu_ker<T,Ti,layout_base> toKernel()
		{
			vector_sparse_gpu_ker<T,Ti,layout_base> mvsck(vct_index.toKernel(),vct_data.toKernel(),
														  vct_add_index.toKernel(),
														  vct_rem_index.toKernel(),vct_add_data.toKernel(),
														  vct_nadd_index.toKernel(),
														  vct_nrem_index.toKernel(),
														  bck,
														  n_gpu_add_block_slot,
														  n_gpu_rem_block_slot);

			return mvsck;
		}

		/*! \brief set the gpu insert buffer for every block
		 *
		 * \param nblock number of blocks
		 * \param nslot number of slots free for each block
		 *
		 */
		void setGPUInsertBuffer(int nblock, int nslot)
		{
			vct_add_index.resize(nblock*nslot);
			vct_nadd_index.resize(nblock);
			vct_add_data.resize(nblock*nslot);
			n_gpu_add_block_slot = nslot;
			vct_nadd_index.template fill<0>(0);
		}

        /*! \brief Get the GPU insert buffer
         *
         * \return the reference to the GPU insert buffer
         */
        auto getGPUInsertBuffer() -> decltype(vct_add_data)&
        {
            return vct_add_data;
        }

		/*! \brief set the gpu remove buffer for every block
		 *
		 * \param nblock number of blocks
		 * \param nslot number of slots free for each block
		 *
		 */
		void setGPURemoveBuffer(int nblock, int nslot)
		{
			vct_rem_index.resize(nblock*nslot);
			vct_nrem_index.resize(nblock);
			n_gpu_rem_block_slot = nslot;
			vct_nrem_index.template fill<0>(0);
		}

#ifdef CUDA_GPU

		/*! \brief Get iterator over the stored elements
		 *
		 * \return an iterator
		 *
		 */
		auto getGPUIterator() -> decltype(vct_index.getGPUIterator())
		{
			return vct_index.getGPUIterator();
		}

#endif

		/*! \brief Clear all from all the elements
		 *
		 *
		 */
		void clear()
		{
			vct_data.clear();
			vct_index.clear();
			vct_add_index.clear();
			vct_add_data.clear();

			max_ele = 0;
			n_gpu_add_block_slot = 0;
			n_gpu_rem_block_slot = 0;
		}

		void swap(vector_sparse<T,Ti,Memory,layout,layout_base,grow_p,impl> & sp)
		{
			vct_data.swap(sp.vct_data);
			vct_index.swap(sp.vct_index);
			vct_add_index.swap(sp.vct_add_index);
			vct_add_data.swap(sp.vct_add_data);

			size_t max_ele_ = sp.max_ele;
			sp.max_ele = max_ele;
			this->max_ele = max_ele_;

			decltype(bck) bck_tmp = sp.bck;
			sp.bck = bck;
			this->bck = bck_tmp;
		}
	};

	template<typename T, unsigned int blockSwitch = VECTOR_SPARSE_STANDARD, typename block_functor = stub_block_functor>
	using vector_sparse_gpu = openfpm::vector_sparse<
	        T,
	        int,
	        CudaMemory,
	        typename memory_traits_inte<T>::type,
	        memory_traits_inte,
            grow_policy_double,
            vect_isel<T>::value,
            blockSwitch,
            block_functor
            >;
    template<typename T, unsigned int blockSwitch = VECTOR_SPARSE_STANDARD, typename block_functor = stub_block_functor>
    using vector_sparse_u_gpu = openfpm::vector_sparse<
            T,
            unsigned int,
            CudaMemory,
            typename memory_traits_inte<T>::type,
            memory_traits_inte,
            grow_policy_double,
            vect_isel<T>::value,
            blockSwitch,
            block_functor
    >;
}



#endif /* MAP_VECTOR_SPARSE_HPP_ */
