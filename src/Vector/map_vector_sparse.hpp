/*
 * map_vector_sparse.hpp
 *
 *  Created on: Jan 22, 2019
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_SPARSE_HPP_
#define MAP_VECTOR_SPARSE_HPP_

#include "util/cuda_launch.hpp"
#include "Vector/map_vector.hpp"
#include "Vector/cuda/map_vector_sparse_cuda_ker.cuh"
#include "Vector/cuda/map_vector_sparse_cuda_kernels.cuh"
#include "util/ofp_context.hpp"
#include <iostream>
#include <limits>

#if defined(__NVCC__)
 #include "util/cuda/kernels.cuh"
#endif

#include "util/cuda/scan_ofp.cuh"
#include "util/cuda/sort_ofp.cuh"
#include "util/cuda/segreduce_ofp.cuh"
#include "util/cuda/merge_ofp.cuh"

enum flush_type
{
	FLUSH_ON_HOST = 0,
	FLUSH_ON_DEVICE = 1,
	FLUSH_NO_DATA = 2
};

template<typename OfpmVectorT>
using ValueTypeOf = typename std::remove_reference<OfpmVectorT>::type::value_type;

namespace openfpm
{
	// All props
	template<typename sg_type>
	struct htoD
	{
		//! encapsulated source object
		sg_type & sg;

		unsigned int lele;

		htoD(sg_type & sg, unsigned int lele)
		:sg(sg),lele(lele)
		{};


		//! It call the copy function for each property
		template<typename T>
		__device__ __host__ inline void operator()(T& t) const
		{
			sg.template hostToDevice<T::value>(lele,lele);
		}
	};

    constexpr int VECTOR_SPARSE_STANDARD = 1;
    constexpr int VECTOR_SPARSE_BLOCK = 2;

    template<typename reduction_type, unsigned int impl>
    struct cpu_block_process
    {
    	template<typename encap_src, typename encap_dst>
    	static inline void process(encap_src & src, encap_dst & dst)
    	{
    		dst = reduction_type::red(dst,src);
    	}
    };

    template<typename reduction_type>
    struct cpu_block_process<reduction_type,VECTOR_SPARSE_BLOCK>
    {
    	template<typename encap_src, typename encap_dst>
    	static inline void process(encap_src & src, encap_dst & dst)
    	{
    		for (size_t i = 0 ; i < encap_src::size ; i++)
    		{
    			dst[i] = reduction_type::red(dst[i],src[i]);
    		}
    	}
    };

    template<typename reduction_type>
    struct cpu_block_process<reduction_type,3>
    {
    	template<typename encap_src, typename encap_dst,unsigned int N1>
    	static inline void process(encap_src & src, encap_dst (& dst)[N1])
    	{
    		for (unsigned int j = 0 ; j < N1 ; j++)
    		{
				for (size_t i = 0 ; i < encap_dst::size ; i++)
				{
					dst[i] = reduction_type::red(dst[i][j],src[j][i]);
				}
    		}
    	}

    	template<unsigned int N1, unsigned int blockSize, typename encap_src, typename encap_dst>
    	static inline void process_e(encap_src & src, encap_dst & dst)
    	{
    		for (unsigned int j = 0 ; j < N1 ; j++)
    		{
				for (size_t i = 0 ; i < blockSize ; i++)
				{
					dst[i] = reduction_type::red(dst[i][j],src[i][j]);
				}
    		}
    	}
    };

    /*! \brief Functor switch to select the vector sparse for standars scalar and blocked implementation
     *
     *
     */
    template<unsigned int impl, typename block_functor>
    struct scalar_block_implementation_switch // Case for scalar implementations
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

        template <unsigned int pSegment, typename vector_reduction, typename T, typename vector_data_type, typename vector_index_type , typename vector_index_type2>
        static void segreduce(vector_data_type & vector_data,
        		vector_data_type & vector_data_unsorted,
        		vector_index_type & vector_data_map,
                vector_index_type2 & segment_offset,
                vector_data_type & vector_data_red,
                block_functor & blf,
                gpu::ofp_context_t & context)
        {
#ifdef __NVCC__
            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<typename vector_data_type::value_type::type,typename reduction_type::prop>::type red_type;
            typedef typename reduction_type::template op_red<red_type> red_op;
            typedef typename boost::mpl::at<typename vector_index_type::value_type::type,boost::mpl::int_<0>>::type seg_type;
            typename reduction_type::template op_initial_value<red_type> initial_value_functor;

            assert((std::is_same<seg_type,int>::value == true));

            openfpm::segreduce(
                    (red_type *)vector_data.template getDeviceBuffer<reduction_type::prop::value>(), vector_data.size(),
                    (int *)segment_offset.template getDeviceBuffer<1>(), segment_offset.size()-1,
                    (red_type *)vector_data_red.template getDeviceBuffer<reduction_type::prop::value>(),
                    red_op(), initial_value_functor(), context);
#else // __NVCC__
    std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif // __NVCC__
        }


        /*! \briefMerge all datas
         *
         * \param vct_index_old sorted vector of the old indexes
         * \param vct_data_old vector of old data
         * \param vct_index_out output indexes merged new and old indexes
         * \param vct_index_merge_id indicate from where it come from the merged index (if the number is bigger than vct_index_old.size()
         *                                                          the merged index come from the new data)
         * \param vct_index_merge indexes old and new merged with conflicts
         * \param vct_add_data_unique data to add (has been already reduced)
         * \param vct_data_old old data
         * \param vct_add_data data to add in its original form in the insert buffer
         * \param vct_data_out reduced data vector new + old
         * \param vct_index_dtmp temporal buffer vector used for intermediate calculation
         *
         */
        template <
                typename vector_data_type,
                typename vector_index_type,
                typename vector_index_type2,
                typename vector_index_dtmp_type,
                typename Ti,
                typename ... v_reduce>
        static void solveConflicts(
                vector_index_type & vct_index_old,
                vector_index_type & vct_index_merge,
                vector_index_type & vct_index_merge_id,
                vector_index_type & vct_index_out,
                vector_index_dtmp_type & vct_index_dtmp,
                vector_index_type & data_map,
                vector_index_type2 & segments_new,
                vector_data_type & vct_data_old,
                vector_data_type & vct_add_data,
                vector_data_type & vct_add_data_unique,
                vector_data_type & vct_data_out,
                ite_gpu<1> & itew,
                block_functor & blf,
                gpu::ofp_context_t & context
                )
        {
#ifdef __NVCC__

            CUDA_LAUNCH((solve_conflicts<
                        decltype(vct_index_merge.toKernel()),
                        decltype(vct_data_old.toKernel()),
                        decltype(vct_index_dtmp.toKernel()),
                        128,
                        v_reduce ...
                        >),
            			itew,
                                              vct_index_merge.toKernel(),vct_data_old.toKernel(),
                                              vct_index_merge_id.toKernel(),vct_add_data_unique.toKernel(),
                                              vct_index_out.toKernel(),vct_data_out.toKernel(),
                                              vct_index_dtmp.toKernel(),
                                              vct_index_old.size());

                // we scan tmp3
                openfpm::scan(
                        (Ti*)vct_index_dtmp.template getDeviceBuffer<0>(),
                        vct_index_dtmp.size(),
                        (Ti *)vct_index_dtmp.template getDeviceBuffer<1>(),
                        context);

                // get the size to resize vct_index and vct_data
                vct_index_dtmp.template deviceToHost<0,1>(vct_index_dtmp.size()-1,vct_index_dtmp.size()-1);
                int size = vct_index_dtmp.template get<1>(vct_index_dtmp.size()-1) + vct_index_dtmp.template get<0>(vct_index_dtmp.size()-1);

                vct_index_old.resize(size);
                vct_data_old.resize(size);

                CUDA_LAUNCH(realign,itew,vct_index_out.toKernel(),vct_data_out.toKernel(),
                                      vct_index_old.toKernel(), vct_data_old.toKernel(),
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

        template <unsigned int pSegment, typename vector_reduction, typename T, typename vector_data_type, typename vector_index_type ,typename vector_index_type2>
        static void segreduce(vector_data_type & vector_data,
        					  vector_data_type & vector_data_unsorted,
        					  vector_index_type & vector_data_map,
        					  vector_index_type2 & segment_offset,
        					  vector_data_type & vector_data_red,
        					  block_functor & blf,
						  gpu::ofp_context_t & context)
        {

        }

        template <
                typename vector_data_type,
                typename vector_index_type,
                typename vector_index_type2,
                typename vector_index_dtmp_type,
                typename Ti,
                typename ... v_reduce>
        static void solveConflicts(
                vector_index_type & vct_index_old,
                vector_index_type & vct_index_merge,
                vector_index_type & vct_index_merge_id,
                vector_index_type & vct_index_out,
                vector_index_dtmp_type & vct_index_dtmp,
                vector_index_type & data_map,
                vector_index_type2 & segments_new,
                vector_data_type & vct_data,
                vector_data_type & vct_add_data,
                vector_data_type & vct_add_data_unique,
                vector_data_type & vct_data_out,
                ite_gpu<1> & itew,
                block_functor & blf,
                gpu::ofp_context_t & context
        )
        {
#ifdef __NVCC__
            blf.template solve_conflicts<1,
			            decltype(vct_index_merge),
			            decltype(segments_new),
			            decltype(vct_data),
			            v_reduce ...>
			            (vct_index_merge, vct_index_merge_id, segments_new, data_map,
			            vct_data, vct_add_data,
			            vct_index_old, vct_data_out,
			            context);
                vct_data_out.swap(vct_data);

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

	template<typename reduction_type, typename vector_reduction, typename T,unsigned int impl, typename red_type>
	struct sparse_vector_reduction_cpu_impl
	{
		template<typename vector_data_type, typename vector_index_type,typename vector_index_type_reo>
		static inline void red(size_t & i, vector_data_type & vector_data_red,
				   vector_data_type & vector_data,
				   vector_index_type & vector_index,
				   vector_index_type_reo & reorder_add_index_cpu)
		{
			size_t start = reorder_add_index_cpu.get(i).id;
			red_type red = vector_data.template get<reduction_type::prop::value>(i);

			size_t j = 1;
			for ( ; i+j < reorder_add_index_cpu.size() && reorder_add_index_cpu.get(i+j).id == start ; j++)
			{
				cpu_block_process<reduction_type,impl>::process(vector_data.template get<reduction_type::prop::value>(i+j),red);
				//reduction_type::red(red,vector_data.template get<reduction_type::prop::value>(i+j));
			}
			vector_data_red.add();
			vector_data_red.template get<reduction_type::prop::value>(vector_data_red.size()-1) = red;

			if (T::value == 0)
			{
				vector_index.add();
				vector_index.template get<0>(vector_index.size() - 1) = reorder_add_index_cpu.get(i).id;
			}

			i += j;
		}
	};


	template<typename reduction_type, typename vector_reduction, typename T,unsigned int impl, typename red_type, unsigned int N1>
	struct sparse_vector_reduction_cpu_impl<reduction_type,vector_reduction,T,impl,red_type[N1]>
	{
		template<typename vector_data_type, typename vector_index_type,typename vector_index_type_reo>
		static inline void red(size_t & i, vector_data_type & vector_data_red,
				   vector_data_type & vector_data,
				   vector_index_type & vector_index,
				   vector_index_type_reo & reorder_add_index_cpu)
		{
			size_t start = reorder_add_index_cpu.get(i).id;
			red_type red[N1];

			for (size_t k = 0 ; k < N1 ; k++)
			{
				red[k] = vector_data.template get<reduction_type::prop::value>(i)[k];
			}

			size_t j = 1;
			for ( ; i+j < reorder_add_index_cpu.size() && reorder_add_index_cpu.get(i+j).id == start ; j++)
			{
				auto ev = vector_data.template get<reduction_type::prop::value>(i+j);
				cpu_block_process<reduction_type,impl+1>::process(ev,red);
				//reduction_type::red(red,vector_data.template get<reduction_type::prop::value>(i+j));
			}

			vector_data_red.add();

			for (size_t k = 0 ; k < N1 ; k++)
			{
				vector_data_red.template get<reduction_type::prop::value>(vector_data_red.size()-1)[k] = red[k];
			}

			if (T::value == 0)
			{
				vector_index.add();
				vector_index.template get<0>(vector_index.size() - 1) = reorder_add_index_cpu.get(i).id;
			}

			i += j;
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
	        typename vector_index_type_reo,
	        typename vector_reduction,
	        unsigned int impl>
	struct sparse_vector_reduction_cpu
	{
		//! Vector in which to the reduction
		vector_data_type & vector_data_red;

		//! Vector in which to the reduction
		vector_data_type & vector_data;

		//! reorder vector index
		vector_index_type_reo & reorder_add_index_cpu;

		//! Index type vector
		vector_index_type & vector_index;

		/*! \brief constructor
		 *
		 * \param src source encapsulated object
		 * \param dst source encapsulated object
		 *
		 */
		inline sparse_vector_reduction_cpu(vector_data_type & vector_data_red,
									   vector_data_type & vector_data,
									   vector_index_type & vector_index,
									   vector_index_type_reo & reorder_add_index_cpu)
		:vector_data_red(vector_data_red),vector_data(vector_data),vector_index(vector_index),reorder_add_index_cpu(reorder_add_index_cpu)
		{};

		//! It call the copy function for each property
		template<typename T>
		inline void operator()(T& t) const
		{
            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<typename ValueTypeOf<vector_data_type>::type,typename reduction_type::prop>::type red_type;

            if (reduction_type::is_special() == false)
			{
    			for (size_t i = 0 ; i < reorder_add_index_cpu.size() ; )
    			{
    				sparse_vector_reduction_cpu_impl<reduction_type,vector_reduction,T,impl,red_type>::red(i,vector_data_red,vector_data,vector_index,reorder_add_index_cpu);

/*    				size_t start = reorder_add_index_cpu.get(i).id;
    				red_type red = vector_data.template get<reduction_type::prop::value>(i);

    				size_t j = 1;
    				for ( ; i+j < reorder_add_index_cpu.size() && reorder_add_index_cpu.get(i+j).id == start ; j++)
    				{
    					cpu_block_process<reduction_type,impl>::process(vector_data.template get<reduction_type::prop::value>(i+j),red);
    					//reduction_type::red(red,vector_data.template get<reduction_type::prop::value>(i+j));
    				}
    				vector_data_red.add();
    				vector_data_red.template get<reduction_type::prop::value>(vector_data_red.size()-1) = red;

    				if (T::value == 0)
    				{
    					vector_index.add();
    					vector_index.template get<0>(vector_index.size() - 1) = reorder_add_index_cpu.get(i).id;
    				}

    				i += j;*/
    			}
			}
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
	template<typename encap_src,
			typename encap_dst,
	        typename vector_reduction>
	struct sparse_vector_reduction_solve_conflict_assign_cpu
	{
		//! source
		encap_src & src;

		//! destination
		encap_dst & dst;


		/*! \brief constructor
		 *
		 * \param src source encapsulated object
		 * \param dst source encapsulated object
		 *
		 */
		inline sparse_vector_reduction_solve_conflict_assign_cpu(encap_src & src, encap_dst & dst)
		:src(src),dst(dst)
		{};

		//! It call the copy function for each property
		template<typename T>
		inline void operator()(T& t) const
		{
            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;

            dst.template get<reduction_type::prop::value>() = src.template get<reduction_type::prop::value>();
		}
	};


	template<unsigned int impl,typename vector_reduction, typename T,typename red_type>
	struct sparse_vector_reduction_solve_conflict_reduce_cpu_impl
	{
		template<typename encap_src, typename encap_dst>
		static inline void red(encap_src & src, encap_dst & dst)
		{
			typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;

			cpu_block_process<reduction_type,impl>::process(src.template get<reduction_type::prop::value>(),dst.template get<reduction_type::prop::value>());
		}
	};

	template<unsigned int impl, typename vector_reduction, typename T,typename red_type, unsigned int N1>
	struct sparse_vector_reduction_solve_conflict_reduce_cpu_impl<impl,vector_reduction,T,red_type[N1]>
	{
		template<typename encap_src, typename encap_dst>
		static inline void red(encap_src & src, encap_dst & dst)
		{
            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;

			auto src_e = src.template get<reduction_type::prop::value>();
			auto dst_e = dst.template get<reduction_type::prop::value>();

			cpu_block_process<reduction_type,impl+1>::template process_e<N1,red_type::size>(src_e,dst_e);
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
	template<typename encap_src,
			typename encap_dst,
	        typename vector_reduction,
	        unsigned int impl>
	struct sparse_vector_reduction_solve_conflict_reduce_cpu
	{
		//! source
		encap_src & src;

		//! destination
		encap_dst & dst;


		/*! \brief constructor
		 *
		 * \param src source encapsulated object
		 * \param dst source encapsulated object
		 *
		 */
		inline sparse_vector_reduction_solve_conflict_reduce_cpu(encap_src & src, encap_dst & dst)
		:src(src),dst(dst)
		{};

		//! It call the copy function for each property
		template<typename T>
		inline void operator()(T& t) const
		{
            typedef typename boost::mpl::at<vector_reduction, T>::type reduction_type;
            typedef typename boost::mpl::at<typename encap_src::T_type::type, typename reduction_type::prop>::type red_type;

            sparse_vector_reduction_solve_conflict_reduce_cpu_impl<impl,vector_reduction,T,red_type>::red(src,dst);

//            cpu_block_process<reduction_type,impl>::process(src.template get<reduction_type::prop::value>(),dst.template get<reduction_type::prop::value>());
//            reduction_type::red(dst.template get<reduction_type::prop::value>(),src.template get<reduction_type::prop::value>());
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
	        typename vector_index_type2,
	        typename vector_reduction,
	        typename block_functor,
	        unsigned int impl2, unsigned int pSegment=1>
	struct sparse_vector_reduction
	{
		//! Vector in which to the reduction
		vector_data_type & vector_data_red;

		//! new datas
		vector_data_type & vector_data;

		//! new data in an unsorted way
		vector_data_type & vector_data_unsorted;

		//! segment of offsets
		vector_index_type2 & segment_offset;

		//! map of the data
		vector_index_type & vector_data_map;

		//! block functor
		block_functor & blf;

		//! gpu context
		gpu::ofp_context_t & context;

		/*! \brief constructor
		 *
		 * \param src source encapsulated object
		 * \param dst source encapsulated object
		 *
		 */
		inline sparse_vector_reduction(vector_data_type & vector_data_red,
									   vector_data_type & vector_data,
									   vector_data_type & vector_data_unsorted,
									   vector_index_type & vector_data_map,
									   vector_index_type2 & segment_offset,
									   block_functor & blf,
									   gpu::ofp_context_t & context)
		:vector_data_red(vector_data_red),
		 vector_data(vector_data),
		 vector_data_unsorted(vector_data_unsorted),
		 segment_offset(segment_offset),
		 vector_data_map(vector_data_map),
		 blf(blf),
		 context(context)
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
			            vector_data_unsorted,
			            vector_data_map,
			            segment_offset,
			            vector_data_red,
			            blf,
			            context);
			}
#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: this file is supposed to be compiled with nvcc" << std::endl;
#endif
		}
	};


	struct stub_block_functor
	{
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
                                    gpu::ofp_context_t & context)
		{
			return true;
		}

		openfpm::vector_gpu<aggregate<unsigned int>> outputMap;

        openfpm::vector_gpu<aggregate<unsigned int>> & get_outputMap()
		{
        	return outputMap;
		}

        const openfpm::vector_gpu<aggregate<unsigned int>> & get_outputMap() const
		{
        	return outputMap;
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
		gpu::ofp_context_t & context;

		/*! \brief constructor
		 *
		 * \param src source encapsulated object
		 * \param dst source encapsulated object
		 *
		 */
		inline sparse_vector_special(vector_data_type & vector_data_red,
									   vector_data_type & vector_data,
									   vector_index_type & segment_offset,
									   gpu::ofp_context_t & context)
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

				CUDA_LAUNCH((reduce_from_offset<decltype(segment_offset.toKernel()),decltype(vector_data_red.toKernel()),reduction_type>),
															ite,segment_offset.toKernel(),vector_data_red.toKernel(),vector_data.size());
			}

#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: this file si supposed to be compiled with nvcc" << std::endl;
#endif
		}
	};

	template<typename T,
			 typename Ti = long int,
			 typename Memory=HeapMemory,
			 typename layout=typename memory_traits_lin<T>::type,
			 template<typename> class layout_base=memory_traits_lin ,
			 typename grow_p=grow_policy_double,
			 unsigned int impl=vect_isel<T>::value,
			 unsigned int impl2 = VECTOR_SPARSE_STANDARD,
			 typename block_functor = stub_block_functor>
	class vector_sparse
	{
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_index;
		vector<T,Memory,layout_base,grow_p,impl> vct_data;
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_m_index;

		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_add_index;
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_rem_index;
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_nadd_index;
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_nrem_index;
		vector<T,Memory,layout_base,grow_p> vct_add_data;
		vector<T,Memory,layout_base,grow_p> vct_add_data_reord;

		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_add_index_cont_0;
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_add_index_cont_1;
		vector<T,Memory,layout_base,grow_p> vct_add_data_cont;
		vector<aggregate<Ti,Ti>,Memory,layout_base,grow_p> vct_add_index_unique;
		vector<aggregate<int,int>,Memory,layout_base,grow_p> segments_int;

		vector<T,Memory,layout_base,grow_p,impl> vct_add_data_unique;

		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_index_tmp4;
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_index_tmp;
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_index_tmp2;
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_index_tmp3;
		vector<aggregate<Ti,Ti,Ti>,Memory,layout_base,grow_p> vct_index_dtmp;

		// segments map (This is used only in case of Blocked data)
		vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_segment_index_map;

		block_functor blf;

		T bck;

		CudaMemory mem;

		openfpm::vector<reorder<Ti>> reorder_add_index_cpu;

		size_t max_ele;

		int n_gpu_add_block_slot = 0;
		int n_gpu_rem_block_slot = 0;

		/*! \brief get the element i
		 *
		 * search the element x
		 *
		 * \param i element i
		 */
		template<bool prefetch>
		inline Ti _branchfree_search_nobck(Ti x, Ti & id) const
		{
			if (vct_index.size() == 0)	{id = 0; return -1;}
			const Ti *base = &vct_index.template get<0>(0);
			const Ti *end = (const Ti *)vct_index.template getPointer<0>() + vct_index.size();
			Ti n = vct_data.size()-1;
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
			return (base + off != end)?*(base + off):-1;
		}

		/*! \brief get the element i
		 *
		 * search the element x
		 *
		 * \param i element i
		 */
		template<bool prefetch>
		inline void _branchfree_search(Ti x, Ti & id) const
		{
			Ti v = _branchfree_search_nobck<prefetch>(x,id);
			id = (x == v)?id:vct_data.size()-1;
		}


		/* \brief take the indexes for the insertion pools and create a continuos array
		 *
		 * \param vct_nadd_index number of insertions of each pool
		 * \param vct_add_index pool of inserts
		 * \param vct_add_cont_index output continuos array of inserted indexes
		 * \param vct_add_data array of added data
		 * \param vct_add_data_cont continuos array of inserted data
		 * \param contect gpu context
		 *
		 */
		size_t make_continuos(vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_nadd_index,
							  vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_index,
							  vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_cont_index,
							  vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_cont_index_map,
							  vector<T,Memory,layout_base,grow_p> & vct_add_data,
							  vector<T,Memory,layout_base,grow_p> & vct_add_data_cont,
							  gpu::ofp_context_t & context)
		{
#ifdef __NVCC__

			// Add 0 to the last element to vct_nadd_index
			vct_nadd_index.resize(vct_nadd_index.size()+1);
			vct_nadd_index.template get<0>(vct_nadd_index.size()-1) = 0;
			vct_nadd_index.template hostToDevice<0>(vct_nadd_index.size()-1,vct_nadd_index.size()-1);

			// Merge the list of inserted points for each block
			vct_index_tmp4.resize(vct_nadd_index.size());

			openfpm::scan((Ti *)vct_nadd_index.template getDeviceBuffer<0>(),
			            vct_nadd_index.size(),
			            (Ti *)vct_index_tmp4.template getDeviceBuffer<0>() ,
                        context);

			vct_index_tmp4.template deviceToHost<0>(vct_index_tmp4.size()-1,vct_index_tmp4.size()-1);
			size_t n_ele = vct_index_tmp4.template get<0>(vct_index_tmp4.size()-1);

			// we reuse vct_nadd_index
			vct_add_cont_index.resize(n_ele);
			vct_add_cont_index_map.resize(n_ele);

			if (impl2 == VECTOR_SPARSE_STANDARD)
			{
				vct_add_data_cont.resize(n_ele);
			}
			else
			{
				vct_segment_index_map.resize(n_ele);
			}

			if (n_gpu_add_block_slot >= 128)
			{
				ite_gpu<1> itew;
				itew.wthr.x = vct_nadd_index.size()-1;
				itew.wthr.y = 1;
				itew.wthr.z = 1;
				itew.thr.x = 128;
				itew.thr.y = 1;
				itew.thr.z = 1;

				CUDA_LAUNCH(construct_insert_list_key_only,itew,vct_add_index.toKernel(),
									vct_nadd_index.toKernel(),
									vct_index_tmp4.toKernel(),
									vct_add_cont_index.toKernel(),
									vct_add_cont_index_map.toKernel(),
									n_gpu_add_block_slot);
			}
			else
			{
				auto itew = vct_add_index.getGPUIterator();

				CUDA_LAUNCH(construct_insert_list_key_only_small_pool,itew,vct_add_index.toKernel(),
									vct_nadd_index.toKernel(),
									vct_index_tmp4.toKernel(),
									vct_add_cont_index.toKernel(),
									vct_add_cont_index_map.toKernel(),
									n_gpu_add_block_slot);
			}

			return n_ele;
#endif
			return 0;
		}

		/*! \brief sort the continuos array of inserted key
		 *
		 * \param context modern gpu context
		 * \param vct_add_cont_index array of indexes (unsorted), as output will be sorted
		 * \param vct_add_cont_index_map reference to the original indexes
		 * \param vct_add_data_reord sorted data output
		 * \param vct_add_data_cont added data in a continuos unsorted array
		 *
		 */
		void reorder_indexes(vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_cont_index,
							 vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_cont_index_map,
							 vector<T,Memory,layout_base,grow_p> & vct_add_data_reord,
							 vector<T,Memory,layout_base,grow_p> & vct_add_data_cont,
							 gpu::ofp_context_t & context)
		{
#ifdef __NVCC__
			ite_gpu<1> itew;
			itew.wthr.x = vct_nadd_index.size()-1;
			itew.wthr.y = 1;
			itew.wthr.z = 1;
			itew.thr.x = 128;
			itew.thr.y = 1;
			itew.thr.z = 1;

			size_t n_ele = vct_add_cont_index.size();

			n_gpu_add_block_slot = 0;

			// now we sort
			openfpm::sort(
			        (Ti *)vct_add_cont_index.template getDeviceBuffer<0>(),
                    (Ti *)vct_add_cont_index_map.template getDeviceBuffer<0>(),
					vct_add_cont_index.size(),
					gpu::template less_t<Ti>(),
                    context);

			auto ite = vct_add_cont_index.getGPUIterator();

			// Now we reorder the data vector accordingly to the indexes

			if (impl2 == VECTOR_SPARSE_STANDARD)
			{
				vct_add_data_reord.resize(n_ele);
				CUDA_LAUNCH(reorder_vector_data,ite,vct_add_cont_index_map.toKernel(),vct_add_data_cont.toKernel(),vct_add_data_reord.toKernel());
			}

#endif
		}

		/*! \brief Merge indexes
		 *
		 * \param
		 *
		 *
		 */
		template<typename ... v_reduce>
		void merge_indexes(vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_index_sort,
						   vector<aggregate<Ti,Ti>,Memory,layout_base,grow_p> & vct_add_index_unique,
				  	  	   vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_merge_index,
				  	  	   vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_merge_index_map,
						   gpu::ofp_context_t & context)
		{
#ifdef __NVCC__

			typedef boost::mpl::vector<v_reduce...> vv_reduce;

			auto ite = vct_add_index_sort.getGPUIterator();

			mem.allocate(sizeof(int));
			mem.fill(0);
			vct_add_index_unique.resize(vct_add_index_sort.size()+1);

			ite = vct_add_index_sort.getGPUIterator();

			vct_index_tmp4.resize(vct_add_index_sort.size()+1);

			CUDA_LAUNCH(
					(
						find_buffer_offsets_for_scan
								<0,
								decltype(vct_add_index_sort.toKernel()),
								decltype(vct_index_tmp4.toKernel())
								>
					),
					ite,
					vct_add_index_sort.toKernel(),
					vct_index_tmp4.toKernel());

			openfpm::scan((Ti *)vct_index_tmp4.template getDeviceBuffer<0>(),vct_index_tmp4.size(),(Ti *)vct_index_tmp4.template getDeviceBuffer<0>(),context);

			vct_index_tmp4.template deviceToHost<0>(vct_index_tmp4.size()-1,vct_index_tmp4.size()-1);
			int n_ele_unique = vct_index_tmp4.template get<0>(vct_index_tmp4.size()-1);

			vct_add_index_unique.resize(n_ele_unique);

			if (impl2 == VECTOR_SPARSE_STANDARD)
			{
				vct_add_data_unique.resize(n_ele_unique);
			}

			CUDA_LAUNCH(
					(construct_index_unique<0>),
					ite,
					vct_add_index_sort.toKernel(),
					vct_index_tmp4.toKernel(),
					vct_add_index_unique.toKernel());

			typedef boost::mpl::vector<v_reduce...> vv_reduce;

			// Then we merge the two list vct_index and vct_add_index_unique

			// index to get merge index
			vct_m_index.resize(vct_index.size());

			if (vct_m_index.size() != 0)
			{
				ite = vct_m_index.getGPUIterator();
				CUDA_LAUNCH((set_indexes<0>),ite,vct_m_index.toKernel(),0);
			}

			// after merge we solve the last conflicts, running across the vector again and spitting 1 when there is something to merge
			// we reorder the data array also

			vct_merge_index.resize(vct_index.size() + vct_add_index_unique.size());
			vct_merge_index_map.resize(vct_index.size() + vct_add_index_unique.size());
			vct_index_tmp3.resize(vct_index.size() + vct_add_index_unique.size());

			// Do not delete this reserve
			// Unfortunately all resize with DataBlocks are broken
			if (impl2 == VECTOR_SPARSE_STANDARD)
			{
				vct_add_data_cont.reserve(vct_index.size() + vct_add_index_unique.size()+1);
				vct_add_data_cont.resize(vct_index.size() + vct_add_index_unique.size());
			}

			ite = vct_add_index_unique.getGPUIterator();
			vct_index_tmp4.resize(vct_add_index_unique.size());
			CUDA_LAUNCH((set_indexes<0>),ite,vct_index_tmp4.toKernel(),vct_index.size());

			ite_gpu<1> itew;

			itew.wthr.x = vct_merge_index.size() / 128 + (vct_merge_index.size() % 128 != 0);
			itew.wthr.y = 1;
			itew.wthr.z = 1;
			itew.thr.x = 128;
			itew.thr.y = 1;
			itew.thr.z = 1;

			vct_index_dtmp.resize(itew.wthr.x);

			// we merge with vct_index with vct_add_index_unique in vct_merge_index, vct_merge_index contain the merged indexes
			//

			openfpm::merge((Ti *)vct_index.template getDeviceBuffer<0>(),(Ti *)vct_m_index.template getDeviceBuffer<0>(),vct_index.size(),
						(Ti *)vct_add_index_unique.template getDeviceBuffer<0>(),(Ti *)vct_index_tmp4.template getDeviceBuffer<0>(),vct_add_index_unique.size(),
						(Ti *)vct_merge_index.template getDeviceBuffer<0>(),(Ti *)vct_merge_index_map.template getDeviceBuffer<0>(),gpu::less_t<Ti>(),context);


#endif
		}



		template<typename ... v_reduce>
		void merge_datas(vector<T,Memory,layout_base,grow_p> & vct_add_data_reord,
						 vector<aggregate<Ti,Ti>,Memory,layout_base,grow_p> & segments_new,
						 vector<T,Memory,layout_base,grow_p> & vct_add_data,
						 vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_data_reord_map,
						 gpu::ofp_context_t & context)
		{
#ifdef __NVCC__
			ite_gpu<1> itew;
			itew.wthr.x = vct_index_tmp.size() / 128 + (vct_index_tmp.size() % 128 != 0);
			itew.wthr.y = 1;
			itew.wthr.z = 1;
			itew.thr.x = 128;
			itew.thr.y = 1;
			itew.thr.z = 1;

			typedef boost::mpl::vector<v_reduce...> vv_reduce;

			////////////////////////////////////////////////////////////////////////////////////////////////////

			// Now we can do a segmented reduction
			scalar_block_implementation_switch<impl2, block_functor>
			        ::template extendSegments<1>(vct_add_index_unique, vct_add_data_reord_map.size());

			if (impl2 == VECTOR_SPARSE_STANDARD)
			{
				sparse_vector_reduction<typename std::remove_reference<decltype(vct_add_data)>::type,
								    decltype(vct_add_data_reord_map),
								    decltype(vct_add_index_unique),vv_reduce,block_functor,impl2>
			        svr(
			                vct_add_data_unique,
			                vct_add_data_reord,
			                vct_add_data,
			                vct_add_data_reord_map,
			                vct_add_index_unique,
			                blf,
			                context);

				boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(svr);
				vct_add_index_unique.remove(vct_add_index_unique.size()-1);
			}

			sparse_vector_special<typename std::remove_reference<decltype(vct_add_data)>::type,
								  decltype(vct_add_index_unique),
								  vv_reduce> svr2(vct_add_data_unique,vct_add_data_reord,vct_add_index_unique,context);
			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(svr2);

			//////////////////////////////////////////////////////////////////////////////////////////////////////

            // Now perform the right solve_conflicts according to impl2
            scalar_block_implementation_switch<impl2, block_functor>::template solveConflicts<
                    decltype(vct_data),
                    decltype(vct_index),
                    decltype(segments_new),
                    decltype(vct_index_dtmp),
                    Ti,
                    v_reduce ...
                    >
                    (
                        vct_index,
                        vct_index_tmp,
                        vct_index_tmp2,
                        vct_index_tmp3,
                        vct_index_dtmp,
                        vct_add_data_reord_map,
                        segments_new,
                        vct_data,
                        vct_add_data,
                        vct_add_data_unique,
                        vct_add_data_cont,
                        itew,
                        blf,
                        context
                    );


#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif
		}

		template<typename ... v_reduce>
		void flush_on_gpu_insert(vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_index_cont_0,
				  vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_index_cont_1,
				  vector<T,Memory,layout_base,grow_p> & vct_add_data_reord,
				  gpu::ofp_context_t & context)
		{
#ifdef __NVCC__

			// To avoid the case where you never called setGPUInsertBuffer
			if (n_gpu_add_block_slot == 0 || vct_add_index.size() == 0)
			{
				return;
			}

			size_t n_ele = make_continuos(vct_nadd_index,vct_add_index,vct_add_index_cont_0,vct_add_index_cont_1,
										  vct_add_data,vct_add_data_cont,context);

            // At this point we can check whether we have not inserted anything actually,
            // in this case, return without further ado...
			if (vct_add_index_cont_0.size() == 0)
            {return;}

			reorder_indexes(vct_add_index_cont_0,vct_add_index_cont_1,vct_add_data_reord,vct_add_data,context);

			merge_indexes<v_reduce ... >(vct_add_index_cont_0,vct_add_index_unique,
										 vct_index_tmp,vct_index_tmp2,
										 context);

			merge_datas<v_reduce ... >(vct_add_data_reord,vct_add_index_unique,vct_add_data,vct_add_index_cont_1,context);

#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif
		}


		void flush_on_gpu_remove(
				  gpu::ofp_context_t & context)
		{
#ifdef __NVCC__

			// Add 0 to the last element to vct_nadd_index
			vct_nrem_index.resize(vct_nrem_index.size()+1);
			vct_nrem_index.template get<0>(vct_nrem_index.size()-1) = 0;
			vct_nrem_index.template hostToDevice<0>(vct_nrem_index.size()-1,vct_nrem_index.size()-1);

			// Merge the list of inserted points for each block
			vct_index_tmp4.resize(vct_nrem_index.size());

			openfpm::scan((Ti *)vct_nrem_index.template getDeviceBuffer<0>(), vct_nrem_index.size(), (Ti *)vct_index_tmp4.template getDeviceBuffer<0>() , context);

			vct_index_tmp4.template deviceToHost<0>(vct_index_tmp4.size()-1,vct_index_tmp4.size()-1);
			size_t n_ele = vct_index_tmp4.template get<0>(vct_index_tmp4.size()-1);

			// we reuse vct_nadd_index
			vct_add_index_cont_0.resize(n_ele);
			vct_add_index_cont_1.resize(n_ele);

			ite_gpu<1> itew;
			itew.wthr.x = vct_nrem_index.size()-1;
			itew.wthr.y = 1;
			itew.wthr.z = 1;
			itew.thr.x = 128;
			itew.thr.y = 1;
			itew.thr.z = 1;

			CUDA_LAUNCH(construct_remove_list,itew,vct_rem_index.toKernel(),
										vct_nrem_index.toKernel(),
										vct_index_tmp4.toKernel(),
										vct_add_index_cont_0.toKernel(),
										vct_add_index_cont_1.toKernel(),
										n_gpu_rem_block_slot);

			// now we sort
			openfpm::sort((Ti *)vct_add_index_cont_0.template getDeviceBuffer<0>(),(Ti *)vct_add_index_cont_1.template getDeviceBuffer<0>(),
					vct_add_index_cont_0.size(), gpu::template less_t<Ti>(), context);

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

			openfpm::sort((Ti *)vct_add_index_unique.template getDeviceBuffer<1>(),(Ti *)vct_add_index_unique.template getDeviceBuffer<0>(),
							vct_add_index_unique.size(),gpu::template less_t<Ti>(),context);

			// Then we merge the two list vct_index and vct_add_index_unique

			// index to get merge index
			vct_m_index.resize(vct_index.size() + vct_add_index_unique.size());

			ite = vct_m_index.getGPUIterator();
			CUDA_LAUNCH((set_indexes<0>),ite,vct_m_index.toKernel(),0);

			ite = vct_add_index_unique.getGPUIterator();
			CUDA_LAUNCH((set_indexes<1>),ite,vct_add_index_unique.toKernel(),vct_index.size());

			// after merge we solve the last conflicts, running across the vector again and spitting 1 when there is something to merge
			// we reorder the data array also

			vct_index_tmp.resize(vct_index.size() + vct_add_index_unique.size());
			vct_index_tmp2.resize(vct_index.size() + vct_add_index_unique.size());

			itew.wthr.x = vct_index_tmp.size() / 128 + (vct_index_tmp.size() % 128 != 0);
			itew.wthr.y = 1;
			itew.wthr.z = 1;
			itew.thr.x = 128;
			itew.thr.y = 1;
			itew.thr.z = 1;

			vct_index_dtmp.resize(itew.wthr.x);

			// we merge with vct_index with vct_add_index_unique in vct_index_tmp, vct_intex_tmp2 contain the merging index
			//
			openfpm::merge((Ti *)vct_index.template getDeviceBuffer<0>(),(Ti *)vct_m_index.template getDeviceBuffer<0>(),vct_index.size(),
						(Ti *)vct_add_index_unique.template getDeviceBuffer<0>(),(Ti *)vct_add_index_unique.template getDeviceBuffer<1>(),vct_add_index_unique.size(),
						(Ti *)vct_index_tmp.template getDeviceBuffer<0>(),(Ti *)vct_index_tmp2.template getDeviceBuffer<0>(),gpu::less_t<Ti>(),context);

			vct_index_tmp3.resize(128*itew.wthr.x);

			CUDA_LAUNCH((solve_conflicts_remove<decltype(vct_index_tmp.toKernel()),decltype(vct_index_dtmp.toKernel()),128>),
										itew,
										vct_index_tmp.toKernel(),
										 vct_index_tmp2.toKernel(),
										 vct_index_tmp3.toKernel(),
										 vct_m_index.toKernel(),
										  vct_index_dtmp.toKernel(),
										  vct_index.size());

			// we scan tmp3
			openfpm::scan((Ti*)vct_index_dtmp.template getDeviceBuffer<0>(),vct_index_dtmp.size(),(Ti *)vct_index_dtmp.template getDeviceBuffer<1>(),context);

			// get the size to resize vct_index and vct_data
			vct_index_dtmp.template deviceToHost<0,1>(vct_index_dtmp.size()-1,vct_index_dtmp.size()-1);
			int size = vct_index_dtmp.template get<1>(vct_index_dtmp.size()-1) + vct_index_dtmp.template get<0>(vct_index_dtmp.size()-1);

			vct_add_data_cont.resize(size);
			vct_index.resize(size);

			CUDA_LAUNCH(realign_remove,itew,vct_index_tmp3.toKernel(),vct_m_index.toKernel(),vct_data.toKernel(),
								  vct_index.toKernel(),vct_add_data_cont.toKernel(),
								  vct_index_dtmp.toKernel());

			vct_data.swap(vct_add_data_cont);

#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: you are suppose to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif
		}

		void resetBck()
		{
			// re-add background
			vct_data.resize(vct_data.size()+1);
			vct_data.get(vct_data.size()-1) = bck;

			htoD<decltype(vct_data)> trf(vct_data,vct_data.size()-1);
			boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(trf);
		}

		template<typename ... v_reduce>
		void flush_on_gpu(vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_index_cont_0,
						  vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_index_cont_1,
						  vector<T,Memory,layout_base,grow_p> & vct_add_data_reord,
						  gpu::ofp_context_t & context)
		{
			flush_on_gpu_insert<v_reduce ... >(vct_add_index_cont_0,vct_add_index_cont_1,vct_add_data_reord,context);
		}

		template<typename ... v_reduce>
		void flush_on_cpu()
		{
			if (vct_add_index.size() == 0)
			{return;}

			// First copy the added index to reorder
			reorder_add_index_cpu.resize(vct_add_index.size());
			vct_add_data_cont.resize(vct_add_index.size());

			for (size_t i = 0 ; i < reorder_add_index_cpu.size() ; i++)
			{
				reorder_add_index_cpu.get(i).id = vct_add_index.template get<0>(i);
				reorder_add_index_cpu.get(i).id2 = i;
			}

			reorder_add_index_cpu.sort();

			// Copy the data
			for (size_t i = 0 ; i < reorder_add_index_cpu.size() ; i++)
			{
				vct_add_data_cont.get(i) = vct_add_data.get(reorder_add_index_cpu.get(i).id2);
			}

			typedef boost::mpl::vector<v_reduce...> vv_reduce;

			sparse_vector_reduction_cpu<decltype(vct_add_data),
										decltype(vct_add_index_unique),
										decltype(reorder_add_index_cpu),
										vv_reduce,
										impl2>
			        svr(vct_add_data_unique,
			        	vct_add_data_cont,
			        	vct_add_index_unique,
			        	reorder_add_index_cpu);

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(svr);

			// merge the the data

			vector<T,Memory,layout_base,grow_p,impl> vct_data_tmp;
			vector<aggregate<Ti>,Memory,layout_base,grow_p> vct_index_tmp;

			vct_data_tmp.resize(vct_data.size() + vct_add_data_unique.size());
			vct_index_tmp.resize(vct_index.size() + vct_add_index_unique.size());

			Ti di = 0;
			Ti ai = 0;
			size_t i = 0;

			for ( ; i < vct_data_tmp.size() ; i++)
			{
				Ti id_a = (ai < vct_add_index_unique.size())?vct_add_index_unique.template get<0>(ai):std::numeric_limits<Ti>::max();
				Ti id_d = (di < vct_index.size())?vct_index.template get<0>(di):std::numeric_limits<Ti>::max();

				if (  id_a <= id_d )
				{
					vct_index_tmp.template get<0>(i) = id_a;

					if (id_a == id_d)
					{
						auto dst = vct_data_tmp.get(i);
						auto src = vct_add_data_unique.get(ai);

						sparse_vector_reduction_solve_conflict_assign_cpu<decltype(vct_data_tmp.get(i)),
																		  decltype(vct_add_data.get(ai)),
																		  vv_reduce>
						sva(src,dst);

						boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(sva);
						ai++;

						dst = vct_data_tmp.get(i);
						src = vct_data.get(di);

						sparse_vector_reduction_solve_conflict_reduce_cpu<decltype(vct_data_tmp.get(i)),
								  	  	  	  	  	  	  	  	  	  	  decltype(vct_data.get(di)),
								  	  	  	  	  	  	  	  	  	  	  vv_reduce,
								  	  	  	  	  	  	  	  	  	  	  impl2>
						svr(src,dst);
						boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(svr);

						di++;

						vct_data_tmp.resize(vct_data_tmp.size()-1);
						vct_index_tmp.resize(vct_index_tmp.size()-1);
					}
					else
					{
						vct_index_tmp.template get<0>(i) = vct_add_index_unique.template get<0>(ai);
						vct_data_tmp.get(i) = vct_add_data_unique.get(ai);
						ai++;
					}
				}
				else
				{
					vct_index_tmp.template get<0>(i) = vct_index.template get<0>(di);
					vct_data_tmp.get(i) = vct_data.get(di);
					di++;
				}
			}

			vct_index.swap(vct_index_tmp);
			vct_data.swap(vct_data_tmp);

			vct_add_data.clear();
			vct_add_index.clear();
			vct_add_index_unique.clear();
			vct_add_data_unique.clear();
		}

	public:

		vector_sparse()
		:max_ele(0)
		{
			vct_data.resize(1);
		}

        /*! \brief Get the indices buffer
        *
        * \return the reference to the indices buffer
        */
        auto getIndexBuffer() -> decltype(vct_index)&
        {
            return vct_index;
        }

        /*! \brief Get the data buffer
         *
         * \return the reference to the data buffer
         */
        auto getDataBuffer() -> decltype(vct_data)&
        {
            return vct_data;
        }

        /*! \brief Get the indices buffer
        *
        * \return the reference to the indices buffer
        */
        auto getIndexBuffer() const -> const decltype(vct_index)&
        {
            return vct_index;
        }

        /*! \brief Get the data buffer
         *
         * \return the reference to the data buffer
         */
        auto getDataBuffer() const -> const decltype(vct_data)&
        {
            return vct_data;
        }

		/*! \brief Get the sparse index
		 *
		 * Get the sparse index of the element id
		 *
		 * \note use get_index and get to retrieve the value index associated to the sparse index
		 *
		 * \param id Element to get
		 *
		 * \return the element value requested
		 *
		 */
		inline openfpm::sparse_index<Ti> get_sparse(Ti id) const
		{
			Ti di;
			this->_branchfree_search<false>(id,di);
			openfpm::sparse_index<Ti> sid;
			sid.id = di;

			return sid;
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
			this->_branchfree_search<false>(id,di);
			return vct_data.template get<p>(di);
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
		inline auto get(Ti id) const -> decltype(vct_data.get(id))
		{
			Ti di;
			this->_branchfree_search<false>(id,di);
			return vct_data.get(di);
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
		void swapIndexVector(vector<aggregate<Ti>,Memory,layout_base,grow_p> & iv)
		{
			vct_index.swap(iv);
		}

		/*! \brief Set the background to bck (which value get must return when the value is not find)
		 *
		 * \param bck
		 *
		 */
		template <unsigned int p>
		auto getBackground() const -> decltype(vct_data.template get<p>(vct_data.size()-1))
		{
			return vct_data.template get<p>(vct_data.size()-1);
		}

		/*! \brief Set the background to bck (which value get must return when the value is not find)
		 *
		 * \param bck
		 *
		 */
		auto getBackground() const -> decltype(vct_data.get(vct_data.size()-1))
		{
			return vct_data.get(vct_data.size()-1);
		}

	    template<unsigned int p>
	    void setBackground(const typename boost::mpl::at<typename T::type, boost::mpl::int_<p>>::type & bck_)
	    {
	    	meta_copy_d<typename boost::mpl::at<typename T::type, boost::mpl::int_<p>>::type,
	    				typename std::remove_reference<decltype(vct_data.template get<p>(vct_data.size()-1))>::type>
	    				::meta_copy_d_(bck_,vct_data.template get<p>(vct_data.size()-1));

	    	vct_data.template hostToDevice<p>(vct_data.size()-1,vct_data.size()-1);

	    	meta_copy<typename boost::mpl::at<typename T::type, boost::mpl::int_<p>>::type>
	    				::meta_copy_(bck_,bck.template get<p>());
	    }

		/*! \brief It insert an element in the sparse vector
		 *
		 * \tparam p property id
		 *
		 * \param ele element id
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

		/*! \brief It insert an element in the sparse vector
		 *
		 * \tparam p property id
		 *
		 * \param ele element id
		 *
		 */
		template <unsigned int p>
		auto insertFlush(Ti ele, bool & is_new) -> decltype(vct_data.template get<p>(0))
		{
			is_new = false;
			size_t di;

			// first we have to search if the block exist
			Ti v = _branchfree_search_nobck(ele,di);

			if (v == ele)
			{
				// block exist
				return vct_data.template get<p>(di);
			}
			is_new = true;
			
			// It does not exist, we create it di contain the index where we have to create the new block
			vct_index.insert(di);
			vct_data.isert(di);

			return vct_data.template get<p>(di);
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 * \param ele element id
		 *
		 */
		auto insertFlush(Ti ele, bool & is_new) -> decltype(vct_data.get(0))
		{
			is_new = false;
			Ti di;

			// first we have to search if the block exist
			Ti v = _branchfree_search_nobck<true>(ele,di);

			if (v == ele)
			{
				// block exist
				return vct_data.get(di);
			}

			// It does not exist, we create it di contain the index where we have to create the new block
			vct_index.insert(di);
			vct_data.insert(di);
			is_new = true;

			vct_index.template get<0>(di) = ele;

			return vct_data.get(di);
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 * \param ele element id
		 *
		 */
		auto insert(Ti ele) -> decltype(vct_data.get(0))
		{
			vct_add_index.add();
			vct_add_index.template get<0>(vct_add_index.size()-1) = ele;
			vct_add_data.add();
			return vct_add_data.get(vct_add_data.size()-1);
		}

		/*! \brief merge the added element to the main data array but save the insert buffer in v
		 *
		 * \param v insert buffer
		 *
		 * \param opt options
		 *
		 */
		template<typename ... v_reduce>
		void flush_v(vector<aggregate<Ti>,Memory,layout_base,grow_p> & vct_add_index_cont_0,
				     gpu::ofp_context_t & context,
				     flush_type opt = FLUSH_ON_HOST,
				     int i = 0)
		{
			// Eliminate background
			vct_data.resize(vct_index.size());

			if (opt & flush_type::FLUSH_ON_DEVICE)
			{this->flush_on_gpu<v_reduce ... >(vct_add_index_cont_0,vct_add_index_cont_1,vct_add_data_reord,context,i);}
			else
			{this->flush_on_cpu<v_reduce ... >();}

			resetBck();
		}

		/*! \brief merge the added element to the main data array but save the insert buffer in v
		 *
		 * \param v insert buffer
		 *
		 * \param opt options
		 *
		 */
		template<typename ... v_reduce>
		void flush_vd(vector<T,Memory,layout_base,grow_p> & vct_add_data_reord,
				     gpu::ofp_context_t & context,
				     flush_type opt = FLUSH_ON_HOST)
		{
			// Eliminate background
			vct_data.resize(vct_index.size());

			if (opt & flush_type::FLUSH_ON_DEVICE)
			{this->flush_on_gpu<v_reduce ... >(vct_add_index_cont_0,vct_add_index_cont_1,vct_add_data_reord,context);}
			else
			{this->flush_on_cpu<v_reduce ... >();}

			resetBck();
		}

		/*! \brief merge the added element to the main data array
		 *
		 * \param opt options
		 *
		 */
		template<typename ... v_reduce>
		void flush(gpu::ofp_context_t & context, flush_type opt = FLUSH_ON_HOST)
		{
			// Eliminate background
			vct_data.resize(vct_index.size());

			if (opt & flush_type::FLUSH_ON_DEVICE)
			{this->flush_on_gpu<v_reduce ... >(vct_add_index_cont_0,vct_add_index_cont_1,vct_add_data_reord,context);}
			else
			{this->flush_on_cpu<v_reduce ... >();}

			resetBck();
		}

		/*! \brief merge the added element to the main data array
		 *
		 * \param opt options
		 *
		 */
		void flush_remove(gpu::ofp_context_t & context, flush_type opt = FLUSH_ON_HOST)
		{
			vct_data.resize(vct_data.size()-1);

			if (opt & flush_type::FLUSH_ON_DEVICE)
			{this->flush_on_gpu_remove(context);}
			else
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " error, flush_remove on CPU has not implemented yet";
			}

			resetBck();
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
		vector<aggregate<Ti>,Memory,layout_base,grow_p> &
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

        /*! \brief Transfer from host to device
         *
         * \tparam set of parameters to transfer to device
         *
         */
        template<unsigned int ... prp>
        void hostToDevice()
        {
            vct_index.template hostToDevice<0>();
            vct_data.template hostToDevice<prp...>();
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

		/*! \brief In case we manually set the added index buffer and the add data buffer we have to call this
		 *         function before flush
		 *
		 *
		 */
		void preFlush()
		{
#ifdef __NVCC__
			vct_nadd_index.resize(vct_add_index.size());

			if (vct_nadd_index.size() != 0)
			{
				auto ite = vct_nadd_index.getGPUIterator();
				CUDA_LAUNCH((set_one_insert_buffer),ite,vct_nadd_index.toKernel());
			}
			n_gpu_add_block_slot = 1;
#endif
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

			// re-add background
			vct_data.resize(vct_data.size()+1);
			vct_data.get(vct_data.size()-1) = bck;

			htoD<decltype(vct_data)> trf(vct_data,vct_data.size()-1);
			boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(trf);

			max_ele = 0;
			n_gpu_add_block_slot = 0;
			n_gpu_rem_block_slot = 0;
		}

		void swap(vector_sparse<T,Ti,Memory,layout,layout_base,grow_p,impl,impl2,block_functor> & sp)
		{
			vct_data.swap(sp.vct_data);
			vct_index.swap(sp.vct_index);
			vct_add_index.swap(sp.vct_add_index);
			vct_add_data.swap(sp.vct_add_data);

			size_t max_ele_ = sp.max_ele;
			sp.max_ele = max_ele;
			this->max_ele = max_ele_;
		}

		vector<T,Memory,layout_base,grow_p> & private_get_vct_add_data()
		{
			return vct_add_data;
		}

		vector<aggregate<Ti>,Memory,layout_base,grow_p> & private_get_vct_add_index()
		{
			return vct_add_index;
		}

		const vector<aggregate<Ti>,Memory,layout_base,grow_p> & private_get_vct_add_index() const
		{
			return vct_add_index;
		}

		vector<aggregate<Ti>,Memory,layout_base,grow_p> & private_get_vct_nadd_index()
		{
			return vct_nadd_index;
		}

		const vector<aggregate<Ti>,Memory,layout_base,grow_p> & private_get_vct_nadd_index() const
		{
			return vct_nadd_index;
		}

		auto getSegmentToOutMap() -> decltype(blf.get_outputMap())
		{
			return blf.get_outputMap();
		}

		auto getSegmentToOutMap() const -> decltype(blf.get_outputMap())
		{
			return blf.get_outputMap();
		}

		/*! \brief Eliminate many internal temporary buffer you can use this between flushes if you get some out of memory
		 *
		 *
		 */
		void removeUnusedBuffers()
		{
			vct_add_data.resize(0);
			vct_add_data.shrink_to_fit();

			vct_add_data.resize(0);
			vct_add_data.shrink_to_fit();

			vct_add_data_reord.resize(0);
			vct_add_data_reord.shrink_to_fit();

			vct_add_data_cont.resize(0);
			vct_add_data_cont.shrink_to_fit();

			vct_add_data_unique.resize(0);
			vct_add_data_unique.shrink_to_fit();
		}

		/* \brief Return the offsets of the segments for the merge indexes
		 *
		 *
		 */
		vector<aggregate<Ti,Ti>,Memory,layout_base,grow_p> & getSegmentToMergeIndexMap()
		{
			return vct_add_index_unique;
		}

		vector<aggregate<Ti,Ti>,Memory,layout_base,grow_p> & getSegmentToMergeIndexMap() const
		{
			return vct_add_index_unique;
		}

		/*! \brief Return the mapping vector
		 *
		 * When we add new elements this vector contain the merged old elements and new elements position
		 *
		 * For example the old vector contain
		 *
		 * Old: 5 10 35 50 66 79 (6 elements)
		 * New: 7 44 7 9 44      (5 elements)  (in order are 7 7 9 44 44)
		 *
		 * The merged indexes are (when reordered)
		 *
		 * 5 7 7 9 10 35 44 44 50 66 79
		 *
		 * The returned map contain 5 elements indicating the position of the reordered elements:
		 *
		 *  0  2  3   1   4
		 * (7)(7)(9)(44)(44)
		 */
		vector<aggregate<Ti>,Memory,layout_base,grow_p> & getMappingVector()
		{
			return vct_add_index_cont_1;
		}

		/*! \brief Return the merge mapping vector
		 *
		 * When we add new elements this vector contain the merged old elements and new elements position
		 *
		 * For example the old vector contain
		 *
		 * Old: 5 10 35 50 66 79 (6 elements)
		 * New: 7 44 7 9 44      (5 elements)  (in order are 7 7 9 44 44)
		 *
		 * The merged indexes are (when reordered)
		 *
		 * 5 7 7 9 10 35 44 44 50 66 79
		 *
		 * The returned map contain 5 elements indicating the position of the reordered elements:
		 *
		 *  0  6  7  8   1   2   9  10   3   4   5
		 * (5)(7)(7)(9)(10)(35)(44)(44)(50)(66)(79)
		 */
		vector<aggregate<Ti>,Memory,layout_base,grow_p> & getMergeIndexMapVector()
		{
			return vct_index_tmp2;
		}
	};


	template<typename T, unsigned int blockSwitch = VECTOR_SPARSE_STANDARD, typename block_functor = stub_block_functor, typename indexT = int>
	using vector_sparse_gpu = openfpm::vector_sparse<
	        T,
	        indexT,
	        CudaMemory,
	        typename memory_traits_inte<T>::type,
	        memory_traits_inte,
            grow_policy_double,
            vect_isel<T>::value,
            blockSwitch,
            block_functor
            >;

	template<typename T, typename block_functor = stub_block_functor, typename indexT = long int>
	using vector_sparse_gpu_block = openfpm::vector_sparse<
	        T,
	        indexT,
	        CudaMemory,
	        typename memory_traits_inte<T>::type,
	        memory_traits_inte,
            grow_policy_double,
            vect_isel<T>::value,
            VECTOR_SPARSE_BLOCK,
            block_functor
            >;
}



#endif /* MAP_VECTOR_SPARSE_HPP_ */
