/*
 * map_vector.hpp
 *
 *  Created on: Aug 30, 2014
 *      Author: Pietro Incardona
 */

#ifndef MAP_VECTOR_HPP
#define MAP_VECTOR_HPP

#include <iostream>
#include <typeinfo>
#include "util/common.hpp"
#include "memory/PtrMemory.hpp"
#include "util/object_util.hpp"
#include "Grid/util.hpp"
#include "Vector/util.hpp"
#include "Vector/map_vector_grow_p.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "util/util_debug.hpp"
#include "util/Pack_stat.hpp"
#include "Grid/map_grid.hpp"
#include "memory/HeapMemory.hpp"
#include "vect_isel.hpp"
#include "util/object_s_di.hpp"
#include "util.hpp"
#include "util/Pack_stat.hpp"
#include "memory/ExtPreAlloc.hpp"
#include <string.h>
#include "Packer_Unpacker/Unpacker.hpp"
#include "Packer_Unpacker/Packer.hpp"
#include <fstream>
#include "Packer_Unpacker/Packer_util.hpp"
#include "Packer_Unpacker/has_pack_agg.hpp"
#include "timer.hpp"
#include "map_vector_std_util.hpp"
#include "data_type/aggregate.hpp"
#include "vector_map_iterator.hpp"
#include "util/cuda_util.hpp"
#include "cuda/map_vector_cuda_ker.cuh"
#include "map_vector_printers.hpp"

namespace openfpm
{

	template<bool active>
	struct copy_two_vectors_activate_impl
	{
		template<typename vector_type1, typename vector_type2>
		static void copy(vector_type1 & v1, vector_type2 & v2)
		{

		}

		template<typename vector_type1, typename vector_type2>
		static void copy2(vector_type1 & v1, vector_type2 & v2)
		{

		}
	};

	template<>
	struct copy_two_vectors_activate_impl<true>
	{
		template<typename vector_type1, typename vector_type2>
		static void copy(vector_type1 & v1, vector_type2 & v2)
		{
#ifdef __NVCC__
			if (v1.size() != 0)
			{
				auto it = v1.getGPUIterator();
				CUDA_LAUNCH(copy_two_vectors,it,v1.toKernel(),v2.toKernel());
			}
#endif
		}

		template<typename vector_type1, typename vector_type2>
		static void copy2(vector_type1 & v1, vector_type2 & v2)
		{
#ifdef __NVCC__
			if (v2.size() != 0)
			{
				auto it = v1.getGPUIterator();
				CUDA_LAUNCH(copy_two_vectors,it,v1.toKernel(),v2.toKernel());
			}
#endif
		}
	};

	template<bool is_ok_cuda,typename T, typename Memory,
			 template<typename> class layout_base,
			 typename grow_p>
	struct add_prp_device_impl
	{
		template <typename S,
				  typename M,
				  typename gp,
				  unsigned int impl,
				  template <typename> class layout_base2,
				  unsigned int ...args>
		static void run(openfpm::vector<T,Memory,layout_base,grow_p,impl> & this_ ,const openfpm::vector<S,M,layout_base2,gp,impl> & v)
		{
			std::cout << __FILE__ << ":" << __LINE__ << " Error the function add_prp_device only work with cuda enabled vector" << std::endl;
		}
	};

	template<bool is_ok_cuda,typename T, typename Memory,
			 template<typename> class layout_base,
			 typename grow_p>
	struct merge_prp_device_impl
	{
		template <typename S,
				  typename M,
				  typename gp,
				  unsigned int impl,
				  template <typename> class layout_base2,
				  unsigned int ...args>
		static void run(openfpm::vector<T,Memory,layout_base,grow_p,impl> & this_ ,
				        const openfpm::vector<S,M,layout_base2,gp,impl> & v,
				        unsigned int offset)
		{
			std::cout << __FILE__ << ":" << __LINE__ << " Error the function merge_prp_device only work with cuda enabled vector" << std::endl;
		}
	};

	template<typename T, typename Memory,
			 template<typename> class layout_base,
			 typename grow_p>
	struct add_prp_device_impl<true,T,Memory,layout_base,grow_p>
	{
		template <typename S,
				  typename M,
				  typename gp,
				  unsigned int impl,
				  template <typename> class layout_base2,
				  unsigned int ...args>
		static void run(vector<T,Memory,layout_base,grow_p,impl> & this_ ,const vector<S,M,layout_base2,gp,impl> & v)
		{
				// merge the data on device

	#if defined(CUDA_GPU) && defined(__NVCC__)

				size_t old_sz = this_.size();
				this_.resize(this_.size() + v.size(),DATA_ON_DEVICE);

				auto ite = v.getGPUIterator();

				CUDA_LAUNCH((merge_add_prp_device_impl<decltype(v.toKernel()),decltype(this_.toKernel()),args...>),ite,v.toKernel(),this_.toKernel(),(unsigned int)old_sz);

	#else
				std::cout << __FILE__ << ":" << __LINE__ << " Error the function add_prp_device only work when map_vector is compiled with nvcc" << std::endl;
	#endif
		}
	};

	template<typename T, typename Memory,
			 template<typename> class layout_base,
			 typename grow_p>
	struct merge_prp_device_impl<true,T,Memory,layout_base,grow_p>
	{
		template <typename S,
				  typename M,
				  typename gp,
				  unsigned int impl,
				  template <typename> class layout_base2,
				  unsigned int ...args>
		static void run(vector<T,Memory,layout_base,grow_p,impl> & this_ ,
					    const vector<S,M,layout_base2,gp,impl> & v,
					    unsigned int offset)
		{
				// merge the data on device

	#if defined(CUDA_GPU) && defined(__NVCC__)

				auto ite = v.getGPUIterator();

				CUDA_LAUNCH((merge_add_prp_device_impl<decltype(v.toKernel()),decltype(this_.toKernel()),args...>),ite,v.toKernel(),this_.toKernel(),(unsigned int)offset);

	#else
				std::cout << __FILE__ << ":" << __LINE__ << " Error the function merge_prp_device only work when map_vector is compiled with nvcc" << std::endl;
	#endif
		}
	};


	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * Stub object look at the various implementations
	 *
	 * \tparam T type of object the vector store
	 * \tparam Memory allocator to use
	 * \tparam layout layout to use
	 * \tparam grow_p grow policy for vector in case of reallocation
	 *
	 * \see vector<T,HeapMemory,grow_policy_double,STD_VECTOR>
	 * \see vector<T,Memory,grow_p,OPENFPM_NATIVE>
	 * \see vector<T,Memory,grow_p,OPENFPM_NATIVE>
	 *
	 */
	template<typename T, typename Memory, template<typename> class layout_base, typename grow_p, unsigned int impl>
	class vector
	{
		/*! \brief Stub size
		 *
		 * Report an error
		 *
		 * \return 0
		 *
		 */
		size_t size()
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error stub vector created" << std::endl;
			return 0;
		}
	};

	#include "map_vector_std.hpp"
	#include "map_vector_std_ptr.hpp"

#ifdef CUDA_GPU
	#include "cuda/map_vector_std_cuda.hpp"
#endif

	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * The layout is memory_traits_lin
	 *
	 * ### Add and access elements
	 * \snippet vector_test_util.hpp Create add and access
	 *
	 * \tparam T type of object the vector store
	 * \tparam Memory allocator to use
	 * \tparam layout layout to use
	 * \tparam grow_p grow policy for vector in case of reallocation
	 *
	 * OPENFPM_NATIVE implementation
	 *
	 */
	template<typename T,typename Memory, template <typename> class layout_base, typename grow_p>
	class vector<T,Memory,layout_base,grow_p,OPENFPM_NATIVE>
	{
		typedef vector<T,Memory,layout_base,grow_p,OPENFPM_NATIVE> self_type;

		//! Actual size of the vector, warning: it is not the space allocated in grid
		//! grid size increase by a fixed amount every time we need a vector bigger than
		//! the actually allocated space
		size_t v_size;

		//! 1-D static grid
		grid_base<1,T,Memory,typename layout_base<T>::type> base;

		/*! \brief If the argument is zero return 1 otherwise return the argument
		 *
		 * \param sz output
		 * \param arg argument
		 *
		 */
		void non_zero_one(size_t sz[1], size_t arg)
		{
			if (arg == 0)
			{sz[0] = 1;}
			else
			{sz[0] = arg;}
		}

#ifdef SE_CLASS1

		/*! \brief Check that id is not bigger than the vector size
		 *
		 * \param id element id
		 *
		 */

		void check_overflow(size_t id) const
		{
			if (id >= v_size)
			{
				std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " overflow id: " << id << "\n";
				ACTION_ON_ERROR(VECTOR_ERROR_OBJECT);
			}
		}

#endif

	public:

		//! it define that it is a vector
		typedef int yes_i_am_vector;

		//! it define that it is a vector
		typedef int yes_i_am_vector_native;

		//! Type of the encapsulation memory parameter
		typedef typename layout_base<T>::type layout_type;

		//! Type of the encapsulation memory parameter
		typedef layout_base<T> layout_base_;

		//! iterator for the vector
		typedef vector_key_iterator iterator_key;

		//! Object container for T, it is the return type of get_o it return a object type trough
		// you can access all the properties of T
		typedef typename grid_base<1,T,Memory,typename layout_base<T>::type>::container container;

		//! Type of the value the vector is storing
		typedef T value_type;

		//! Type of memory this vector use
		typedef Memory Memory_type;

		//! growing policy of this vector
		typedef grow_p grow_policy;

		template<typename Tobj>
		struct layout_base__
		{
			typedef layout_base<Tobj> type;
		};

		// Implementation of packer and unpacker for vector
#include "vector_pack_unpack.ipp"

		/*! \brief Return the size of the vector
		 *
		 * \return the size
		 *
		 */
		size_t size() const
		{
			return v_size;
		}

        /*! \brief Return the size of the vector
        *
        * \return the size
        *
        */ //remove host device
        size_t size_local() const
        {
            return v_size;
        }

		/*! \brief return the maximum capacity of the vector before reallocation
		 *
		 * \return the capacity of the vector
		 *
		 */

		size_t capacity()
		{
			return base.size();
		}

		/*! \brief Reserve slots in the vector to avoid reallocation
		 *
		 * Reserve slots in the vector to avoid reallocation
		 *
		 * \param sp number of slot to reserve
		 *
		 */

		void reserve(size_t sp)
		{
			if (sp > base.size())
			{
				//! Resize the memory
				size_t sz[1] = {sp};
				base.resize(sz);
			}

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief Clear the vector
		 *
		 * Eliminate all the elements for from the vector
		 *
		 */
		void clear()
		{
			resize(0);
		}

		/*! \brief Clear the vector
		 *
		 * Eliminate all the elements for from the vector
		 *
		 */
		void shrink_to_fit()
		{
			size_t sz[1] = {size()};
			base.resize(sz);

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief Resize the vector
		 *
		 * Resize the vector and allocate n elements
		 *
		 * \param opt options for resize. In case we know that the data are only on device memory we can use DATA_ONLY_DEVICE,
		 *                                In case we know that the data are only on host memory we can use DATA_ONLY_HOST
		 *
		 * \param blockSize The default is equal to 1. In case of accelerator buffer resize indicate the size of the block of
		 *        threads ( this is used in case of a vector of blocks where the block object override to operator=
		 *                  to distribute threads on each block element )
		 *
		 */
		void resize(size_t slot, size_t opt = DATA_ON_DEVICE | DATA_ON_HOST, unsigned int blockSize = 1)
		{
			// If we need more space than what we allocated, allocate new memory

			if (slot > base.size())
			{
				size_t gr = slot;
				// If you increase by one we smartly resize the internal capacity more than 1
				// This is to make faster patterns like resize(size()+1)
				if (slot - base.size() == 1 && opt && (opt & EXACT_RESIZE) == 0)
				{
					gr = grow_p::grow(base.size(),slot);
				}

				//! Resize the memory
				size_t sz[1] = {gr};

				base.resize(sz,opt,blockSize);
			}

			// update the vector size
			v_size = slot;

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}


		/*! \brief Resize the vector ()
		 *
		 * Resize the vector and allocate n elements
		 *
		 * \param slot number of elements
		 * \param opt options
		 *
		 */
		void resize_no_device(size_t slot)
		{
			// If we need more space than what we allocated, allocate new memory

			if (slot > base.size())
			{
				size_t gr = grow_p::grow(base.size(),slot);

				//! Resize the memory
				size_t sz[1] = {gr};
				base.resize_no_device(sz);
			}

			// update the vector size
			v_size = slot;

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		//! Access key for the vector
		typedef size_t access_key;

		/*! \brief It insert a new emtpy object on the vector, eventually it reallocate the grid
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 */
		void add()
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.resize(sz);
			}

			//! increase the vector size
			v_size++;

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief It insert a new emtpy object on the vector, eventually it reallocate the grid
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 */
		void add_no_device()
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.resize_no_device(sz);
			}

			//! increase the vector size
			v_size++;

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * \param v element to add
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 */
		void add(const T & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.resize(sz);
			}

			//! copy the element
			base.set(v_size,v);

			//! increase the vector size
			v_size++;

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief It insert a new object on the vector, eventually it reallocate the vector
		 *
		 * \param v object (encapsulated)
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 *
		 */
		void add(const typename grid_base<1,T,Memory,typename layout_base<T>::type>::container & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.resize(sz);
			}

			//! copy the added element
			base.set(v_size,v);

			//! increase the vector size
			v_size++;

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief It add the element of another vector to this vector
		 *
		 * \param v from where to take the vector
		 *
		 */
		template <typename M, typename gp> void add(const vector<T, M, layout_base,gp,OPENFPM_NATIVE> & v)
		{
			//! Add the element of v
			for (size_t i = 0 ; i < v.size() ; i++)
				add(v.get(i));
		}

		/*! \brief It merge the elements of a source vector to this vector
		 *
		 * Given 2 vector v1 and v2 of size 7,3. and as merging operation the function add.
		 * Merging the second vector v2 to
		 * the first one v1 starting from the element 2. Mean
		 *
		 * \verbatim
		 *
		 * 6   8  3   2  1   0  3    v1 elements
		 *        |   |  |
		 *       op  op  op
		 *        |   |  |
		 *        5   1  9           v2 elements
		 *
		 *-------------------------------------
		 * 6   8  8   3  10  0   3   updated v1 elements
		 *
		 * This operation is done for each selected property in args
		 *
		 * \endverbatim
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam op merging operation
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 * \param start index from where to start the merging
		 *
		 */
		template <template<typename,typename> class op, typename S, typename M, typename gp, unsigned int ...args>
		void merge_prp(const vector<S,M,layout_base,gp,OPENFPM_NATIVE> & v,
				 	   const openfpm::vector<size_t> & opart)
		{
#ifdef SE_CLASS1

			if (v.size() != opart.size())
				std::cerr << __FILE__ << ":" << __LINE__ << " error merge_prp: v.size()=" << v.size() << " must be the same as o_part.size()" << opart.size() << std::endl;

#endif
			//! Add the element of v
			for (size_t i = 0 ; i < v.size() ; i++)
			{
#ifdef SE_CLASS1

				if (opart.get(i) > size())
					std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " try to access element " << opart.get(i) << " but the vector has size " << size() << std::endl;

#endif
				// write the object in the last element
				object_s_di_op<op,decltype(v.get(i)),decltype(get(size()-1)),OBJ_ENCAP,args...>(v.get(i),get(opart.get(i)));
			}
		}

		/*! \brief It merge the elements of a source vector to this vector (on device)
		 *
		 * Given 2 vector v1 and v2 of size 7,3. and as merging operation the function add.
		 * Merging the second vector v2 to
		 * the first one v1 starting from the element 2. Mean
		 *
		 * \verbatim
		 *
		 * 6   8  3   2  1   0  3    v1 elements
		 *        |   |  |
		 *       op  op  op
		 *        |   |  |
		 *        5   1  9           v2 elements
		 *
		 *-------------------------------------
		 * 6   8  8   3  10  0   3   updated v1 elements
		 *
		 * This operation is done for each selected property in args
		 *
		 * \endverbatim
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam op merging operation
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 * \param start index from where to start the merging
		 *
		 */
		template <template<typename,typename> class op, typename S, typename M, typename gp, unsigned int ...args>
		void merge_prp_device(const vector<S,M,layout_base,gp,OPENFPM_NATIVE> & v,
				 	   unsigned int start)
		{
			merge_prp_device_impl<std::is_same<Memory,CudaMemory>::value,T,Memory,layout_base,grow_p>
			::template run<S,M,gp,OPENFPM_NATIVE,layout_base,args...>(*this,v,start);
		}


		/*! \brief It merge the elements of a source vector to this vector
		 *
		 * Given 2 vector v1 and v2 of size 7,3. and as merging operation the function add.
		 * Merging the second vector v2 to
		 * the first one v1 starting from the element 2. Mean
		 *
		 * \verbarim
		 *
		 * 6   8  3   2  1   0  3    v1 elements
		 *        |   |  |
		 *       op  op  op
		 *        |   |  |
		 *        5   1  9           v2 elements
		 *
		 *-------------------------------------
		 * 6   8  8   3  10  0   3   updated v1 elements
		 *
		 * This operation is done for each selected property in args
		 *
		 * \endverbatim
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam op merging operation
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 * \param start index from where to start the merging
		 *
		 */
		template <template<typename,typename> class op,
		          typename S,
				  typename M,
				  typename gp,
				  template <typename> class layout_base2,
				  typename vector_opart_type,
				  unsigned int ...args>
		void merge_prp_v(const vector<S,M,layout_base2,gp,OPENFPM_NATIVE> & v,
						 const vector_opart_type & opart)
		{
#ifdef SE_CLASS1

			if (v.size() != opart.size())
				std::cerr << __FILE__ << ":" << __LINE__ << " error merge_prp: v.size()=" << v.size() << " must be the same as o_part.size()" << opart.size() << std::endl;

#endif
			//! Add the element of v
			for (size_t i = 0 ; i < v.size() ; i++)
			{
#ifdef SE_CLASS1

				if (i >= opart.size())
					std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " try to access element " << opart.template get<0>(i) << " but the vector has size " << size() << std::endl;

#endif
				// write the object in the last element
				object_s_di_op<op,decltype(v.get(i)),decltype(get(size()-1)),OBJ_ENCAP,args...>(v.get(i),get(opart.template get<0>(i)));
			}
		}

		/*! \brief It merge the elements of a source vector to this vector
		 *
		 * Given 2 vector v1 and v2 of size 7,3. and as merging operation the function add.
		 * Merging the second vector v2 to
		 * the first one v1 starting from the element 2. Mean
		 *
		 * \verbarim
		 *
		 * 6   8  3   2  1   0  3    v1 elements
		 *        |   |  |
		 *       op  op  op
		 *        |   |  |
		 *        5   1  9           v2 elements
		 *
		 *-------------------------------------
		 * 6   8  8   3  10  0   3   updated v1 elements
		 *
		 * This operation is done for each selected property in args
		 *
		 * \endverbatim
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam op merging operation
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 * \param offset offset from where to copy in v
		 * \param start index from where to start the merging
		 *
		 */
		template <template<typename,typename> class op,
		          typename S,
				  typename M,
				  typename gp,
				  template <typename> class layout_base2,
				  typename vector_opart_type,
				  unsigned int ...args>
		void merge_prp_v(const vector<S,M,layout_base2,gp,OPENFPM_NATIVE> & v,
				         unsigned int offset,
						 const vector_opart_type & opart)
		{
			size_t i2 = 0;

			for (size_t i = offset ; i < v.size() ; i++)
			{
				auto dst = v.get(opart.template get<0>(i2));
				auto src = v.get(i);
				copy_cpu_encap_encap_op_prp<op,decltype(v.get(0)),decltype(v.get(0)),args...> cp(src,dst);

				boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(args)> >(cp);
				i2++;
			}
		}

		/*! \brief It merge the elements of a source vector to this vector
		 *
		 * Given 2 vector v1 and v2 of size 7,3. and as merging operation the function add.
		 * Merging the second vector v2 to
		 * the first one v1 starting from the element 2. Mean
		 *
		 * \verbarim
		 *
		 * 6   8  3   2  1   0  3    v1 elements
		 *        |   |  |
		 *       op  op  op
		 *        |   |  |
		 *        5   1  9           v2 elements
		 *
		 *-------------------------------------
		 * 6   8  8   3  10  0   3   updated v1 elements
		 *
		 * This operation is done for each selected property in args
		 *
		 * \endverbatim
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam op merging operation
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 * \param opart merging indexes (property 1)
		 * \param start starting merging index for opart
		 * \param stop stop merging index for opart
		 *
		 */
		template <template<typename,typename> class op,
		          typename S,
				  typename M,
				  typename gp,
				  template <typename> class layout_base2,
				  typename vector_opart_type,
				  unsigned int ...args>
		void merge_prp_v_device(const vector<S,M,layout_base2,gp,OPENFPM_NATIVE> & v,
						 const vector_opart_type & opart,
						 unsigned int start,
						 unsigned int stop)
		{
#ifdef SE_CLASS1

			if (v.size() != stop - start)
				std::cerr << __FILE__ << ":" << __LINE__ << " error merge_prp: v.size()=" << v.size() << " must be the same as stop - start" << stop - start << std::endl;

#endif

#ifdef __NVCC__

			size_t sz[1] = {stop - start};
			grid_sm<1,void> nm(sz);

			auto ite = nm.getGPUIterator();

			// write the object in the last element
			CUDA_LAUNCH((merge_add_prp_device_impl_src_dst_opar_offset<op,
					                                                   decltype(v.toKernel()),
																	   decltype(this->toKernel()),
			                                                           decltype(opart.toKernel()),
																	   args...>),ite,v.toKernel(),this->toKernel(),opart.toKernel(),start);

			// calculate
#else
			std::cout << __FILE__ << ":" << __LINE__ << " Error you have to compile map_vector.hpp with nvcc to make GPU code working" << std::endl;

#endif
		}

		/*! \brief It merge the elements of a source vector to this vector
		 *
		 * Given 2 vector v1 and v2 of size 7,3. and as merging operation the function add.
		 * Merging the second vector v2 to
		 * the first one v1 starting from the element 2. Mean
		 *
		 * \verbarim
		 *
		 * 6   8  3   2  1   0  3    v1 elements
		 *        |   |  |
		 *       op  op  op
		 *        |   |  |
		 *        5   1  9           v2 elements
		 *
		 *-------------------------------------
		 * 6   8  8   3  10  0   3   updated v1 elements
		 *
		 * This operation is done for each selected property in args
		 *
		 * \endverbatim
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam op merging operation
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 * \param opart merging indexes (property 0)
		 * \param i starting mergong indexes
		 *
		 */
		template <template<typename,typename> class op,
		          typename S,
				  typename M,
				  typename gp,
				  template <typename> class layout_base2,
				  typename vector_opart_type,
				  unsigned int ...args>
		void merge_prp_v_device(const vector<S,M,layout_base2,gp,OPENFPM_NATIVE> & v,
						 unsigned int start,
						 const vector_opart_type & opart)
		{
#ifdef SE_CLASS1

			if (v.size() < opart.size() + start)
				std::cerr << __FILE__ << ":" << __LINE__ << " error merge_prp: v.size()=" << v.size() << " must be snaller than o_part.size() + start " << opart.size() + start << std::endl;

#endif

#ifdef __NVCC__

			auto ite = opart.getGPUIterator();

			// write the object in the last element
			CUDA_LAUNCH((merge_add_prp_device_impl_src_offset_dst_opar<op,
					                                                   decltype(v.toKernel()),
																	   decltype(this->toKernel()),
																	   decltype(opart.toKernel()),
																	   args... >),ite,v.toKernel(),this->toKernel(),opart.toKernel(),start);

			// calculate
#else
			std::cout << __FILE__ << ":" << __LINE__ << " Error you have to compile map_vector.hpp with nvcc to make GPU code working" << std::endl;

#endif
		}

		/*! \brief It merge the elements of a source vector to this vector
		 *
		 * Given 2 vector v1 and v2 of size 7,3. and as merging operation the function add.
		 * Merging the second vector v2 to
		 * the first one v1 starting from the element 2. Mean
		 *
		 * \verbarim
		 *
		 * 6   8  3   2  1   0  3    v1 elements
		 *        |   |  |
		 *       op  op  op
		 *        |   |  |
		 *        5   1  9           v2 elements
		 *
		 *-------------------------------------
		 * 6   8  8   3  10  0   3   updated v1 elements
		 *
		 * This operation is done for each selected property in args
		 *
		 * \endverbatim
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam op merging operation
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 * \param start index from where to start the merging
		 *
		 */
		template <template<typename,typename> class op,
		          typename S,
				  typename M,
				  typename gp,
				  template <typename> class layout_base2,
				  unsigned int ...args>
		void merge_prp_v(const vector<S,M,layout_base2,gp,OPENFPM_NATIVE> & v,
				         size_t start)
		{
			//! Add the element of v
			for (size_t i = 0 ; i < v.size() ; i++)
			{
#ifdef SE_CLASS1

				if (start + i >= v_size)
					std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " try to access element " << start+i << " but the vector has size " << size() << std::endl;

#endif
				// write the object in the last element
				object_s_di_op<op,decltype(v.get(0)),decltype(get(0)),OBJ_ENCAP,args...>(v.get(i),get(start+i));
			}
		}

		/*! \brief It add the element of a source vector to this vector
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 *
		 */
		template <typename S,
		          typename M,
				  typename gp,
				  unsigned int impl,
				  template <typename> class layout_base2,
				  unsigned int ...args>
		void add_prp(const vector<S,M,layout_base2,gp,impl> & v)
		{
			//! Add the element of v
			for (size_t i = 0 ; i < v.size() ; i++)
			{
				// Add a new element
				add();

				// write the object in the last element
				object_s_di<decltype(v.get(i)),decltype(get(size()-1)),OBJ_ENCAP,args...>(v.get(i),get(size()-1));
			}
		}

		/*! \brief It add the element of a source vector to this vector
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 *
		 */
		template <typename S,
		          typename M,
				  typename gp,
				  unsigned int impl,
				  template <typename> class layout_base2,
				  unsigned int ...args>
		void add_prp_device(const vector<S,M,layout_base2,gp,impl> & v)
		{
			add_prp_device_impl<std::is_same<Memory,CudaMemory>::value,T,Memory,layout_base,grow_p>
			::template run<S,M,gp,impl,layout_base2,args...>(*this,v);
		}

		/*! \brief Insert an entry in the vector
		 *
		 * \size_t key Where to insert the element
		 *
		 */
		void insert(size_t key)
		{
			add();

			long int d_k = (long int)size()-1;
			long int s_k = (long int)size()-2;

			// keys
			while (s_k >= (long int)key)
			{
				set(d_k,get(s_k));
				d_k--;
				s_k--;
			}
		}


		/*! \brief Remove one entry from the vector
		 *
		 * \param key element to remove
		 *
		 */
		void remove(size_t key)
		{
			size_t d_k = key;
			size_t s_k = key + 1;

			// keys
			while (s_k < size())
			{
				set(d_k,get(s_k));
				d_k++;
				s_k++;
			}

			// re-calculate the vector size

			v_size--;
		}

		/*! \brief Remove several entries from the vector
		 *
		 * \warning the keys in the vector MUST be sorted
		 *
		 * \param keys objects id to remove
		 * \param start key starting point
		 *
		 */
		void remove(openfpm::vector<size_t> & keys, size_t start = 0)
		{
			// Nothing to remove return
			if (keys.size() <= start )
				return;

			size_t a_key = start;
			size_t d_k = keys.get(a_key);
			size_t s_k = keys.get(a_key) + 1;

			// keys
			while (s_k < size())
			{
				// s_k should always point to a key that is not going to be deleted
				while (a_key+1 < keys.size() && s_k == keys.get(a_key+1))
				{
					a_key++;
					s_k = keys.get(a_key) + 1;
				}

				// In case of overflow
				if (s_k >= size())
					break;

				set(d_k,get(s_k));
				d_k++;
				s_k++;
			}

			// re-calculate the vector size

			v_size -= keys.size() - start;
		}

		/*! \brief Remove several entries from the vector
		 *
		 * \warning the keys in the vector MUST be sorted
		 *
		 * \param keys objects id to remove
		 * \param start key starting point
		 *
		 */
		void remove(openfpm::vector<aggregate<int>> & keys, size_t start = 0)
		{
			// Nothing to remove return
			if (keys.size() <= start )
				return;

			size_t a_key = start;
			size_t d_k = keys.template get<0>(a_key);
			size_t s_k = keys.template get<0>(a_key) + 1;

			// keys
			while (s_k < size())
			{
				// s_k should always point to a key that is not going to be deleted
				while (a_key+1 < keys.size() && s_k == keys.template get<0>(a_key+1))
				{
					a_key++;
					s_k = keys.template get<0>(a_key) + 1;
				}

				// In case of overflow
				if (s_k >= size())
					break;

				set(d_k,get(s_k));
				d_k++;
				s_k++;
			}

			// re-calculate the vector size

			v_size -= keys.size() - start;
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
		inline auto get(size_t id) const -> decltype(base.template get<p>(grid_key_dx<1>(0)))
		{
#if defined(SE_CLASS1) && !defined(__NVCC__)
			check_overflow(id);
#endif
			grid_key_dx<1> key(id);


			return base.template get<p>(key);
		}

		/*! \brief Indicate that this class is not a subset
		 *
		 * \return false
		 *
		 */
		bool isSubset() const
		{
			return false;
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */
		inline auto get(size_t id) -> decltype(base.get_o(grid_key_dx<1>(id)))
		{
#if defined(SE_CLASS1) && !defined(__NVCC__)
			check_overflow(id);
#endif
			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * \deprecated
		 *
		 * exactly as get, exist to keep the compatibility with grid
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */

		inline const typename grid_base<1,T,Memory,typename layout_base<T>::type>::container get_o(size_t id) const
		{
#if defined(SE_CLASS1) && !defined(__NVCC__)
			check_overflow(id);
#endif
			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Fill the buffer with a byte
		 *
		 * \param c char to fill the buffer with
		 *
		 */
		template<unsigned int id> void fill(unsigned char c)
		{
			base.template fill<id>(c);
		}

		/*! \brief It return the properties arrays.
		 *
		 * In case of Cuda memory it return the device pointers to pass to the kernels
		 *
		 *
		 */
		template<unsigned int id> void * getDeviceBufferCopy()
		{
			return base.template getDeviceBuffer<id>();
		}

		/*! \brief It return the properties arrays.
		 *
		 * In case of Cuda memory it return the device pointers to pass to the kernels
		 *
		 * This variant does not copy the host memory to the device memory
		 *
		 *
		 */
		template<unsigned int id> void * getDeviceBuffer()
		{
			return base.template getDeviceBuffer<id>();
		}


		/*! \brief Get the last element of the vector
		 *
		 * \return the last element (encapsulated)
		 *
		 */
		inline const typename grid_base<1,T,Memory,layout_type>::container last() const
		{
			grid_key_dx<1> key(size()-1);

			return base.get_o(key);
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
         */ //remove device host
        template <unsigned int p>
        inline auto getProp(const unsigned int & id) -> decltype(base.template get<p>(grid_key_dx<1>(0)))
        {   //uncomment this
            return this->template get<p>(id);
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
		 *///remove host device
        template <unsigned int p,typename KeyType>
        inline auto getProp(const KeyType & id) -> decltype(base.template get<p>(grid_key_dx<1>(0)))
		{
            //uncomment this
            return this->template get<p>(id.getKey());
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
		 */ //remove device host
		template <unsigned int p, typename keyType>
         inline auto getProp(const keyType & id) const -> decltype(base.template get<p>(grid_key_dx<1>(0)))
		{   //uncomment this
			return this->template get<p>(id.getKey());
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
		inline auto get(size_t id) -> decltype(base.template get<p>(grid_key_dx<1>(0)))
		{
#if defined(SE_CLASS1) && !defined(__NVCC__)
			check_overflow(id);
#endif
			grid_key_dx<1> key(id);

			return base.template get<p>(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */
		inline auto get(size_t id) const -> const decltype(base.get_o(grid_key_dx<1>(id)))
		{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
#if defined(SE_CLASS1) && !defined(__NVCC__)
			check_overflow(id);
#endif
			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Get the last element of the vector
		 *
		 * \return the element (encapsulated)
		 *
		 */

		inline typename grid_base<1,T,Memory,typename layout_base<T>::type >::container last()
		{
			grid_key_dx<1> key(size()-1);

			return base.get_o(key);
		}

		//! Destructor
		~vector() THROW
		{
			// Eliminate the pointer
		}

		/*! \brief It duplicate the vector
		 *
		 * \return a duplicated vector
		 *
		 */
		vector<T, Memory,layout_base,grow_p,OPENFPM_NATIVE> duplicate() const
		{
			vector<T, Memory,layout_base,grow_p,OPENFPM_NATIVE> dup;

			dup.v_size = v_size;
			dup.base.swap(base.duplicate());

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			dup.base_gpu.constructor_impl(v_size,dup.base.toKernel());

#endif

			copy_two_vectors_activate_impl<Memory::isDeviceHostSame() == false>::copy(dup,*this);

			return dup;
		}

		/*! \brief Constructor from another temporal vector
		 *
		 * \param v the vector
		 *
		 */
		vector(vector<T, Memory,layout_base,grow_p,OPENFPM_NATIVE> && v)
		:v_size(0)
		{
			swap(v);

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief Constructor from another constant vector
		 *
		 * \param v the vector
		 *
		 */
		vector(const vector<T, Memory,layout_base,grow_p,OPENFPM_NATIVE> & v) THROW
		:v_size(0)
		{
			swap(v.duplicate());

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		//! Constructor, vector of size 0
		vector() THROW
		:v_size(0),base(0)
		{
			base.setMemory();

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		//! Constructor, vector of size sz
		vector(size_t sz) THROW
		:v_size(sz),base(sz)
		{
			base.setMemory();

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief Set the object id to obj
		 *
		 * \param id element
		 * \param obj object (encapsulated)
		 *
		 */
		void set(size_t id, const typename grid_base<1,T,Memory,typename layout_base<T>::type>::container & obj)
		{
#ifdef SE_CLASS1
			check_overflow(id);
#endif
			//! copy the element
			base.set(id,obj);
		}

		/*! \brief It set an element of the vector from a object that is a subset of the vector properties
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam encap_S object that encapsulate the object
		 * \tparam args ids of the properties to map the object to
		 *
		 * \param i element to set
		 * \param obj object that encapsulate the object
		 *
		 * \param v source vector
		 *
		 */
		template <typename encap_S, unsigned int ...args> void set_o(size_t i, const encap_S & obj)
		{
			// write the object in the last element
			object_s_di<encap_S,decltype(get(i)),OBJ_ENCAP,args...>(obj,get(i));
		}

		/*! \brief Set the object id to obj
		 *
		 * \param id
		 * \param obj
		 *
		 */
		void set(size_t id, const T & obj)
		{
#ifdef SE_CLASS1
			check_overflow(id);
#endif
			//! copy the element
			base.set(id,obj);
		}

		/*! \brief Set the element of the vector v from another element of another vector
		 *
		 * \param id element id
		 * \param v vector source
		 * \param src source element
		 *
		 */
		void set(size_t id, vector<T,Memory,layout_base,grow_p,OPENFPM_NATIVE> & v, size_t src)
		{
#ifdef SE_CLASS1
			check_overflow(id);
#endif
			base.set(id,v.base,src);
		}

		/*! \brief Assignment operator
		 *
		 * move semantic movement operator=
		 *
		 * \param mv vector
		 *
		 * \return itself
		 *
		 */
		vector<T, Memory,layout_base,grow_p,OPENFPM_NATIVE> & operator=(vector<T, Memory, layout_base,grow_p,OPENFPM_NATIVE> && mv)
		{
			v_size = mv.v_size;
			base.swap(mv.base);

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif

			return *this;
		}

		/*! \brief Assignment operator
		 *
		 * it copy
		 *
		 * \param mv vector
		 *
		 * \return itself
		 *
		 */
		vector<T, Memory,layout_base,grow_p,OPENFPM_NATIVE> & operator=(const vector<T, Memory, layout_base ,grow_p,OPENFPM_NATIVE> & mv)
		{
			v_size = mv.v_size;
			size_t rsz[1] = {v_size};
			if(rsz[0]>base.size()) {
                base.resize(rsz);
            }
			// copy the object on cpu
			for (size_t i = 0 ; i < v_size ; i++ )
			{
				grid_key_dx<1> key(i);
				base.set(key,mv.base,key);
			}

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif

			copy_two_vectors_activate_impl<Memory::isDeviceHostSame() == false>::copy2(*this,mv);

			return *this;
		}

		/*! \brief Assignment operator
		 *
		 * move semantic movement operator=
		 *
		 * \param mv vector
		 *
		 * \return itself
		 *
		 */
		template<typename Mem, typename gp> vector<T, Memory,layout_base,grow_p,OPENFPM_NATIVE> & operator=(vector<T, Mem, layout_base,gp,OPENFPM_NATIVE> && mv)
		{
			v_size = mv.v_size;
			base.swap(mv.base);

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif

			return *this;
		}

		/*! \brief Assignment operator
		 *
		 * it copy
		 *
		 * \param mv vector
		 *
		 * \return itself
		 *
		 */
		template<typename Mem, typename gp> vector<T, Memory,layout_base,grow_p,OPENFPM_NATIVE> & operator=(const vector<T, Mem, layout_base ,gp,OPENFPM_NATIVE> & mv)
		{
			v_size = mv.getInternal_v_size();
			size_t rsz[1] = {v_size};
			base.resize(rsz);

			// copy the object
			for (size_t i = 0 ; i < v_size ; i++ )
			{
				grid_key_dx<1> key(i);
				base.set(key,mv.getInternal_base(),key);
			}

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif

			copy_two_vectors_activate_impl<Memory::isDeviceHostSame() == false && Mem::isDeviceHostSame() == false>::copy2(*this,mv);

			return *this;
		}

		/*! \brief Assignment operator
		 *
		 * move semantic movement operator=
		 *
		 * \param mv vector
		 *
		 * \return itself
		 *
		 */
		template<typename Mem, template <typename> class layout_base2>
		vector<T, Memory,layout_base2,grow_p,OPENFPM_NATIVE> & operator=(vector<T, Mem, layout_base2,grow_p,OPENFPM_NATIVE> && mv)
		{
			v_size = mv.v_size;
			base.swap(mv.base);

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif

			return *this;
		}

		/*! \brief Assignment operator
		 *
		 *
		 * \param mv vector to copy
		 *
		 * \return itself
		 *
		 *
		 */
		template<typename Mem,
		         template <typename> class layout_base2,
		         typename check = typename std::enable_if<!std::is_same<typename layout_base2<T>::type,typename layout_base<T>::type>::value >::type>
		vector<T, Memory,layout_base,grow_p,OPENFPM_NATIVE> &
		operator=(const vector<T, Mem, layout_base2 ,grow_p,OPENFPM_NATIVE> & mv)
		{
			v_size = mv.getInternal_v_size();
			size_t rsz[1] = {v_size};
			base.resize(rsz);

			// copy the object
			for (size_t i = 0 ; i < v_size ; i++ )
			{
				grid_key_dx<1> key(i);
				base.set_general(key,mv.getInternal_base(),key);
			}

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif

			copy_two_vectors_activate_impl<Memory::isDeviceHostSame() == false && Mem::isDeviceHostSame() == false>::copy2(*this,mv);

			return *this;
		}

		/*! \brief Check that two vectors are equal
		 *
		 * \param vector to compare
		 *
		 */
		bool operator!=(const vector<T, Memory, layout_base,grow_p,OPENFPM_NATIVE> & v) const
		{
			return !this->operator==(v);
		}

		/*! \brief Check that two vectors are not equal
		 *
		 * \param vector to compare
		 *
		 */
		bool operator==(const vector<T, Memory, layout_base, grow_p,OPENFPM_NATIVE> & v) const
		{
			if (v_size != v.v_size)
				return false;

			// check object by object
			for (size_t i = 0 ; i < v_size ; i++ )
			{
				grid_key_dx<1> key(i);

				if (base.get_o(key) != v.base.get_o(key))
					return false;
			}

			return true;
		}

		/*! \brief Swap the memory with another vector
		 *
		 * \note this is different from the normal swap please look a grid_base_implementation.hpp method swap_nomode
		 *       for info
		 *
		 * \param v vector
		 *
		 */
		void swap_nomode(openfpm::vector<T,Memory, layout_base,grow_p,OPENFPM_NATIVE> & v)
		{
			size_t sz_sp = v_size;

			// swap the v_size
			v_size = v.v_size;

			base.swap_nomode(v.base);

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif

			v.v_size = sz_sp;
		}

		/*! \brief Swap the memory with another vector
		 *
		 * \param v vector
		 *
		 */
		void swap(openfpm::vector<T,Memory, layout_base,grow_p,OPENFPM_NATIVE> & v)
		{
			size_t sz_sp = v_size;

			// swap the v_size
			v_size = v.v_size;

			base.swap(v.base);
			v.v_size = sz_sp;


#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());
			v.base_gpu.constructor_impl(v.v_size,v.base.toKernel());

#endif
		}

		/*! \brief Swap the memory with another vector
		 *
		 * \param v vector
		 *
		 */
		void swap(openfpm::vector<T,Memory, layout_base,grow_p,OPENFPM_NATIVE> && v)
		{
			size_t sz_sp = v_size;

			// swap the v_size
			v_size = v.v_size;

			base.swap(v.base);
			v.v_size = sz_sp;

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());
			v.base_gpu.constructor_impl(v.v_size,v.base.toKernel());

#endif
		}

		/*! \brief Get iterator over the particles from a particular index
		 *
		 * \param start starting point
		 *
		 * \return an iterator to iterate from a particular index
		 *
		 */
		vector_key_iterator getIteratorFrom(size_t start) const
		{
			return vector_key_iterator(v_size,start);
		}

		/*! \brief Get iterator over the particles from 0 until a particular index
		 *
		 * \warning stop point is not included
		 *
		 * \param stop stop point
		 *
		 * \return an iterator to iterate until a particular index
		 *
		 */
		vector_key_iterator getIteratorTo(size_t stop) const
		{
			return vector_key_iterator(stop,0);
		}

#ifdef CUDA_GPU

		/*! \brief Get an iterator for the GPU
		 *
		 *
		 */
		ite_gpu<1> getGPUIteratorTo(long int stop, size_t n_thr = default_kernel_wg_threads_) const
		{
			grid_key_dx<1,long int> start(0);
			grid_key_dx<1,long int> stop_(stop);

			return base.getGPUIterator(start,stop_,n_thr);
		}


#endif

		/*! \brief Get the vector elements iterator
		 *
		 *
		 * \return an iterator to iterate through all the elements of the vector
		 *
		 */

		vector_key_iterator getDomainIterator() const
		{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
			return getIterator();
		}

		/*! \brief Get the vector elements iterator
		 *
		 *
		 * \return an iterator to iterate through all the elements of the vector
		 *
		 */

		vector_key_iterator getIterator() const
		{
			return vector_key_iterator(v_size);
		}

        /*! \brief Get the vector elements iterator
 *
 *
 * \return an iterator to iterate through all the elements of the vector
 *
 */
        template<unsigned int p>
        vector_key_iterator_ele<p,self_type> getIteratorElements() const
        {
#ifdef SE_CLASS2
            check_valid(this,8);
#endif
            return vector_key_iterator_ele<p,self_type>(*this,size());
        }

#ifdef CUDA_GPU

		/*! \brief Get an iterator for the GPU
		 *
		 *
		 */
		ite_gpu<1> getGPUIterator(size_t n_thr = default_kernel_wg_threads_) const
		{
			grid_key_dx<1> start(0);
			grid_key_dx<1> stop(size()-1);

			return base.getGPUIterator(start,stop,n_thr);
		}

        /*! \brief Get a domain iterator for the GPU
         *
         *
         */
        ite_gpu<1> getDomainIteratorGPU(size_t n_thr = default_kernel_wg_threads_) const
        {
            return getGPUIterator(n_thr);
        }

#endif
		/*! \brief Return the size of the message needed to pack this object
		 *
		 * \return The size
		 *
		 */

		size_t packObjectSize()
		{
			return base.packObjectSize();
		}

		/*! \brief Pack the object into the given pointer
		 *
		 * \param mem pointer
		 *
		 * \return the size of the packed message
		 *
		 */
		size_t packObject(void * mem)
		{
			return base.packObject(mem);
		}

		/*! \brief Calculate the memory size required to allocate n elements
		 *
		 * Calculate the total size required to store n-elements in a vector
		 *
		 * \param n number of elements
		 * \param e unused
		 *
		 * \return the size of the allocation number e
		 *
		 */
		template<int ... prp> static inline size_t calculateMem(size_t n, size_t e)
		{
			if (n == 0)
			{
				return 0;
			}
			else
			{
				if (sizeof...(prp) == 0)
					return grow_p::grow(0,n) * sizeof(typename T::type);

				typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

				return grow_p::grow(0,n) * sizeof(prp_object);
			}
		}

		/*! \brief Calculate the memory size required to pack n elements
		 *
		 * Calculate the total size required to store n-elements in a vector
		 *
		 * \param n number of elements
		 * \param e unused
		 *
		 * \return the size of the allocation number e
		 *
		 */
		template<int ... prp> static inline size_t packMem(size_t n, size_t e)
		{
			if (sizeof...(prp) == 0)
				return n * sizeof(typename T::type);

			typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

			return n * sizeof(prp_object);
		}

		/*! \brief How many allocation are required to create n-elements
		 *
		 * \param n number of elements
		 *
		 * \return the number of allocations
		 *
		 */
		inline static size_t calculateNMem(size_t n)
		{
			return 1;
		}

		/*! \brief Return the memory object
		*
		* Return the memory object
		*
		* \tparam p array to retrieve
		*
		*/
		template<unsigned int p>
		auto getMemory() -> decltype(base.template getMemory<p>())
		{
			return base.template getMemory<p>();
		}

		/*! \brief Set the memory of the base structure using an object
		 *
		 * \param mem Memory object to use for allocation
		 *
		 */
		template<unsigned int p = 0> void setMemory(Memory & mem)
		{
			base.template setMemory<p>(mem);

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief Set the memory of the base structure using an object
		 *
		 * \param mem Memory object to use for allocation
		 *
		 */
		void setMemoryArray(Memory * mem)
		{
			base.setMemoryArray(mem);

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

			base_gpu.constructor_impl(v_size,this->base.toKernel());

#endif
		}

		/*! \brief Return the pointer that store the data
		 *
		 * \tparam property from which take the pointer
		 *
		 * \return the pointer that store the data
		 *
		 */
		template<unsigned int p = 0> void * getPointer()
		{
			return base.template getPointer<p>();
		}

		/*! \brief Return the pointer that store the data
		 *
		 * \return the pointer that store the data
		 *
		 */
		template<unsigned int p = 0> const void * getPointer() const
		{
			return base.getPointer();
		}

		/*! \brief This class has pointer inside
		 *
		 * \return false
		 *
		 */
		static bool noPointers()
		{
			return false;
		}

		/*! \brief Internal function
		 *
		 * \return the size of the vector
		 *
		 */
		const size_t & getInternal_v_size() const
		{
			return v_size;
		}

		/*! \brief Internal function
		 *
		 * \return the internal 1D grid base
		 *
		 */
		const grid_base<1,T,Memory,layout_type> & getInternal_base() const
		{
			return base;
		}

		/*! \brief Copy the memory from host to device
		 *
		 *
		 */
		template<unsigned int ... prp> void hostToDevice()
		{
			base.template hostToDevice<prp ...>();
		}

		/*! \brief Synchronize the memory buffer in the device with the memory in the host
		 *
		 *
		 */
		template<unsigned int ... prp> void deviceToHost()
		{
			base.template deviceToHost<prp ...>();
		}


		/*! \brief Synchronize the memory buffer in the device with the memory in the host
		 *
		 *
		 */
		template<unsigned int ... prp> void deviceToHost(size_t start, size_t stop)
		{
			base.template deviceToHost<prp ...>(start,stop);
		}

		/*! \brief Synchronize the memory buffer in the device with the memory in the host
		 *
		 *
		 */
		template<unsigned int ... prp> void hostToDevice(size_t start, size_t stop)
		{
			base.template hostToDevice<prp ...>(start,stop);
		}

		/*! \brief Synchronize the memory buffer in the device with the memory in the host respecting NUMA domains
		 *
		 *
		 */
		template<unsigned int ... prp> void hostToDeviceNUMA(size_t start, size_t stop)
		{
			base.template hostToDeviceNUMA<prp ...>(start,stop);
		}

		/*! \brief Synchronize the memory buffer in the device with the memory in the host respecing NUMA domains
		 *
		 *
		 */
		template<unsigned int ... prp> void hostToDeviceNUMA()
		{
			base.template hostToDeviceNUMA<prp ...>();
		}

#if defined(CUDIFY_USE_SEQUENTIAL) || defined(CUDIFY_USE_OPENMP)

		mutable vector_gpu_ker<typename apply_transform<layout_base,T>::type,layout_base> base_gpu;

		/*! \brief Convert the grid into a data-structure compatible for computing into GPU
		 *
		 *  The object created can be considered like a reference of the original
		 *
		 * \return an usable vector in the kernel
		 *
		 */
		inline vector_gpu_ker_ref<typename apply_transform<layout_base,T>::type,layout_base> toKernel()
		{
//			base_gpu.constructor_impl(v_size,this->base.toKernel());

			return vector_gpu_ker_ref<typename apply_transform<layout_base,T>::type,layout_base>(base_gpu);
		}

		/*! \brief Convert the grid into a data-structure compatible for computing into GPU
		 *
		 *  The object created can be considered like a reference of the original
		 *
		 * \return an usable vector in the kernel
		 *
		 */
		inline const vector_gpu_ker_ref<typename apply_transform<layout_base,T>::type,layout_base> toKernel() const
		{
//			base_gpu.constructor_impl(v_size,this->base.toKernel());

			return vector_gpu_ker_ref<typename apply_transform<layout_base,T>::type,layout_base>(base_gpu);
		}

#else

		/*! \brief Convert the grid into a data-structure compatible for computing into GPU
		 *
		 *  The object created can be considered like a reference of the original
		 *
		 * \return an usable vector in the kernel
		 *
		 */
		vector_gpu_ker<typename apply_transform<layout_base,T>::type,layout_base> toKernel()
		{
			vector_gpu_ker<typename apply_transform<layout_base,T>::type,layout_base> v(v_size,this->base.toKernel());

			return v;
		}

		/*! \brief Convert the grid into a data-structure compatible for computing into GPU
		 *
		 *  The object created can be considered like a reference of the original
		 *
		 * \return an usable vector in the kernel
		 *
		 */
		const vector_gpu_ker<typename apply_transform<layout_base,T>::type,layout_base>  toKernel() const
		{
			if (base.size() == 0)
			{std::cout << __FILE__ << ":" << __LINE__ << " Warning you are off-loading with toGPU a vector that seem to be empty or not initialized" << std::endl; }

			vector_gpu_ker<typename apply_transform<layout_base,T>::type,layout_base> v(v_size,this->base.toKernel());

			return v;
		}

#endif

		/*! Convert this vector into a string
		 *
		 * \param prefix prefix to use for printing
		 *
		 * \return a string showing thid vector
		 *
		 */
		template<unsigned int ... prps>
		const std::string toString(std::string prefix = std::string())
		{
			std::stringstream ss;
			auto it = getIterator();

			while (it.isNext())
			{
				auto p = it.get();

				ss << prefix;

				ss << prefix <<  " element[" << p << "]" << " ";

				vector_printer<self_type,prps ...> vp(*this,p,ss);
				boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prps)>>(vp);

				ss << std::endl;

				++it;
			}

			return ss.str();
		}

		void * internal_get_size_pointer()	{return &v_size;}

		void print_size()
		{
#ifndef DISABLE_ALL_RTTI
			std::cout << "the size of: " << demangle(typeid(self_type).name()) << " is " << sizeof(self_type) << std::endl;
			std::cout << "    " << demangle(typeid(decltype(v_size)).name()) << ":" << sizeof(decltype(v_size)) << std::endl;
			std::cout << "    " << demangle(typeid(decltype(base)).name()) << ":" << sizeof(decltype(base)) << std::endl;
#endif
		}

	};

	template <typename T> using vector_std = vector<T, HeapMemory, memory_traits_lin, openfpm::grow_policy_double, STD_VECTOR>;
	template<typename T> using vector_gpu = openfpm::vector<T,CudaMemory,memory_traits_inte>;
    template<typename T> using vector_soa = openfpm::vector<T,HeapMemory,memory_traits_inte>;
	template<typename T> using vector_gpu_lin = openfpm::vector<T,CudaMemory,memory_traits_lin>;
	template<typename T> using vector_gpu_single = openfpm::vector<T,CudaMemory,memory_traits_inte,openfpm::grow_policy_identity>;
	template<typename T> using vector_custd = vector<T, CudaMemory, memory_traits_inte, openfpm::grow_policy_double, STD_VECTOR>;
}

#endif
