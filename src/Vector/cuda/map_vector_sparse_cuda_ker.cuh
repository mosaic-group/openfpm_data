/*
 * map_vector_sparse_cuda_ker.cuh
 *
 *  Created on: Jan 23, 2019
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_SPARSE_CUDA_KER_CUH_
#define MAP_VECTOR_SPARSE_CUDA_KER_CUH_

#include "util/for_each_ref.hpp"

//todo: Check where it's a good place to put the following method...
template<typename dim3Ta, typename dim3Tb>
inline __device__ __host__ int dim3CoordToInt(const dim3Ta & coord, const dim3Tb & dimensions)
{
    int res = coord.z;
    res *= dimensions.y;
    res += coord.y;
    res *= dimensions.x;
    res += coord.x;
    return res;
}
// Specialization allowing transparency
inline __device__ __host__ int dim3CoordToInt(int coord, int dimension)
{
    return coord;
}

namespace openfpm
{
	template<typename index_type>
	struct sparse_index
	{
		index_type id;
	};

#if defined(__NVCC__) && !defined(CUDA_ON_CPU)
	static __shared__ int vct_atomic_add;
	static __shared__ int vct_atomic_rem;
#endif

	template<typename T,
			 typename Ti,
			 template<typename> class layout_base>
	class vector_sparse_gpu_ker
	{
		vector_gpu_ker<aggregate<Ti>,layout_base> vct_index;

		vector_gpu_ker<T,layout_base> vct_data;

		vector_gpu_ker<aggregate<Ti>,layout_base> vct_add_index;

		vector_gpu_ker<aggregate<Ti>,layout_base> vct_rem_index;

		vector_gpu_ker<aggregate<Ti>,layout_base> vct_nadd_index;

		vector_gpu_ker<aggregate<Ti>,layout_base> vct_nrem_index;

		vector_gpu_ker<T,layout_base> vct_add_data;

		// the const is forced by the getter that only return const encap that should not allow the modification of bck
		// this should possible avoid to define an object const_encap
		//mutable vector_gpu_ker<T,layout_base> vct_data_bck;

		int nslot_add;
		int nslot_rem;

		/*! \brief get the element i
		 *
		 * search the element x
		 *
		 * \param i element i
		 */
		inline __device__ void _branchfree_search(Ti x, Ti & id) const
		{
			if (vct_index.size() == 0)	{id = 0; return;}
			const Ti *base = &vct_index.template get<0>(0);
			const Ti *end = (const Ti *)vct_index.template getPointer<0>() + vct_index.size();
			Ti n = vct_data.size()-1;
			while (n > 1)
			{
				Ti half = n / 2;
				base = (base[half] < x) ? base+half : base;
				n -= half;
			}

			int off = (*base < x);
			id = base - &vct_index.template get<0>(0) + off;
			Ti v = (base + off != end)?*(base + off):(Ti)-1;
			id = (x == v)?id:vct_data.size()-1;
		}

	public:

		typedef Ti index_type;

		//! Indicate this structure has a function to check the device pointer
		typedef int yes_has_check_device_pointer;

		vector_sparse_gpu_ker(vector_gpu_ker<aggregate<Ti>,layout_base> vct_index,
							  vector_gpu_ker<T,layout_base> vct_data,
							  vector_gpu_ker<aggregate<Ti>,layout_base> vct_add_index,
							  vector_gpu_ker<aggregate<Ti>,layout_base> vct_rem_index,
							  vector_gpu_ker<T,layout_base> vct_add_data,
							  vector_gpu_ker<aggregate<Ti>,layout_base> vct_nadd_index,
							  vector_gpu_ker<aggregate<Ti>,layout_base> vct_nrem_index,
							  int nslot_add,
							  int nslot_rem)
		:vct_index(vct_index),vct_data(vct_data),
		 vct_add_index(vct_add_index),vct_rem_index(vct_rem_index),vct_add_data(vct_add_data),
		 vct_nadd_index(vct_nadd_index),vct_nrem_index(vct_nrem_index),
		 nslot_add(nslot_add),nslot_rem(nslot_rem)
		{}

		/*! \brief Get the number of elements
		 *
		 * \return the number of elements
		 *
		 */
		__device__ inline int size()
		{
			return vct_index.size();
		}

		/*! \brief This function must be called
		 *
		 */
		__device__ inline void init()
		{
#ifdef __NVCC__
			if (threadIdx.x == 0)
			{
				vct_atomic_add = 0;
				vct_atomic_rem = 0;
			}

			__syncthreads();
#endif
		}

		/*! \brief This function must be called
		 *
		 */
		__device__ inline void init_ins_inc()
		{
#ifdef __NVCC__
			if (threadIdx.x == 0)
			{
			    int blockId = dim3CoordToInt(blockIdx, gridDim);
				vct_atomic_add = vct_nadd_index.template get<0>(blockId);
			}

			__syncthreads();
#endif
		}

		/*! \brief This function must be called
		 *
		 */
		__device__ inline void init_rem_inc()
		{
#ifdef __NVCC__
			if (threadIdx.x == 0)
			{
			    int blockId = dim3CoordToInt(blockIdx, gridDim);
			    vct_atomic_rem = vct_nrem_index.template get<0>(blockId);
			}

			__syncthreads();
#endif
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
		__device__ inline openfpm::sparse_index<Ti> get_sparse(Ti id) const
		{
			Ti di;
			this->_branchfree_search(id,di);
			openfpm::sparse_index<Ti> sid;
			sid.id = di;

			return sid;
		}

        /*! \brief Get the background value
         */
        template <unsigned int p>
        __device__ inline auto getBackground() const -> decltype(vct_data.template get<p>(0)) &
        {
            return vct_data.template get<p>(vct_data.size()-1);
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
		__device__ inline auto get(Ti id) const -> decltype(vct_data.template get<p>(id))
		{
			Ti di;
			this->_branchfree_search(id,di);
			return vct_data.template get<p>(di);
		}

        __device__ inline auto get(Ti id) const -> decltype(vct_data.get(0))
        {
            Ti di;
            Ti v = this->_branchfree_search(id,di);
            return vct_data.get(static_cast<size_t>(di));
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
		__device__ inline auto get(openfpm::sparse_index<Ti> id) const -> decltype(vct_data.template get<p>(id.id))
		{
			return vct_data.template get<p>(id.id);
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
		__device__ inline auto get(openfpm::sparse_index<Ti> id) -> decltype(vct_data.template get<p>(id.id))
		{
			return vct_data.template get<p>(id.id);
		}

		/*! \brief Get the index associated to the element id
		 *
		 *
		 * \return the element value requested
		 *
		 */
		__device__ inline Ti get_index(openfpm::sparse_index<Ti> id) const
		{
			return vct_index.template get<0>(id.id);
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
		__device__ inline auto get(Ti id, Ti & di) const -> decltype(vct_data.template get<p>(id))
		{
			this->_branchfree_search(id,di);
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
		template <unsigned int p>
		__device__ inline auto get_ele(Ti di) const -> decltype(vct_data.template get<p>(di))
		{
			return vct_data.template get<p>(di);
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		template <unsigned int p>
		__device__ auto insert(Ti ele) -> decltype(vct_data.template get<p>(0))
		{
#ifdef __NVCC__

            int blockId = dim3CoordToInt(blockIdx, gridDim);
		    int slot_base = blockId;

			int pos = atomicAdd(&vct_atomic_add,1);
			vct_add_index.template get<0>(slot_base*nslot_add+pos) = ele;
			return vct_add_data.template get<p>(slot_base*nslot_add+pos);
#else

			printf("vector_sparse_gpu_ker.insert[1]: Error, this function in order to work is supposed to be compiled with nvcc\n");

#endif
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 * \param ele element to insert
		 *
		 * \return an object to fill the values
		 *
		 */
		__device__ void remove(Ti ele)
		{
#ifdef __NVCC__

            int blockId = dim3CoordToInt(blockIdx, gridDim);
		    int slot_base = blockId;

			int pos = atomicAdd(&vct_atomic_rem,1);
			vct_rem_index.template get<0>(slot_base*nslot_rem+pos) = ele;

#else
			printf("vector_sparse_gpu_ker.remove: Error, this function in order to work is supposed to be compiled with nvcc\n");
#endif
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 * \param ele element to insert
		 *
		 * \return an object to fill the values
		 *
		 */
		__device__ auto insert(Ti ele) -> decltype(vct_add_data.get(0))
		{
#ifdef __NVCC__

            int blockId = dim3CoordToInt(blockIdx, gridDim);
		    int slot_base = blockId;

			int pos = atomicAdd(&vct_atomic_add,1);
			vct_add_index.template get<0>(slot_base*nslot_add+pos) = ele;

			return vct_add_data.get(slot_base*nslot_add+pos);
#else
			printf("vector_sparse_gpu_ker.insert[2]: Error, this function in order to work is supposed to be compiled with nvcc\n");
#endif
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		__device__ void remove_b(Ti ele,Ti slot_base)
		{
#ifdef __NVCC__

			int pos = atomicAdd(&vct_atomic_rem,1);
			vct_rem_index.template get<0>(slot_base*nslot_rem+pos) = ele;

#else
			printf("vector_sparse_gpu_ker.remove_b: Error, this function in order to work is supposed to be compiled with nvcc\n");
#endif
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		template <unsigned int p>
		__device__ auto insert_b(Ti ele,Ti slot_base) -> decltype(vct_data.template get<p>(0))
		{
#ifdef __NVCC__

			int pos = atomicAdd(&vct_atomic_add,1);
			vct_add_index.template get<0>(slot_base*nslot_add+pos) = ele;
			return vct_add_data.template get<p>(slot_base*nslot_add+pos);
#else
			printf("vector_sparse_gpu_ker.insert_b: Error, this function in order to work is supposed to be compiled with nvcc\n");
#endif
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		__device__ auto insert_b(Ti ele,Ti slot_base) -> decltype(vct_add_data.get(0))
		{
#ifdef __NVCC__

			int pos = atomicAdd(&vct_atomic_add,1);
			vct_add_index.template get<0>(slot_base*nslot_add+pos) = ele;
			return vct_add_data.get(slot_base*nslot_add+pos);
#else
			printf("vector_sparse_gpu_ker.insert_b: Error, this function in order to work is supposed to be compiled with nvcc\n");
#endif
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		__device__ void flush_block_insert()
		{
#ifdef __NVCC__

			__syncthreads();

			if (threadIdx.x == 0)
			{
			    int blockId = dim3CoordToInt(blockIdx, gridDim);
			    vct_nadd_index.template get<0>(blockId) = vct_atomic_add;
			}

#else
			printf("vector_sparse_gpu_ker.flush_block_insert: Error, this function in order to work is supposed to be compiled with nvcc\n");
#endif
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		__device__ void flush_block_remove()
		{
#ifdef __NVCC__

			__syncthreads();

			if (threadIdx.x == 0)
			{
			    int blockId = dim3CoordToInt(blockIdx, gridDim);
			    vct_nrem_index.template get<0>(blockId) = vct_atomic_rem;
			}

#else
			printf("vector_sparse_gpu_ker.flush_block_remove: Error, this function in order to work is supposed to be compiled with nvcc\n");
#endif
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		__device__ void flush_block_insert(Ti b, bool flusher)
		{
#ifdef __NVCC__

			__syncthreads();

			if (flusher == true)
			{vct_nadd_index.template get<0>(b) = vct_atomic_add;}


#else
			printf("vector_sparse_gpu_ker.flush_block_insert: Error, this function in order to work is supposed to be compiled with nvcc\n");
#endif
		}

		__device__ auto private_get_data() -> decltype(vct_add_data.getBase().get_data_())
		{
			return vct_add_data.getBase().get_data_();
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		__device__ void flush_block_remove(unsigned int b, bool flusher)
		{
#ifdef __NVCC__

			__syncthreads();

			if (flusher == true)
			{vct_nrem_index.template get<0>(b) = vct_atomic_rem;}

#else
			printf("vector_sparse_gpu_ker.flush_block_remove: Error, this function in order to work is supposed to be compiled with nvcc\n");
#endif
		}

        /*! \brief Get the data buffer
         *
         * \return the reference to the data buffer
         */
        __device__ auto getAddDataBuffer() -> decltype(vct_add_data)&
        {
            return vct_add_data;
        }

        /*! \brief Get the data buffer
         *
         * \return the reference to the data buffer
         */
        __device__ auto getDataBuffer() -> decltype(vct_data)&
        {
            return vct_data;
        }

        /*! \brief Get the indices buffer
        *
        * \return the reference to the indices buffer
        */
        __device__ auto getAddIndexBuffer() const -> const decltype(vct_add_index)&
        {
            return vct_add_index;
        }

        /*! \brief Get the indices buffer
        *
        * \return the reference to the indices buffer
        */
        __device__ auto getIndexBuffer() const -> const decltype(vct_index)&
        {
            return vct_index;
        }

        /*! \brief Get the data buffer
         *
         * \return the reference to the data buffer
         */
        __device__ auto getDataBuffer() const -> const decltype(vct_data)&
        {
            return vct_data;
        }

#ifdef SE_CLASS1

		/*! \brief Check if the device pointer is owned by this structure
		 *
		 * \return a structure pointer check with information about the match
		 *
		 */
		pointer_check check_device_pointer(void * ptr)
		{
			pointer_check pc;

			pc = vct_index.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("Index vector overflow: ") + "\n" + pc.match_str;
				return pc;
			}

			pc = vct_data.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("Data vector overflow: ") + "\n" + pc.match_str;
				return pc;
			}

			pc = vct_add_index.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("Add index vector overflow: ") + "\n" + pc.match_str;
				return pc;
			}

			pc = vct_rem_index.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("Remove index vector overflow: ") + "\n" + pc.match_str;
				return pc;
			}

			pc = vct_nadd_index.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("Add index counter vector overflow: ") + "\n" + pc.match_str;
				return pc;
			}

			pc = vct_nrem_index.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("Remove index counter vector overflow: ") + "\n" + pc.match_str;
				return pc;
			}

			pc = vct_add_data.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("Add data vector overflow: ") + "\n" + pc.match_str;
				return pc;
			}

			return pc;
		}

#endif

	};

	template<typename T,
			 typename Ti,
			 template<typename> class layout_base>
	class vector_sparse_gpu_ker_reduced
	{
			vector_gpu_ker<aggregate<Ti>,layout_base> vct_index;

			vector_gpu_ker<T,layout_base> vct_data;

		public:

			vector_sparse_gpu_ker_reduced(const vector_gpu_ker<aggregate<Ti>,layout_base> & vct_index,
										  const vector_gpu_ker<T,layout_base> & vct_data)
			:vct_index(vct_index),vct_data(vct_data)
			{}

			/*! \brief Return the number of elements
			*
			* \return number of elements
			* 
			*/
			__device__ Ti size()
			{
				return vct_index.size();
			}

			/*! \brief Get the data at element id
			*
			*
			* \return the data element
			*
			*/
			template<unsigned int p>
			__device__ inline auto getData(Ti id) const -> decltype(vct_data.template get<p>(id))
			{
				return vct_data.template get<p>(id);
			}

			/*! \brief Get the data at element id
			*
			*
			* \return the data element
			*
			*/
			template<unsigned int p>
			__device__ inline auto getData(Ti id) -> decltype(vct_data.template get<p>(id))
			{
				return vct_data.template get<p>(id);
			}


			/*! \brief Get the index at element id
			*
			*
			* \return the index element
			*
			*/
			__device__ inline auto getIndex(Ti id) const -> decltype(vct_index.template get<0>(id))
			{
				return vct_index.template get<0>(id);
			}

			/*! \brief Get the index at element id
			*
			*
			* \return the index element
			*
			*/
			__device__ inline auto getIndex(Ti id) -> decltype(vct_index.template get<0>(id))
			{
				return vct_index.template get<0>(id);
			}
	};
}

#endif /* MAP_VECTOR_SPARSE_CUDA_KER_CUH_ */
