/*
 * map_vector_sparse_cuda_ker.cuh
 *
 *  Created on: Jan 23, 2019
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_SPARSE_CUDA_KER_CUH_
#define MAP_VECTOR_SPARSE_CUDA_KER_CUH_

namespace openfpm
{
	template<typename index_type>
	struct sparse_index
	{
		index_type id;

		__device__ sparse_index(index_type id)
		:id(id)
		{}
	};

	static __shared__ int vct_atomic_add;
	static __shared__ int vct_atomic_rem;

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
		mutable T bck;

		int nslot_add;
		int nslot_rem;

		/*! \brief get the element i
		 *
		 * search the element x
		 *
		 * \param i element i
		 */
		inline __device__ Ti _branchfree_search(Ti x, Ti & id) const
		{
			if (vct_index.size() == 0)	{return (Ti)-1;}
			const Ti *base = &vct_index.template get<0>(0);
			Ti n = vct_data.size();
			while (n > 1)
			{
				Ti half = n / 2;
				base = (base[half] < x) ? base+half : base;
				n -= half;
			}

			int off = (*base < x);
			id = base - &vct_index.template get<0>(0) + off;
			return *(base + off);
		}

	public:

		typedef Ti index_type;

		vector_sparse_gpu_ker(vector_gpu_ker<aggregate<Ti>,layout_base> vct_index,
							  vector_gpu_ker<T,layout_base> vct_data,
							  vector_gpu_ker<aggregate<Ti>,layout_base> vct_add_index,
							  vector_gpu_ker<aggregate<Ti>,layout_base> vct_rem_index,
							  vector_gpu_ker<T,layout_base> vct_add_data,
							  vector_gpu_ker<aggregate<Ti>,layout_base> vct_nadd_index,
							  vector_gpu_ker<aggregate<Ti>,layout_base> vct_nrem_index,
							  T & bck,
							  int nslot_add,
							  int nslot_rem)
		:vct_index(vct_index),vct_data(vct_data),
		 vct_add_index(vct_add_index),vct_rem_index(vct_rem_index),vct_add_data(vct_add_data),
		 vct_nadd_index(vct_nadd_index),vct_nrem_index(vct_nrem_index),
		 bck(bck),nslot_add(nslot_add),nslot_rem(nslot_rem)
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
				vct_atomic_add = vct_nadd_index.template get<0>(blockIdx.x);
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
				vct_atomic_rem = vct_nrem_index.template get<0>(blockIdx.x);
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
			Ti v = this->_branchfree_search(id,di);

			openfpm::sparse_index<Ti> sid((v == id)?di:(Ti)-1);

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
		__device__ inline auto get(Ti id) const -> decltype(vct_data.template get<p>(id))
		{
			Ti di;
			Ti v = this->_branchfree_search(id,di);
			return (v == id)?vct_data.template get<p>(di):bck.template get<p>();
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
			Ti v = this->_branchfree_search(id,di);
			di = (v != id)?-1:di;
			return (v == id)?vct_data.template get<p>(di):bck.template get<p>();
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
			return (di != -1)?vct_data.template get<p>(di):bck.template get<p>();
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		template <unsigned int p>
		__device__ auto insert(Ti ele) -> decltype(vct_data.template get<p>(0))
		{
#ifdef __NVCC__

			int slot_base = blockIdx.x;

			int pos = atomicAdd(&vct_atomic_add,1);
			vct_add_index.template get<0>(slot_base*nslot_add+pos) = ele;
			return vct_add_data.template get<p>(slot_base*nslot_add+pos);
#else
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
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

			int slot_base = blockIdx.x;

			int pos = atomicAdd(&vct_atomic_rem,1);
			vct_rem_index.template get<0>(slot_base*nslot_rem+pos) = ele;

#else
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
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

			int slot_base = blockIdx.x;

			int pos = atomicAdd(&vct_atomic_add,1);
			vct_add_index.template get<0>(slot_base*nslot_add+pos) = ele;
			return vct_add_data.get(slot_base*nslot_add+pos);
#else
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
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
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
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
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
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
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
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
			{vct_nadd_index.template get<0>(blockIdx.x) = vct_atomic_add;}

#else
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
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
			{vct_nrem_index.template get<0>(blockIdx.x) = vct_atomic_rem;}

#else
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
#endif
		}

		/*! \brief It insert an element in the sparse vector
		 *
		 *
		 */
		__device__ void flush_block_insert(unsigned int b, bool flusher)
		{
#ifdef __NVCC__

			__syncthreads();

			if (flusher == true)
			{vct_nadd_index.template get<0>(b) = vct_atomic_add;}


#else
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
#endif
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
			std::cout << __FILE__ << ":" << __LINE__ << " Error, this function in order to work is supposed to be compiled with nvcc" << std::endl;
#endif
		}
	};
}

#endif /* MAP_VECTOR_SPARSE_CUDA_KER_CUH_ */
