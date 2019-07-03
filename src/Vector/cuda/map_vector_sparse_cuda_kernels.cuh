/*
 * map_vector_sparse_cuda_kernels.cuh
 *
 *  Created on: Jan 26, 2019
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_SPARSE_CUDA_KERNELS_CUH_
#define MAP_VECTOR_SPARSE_CUDA_KERNELS_CUH_

#ifdef __NVCC__

#include "util/cuda/cub/util_type.cuh"
#include "util/cuda/cub/block/block_scan.cuh"
#include "util/cuda/moderngpu/operators.hxx"

#endif

template<unsigned int prp>
struct sadd_
{
	typedef boost::mpl::int_<prp> prop;

#ifdef __NVCC__
	template<typename red_type> using op_red = mgpu::plus_t<red_type>;
#endif

	template<typename red_type> __device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		return r1 + r2;
	}

	static bool is_special()
	{
		return false;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{}
};

#ifdef __NVCC__

template<typename type_t, unsigned int blockLength>
struct plus_block_t  : public std::binary_function<type_t, type_t, type_t> {
  MGPU_HOST_DEVICE type_t operator()(type_t a, type_t b) const {
  	type_t res;
  	for (int i=0; i<blockLength; ++i)
  	{
  		res[i] = a[i] + b[i];
  	}
    return res;
  }
};

#endif

template<unsigned int prp, unsigned int blockLength>
struct sadd_block_
{
	typedef boost::mpl::int_<prp> prop;

#ifdef __NVCC__
	template<typename red_type> using op_red = plus_block_t<red_type, blockLength>;
#endif

	template<typename red_type> __device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		red_type res;
		for (int i=0; i<blockLength; ++i)
		{
			res[i] = r1[i] + r2[i];
		}
		return res;
	}

	static bool is_special()
	{
		return false;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{}
};

template<unsigned int prp>
struct smax_
{
	typedef boost::mpl::int_<prp> prop;

#ifdef __NVCC__
	template<typename red_type> using op_red = mgpu::maximum_t<red_type>;
#endif

	template<typename red_type>
	__device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		return (r1 < r2)?r2:r1;
	}

	static bool is_special()
	{
		return false;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{}
};

#ifdef __NVCC__

template<typename type_t, unsigned int blockLength>
struct maximum_block_t  : public std::binary_function<type_t, type_t, type_t> {
  MGPU_HOST_DEVICE type_t operator()(type_t a, type_t b) const {
  	type_t res;
  	for (int i=0; i<blockLength; ++i)
  	{
  		res[i] = max(a[i], b[i]);
  	}
    return res;
  }
};

#endif

template<unsigned int prp, unsigned int blockLength>
struct smax_block_
{
	typedef boost::mpl::int_<prp> prop;

#ifdef __NVCC__
	template<typename red_type> using op_red = maximum_block_t<red_type, blockLength>;
#endif

	template<typename red_type>
	__device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		red_type res;
		for (int i=0; i<blockLength; ++i)
		{
			res[i] = (r1[i] < r2[i])?r2[i]:r1[i];
		}
		return res;
	}

	static bool is_special()
	{
		return false;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{}
};



template<unsigned int prp>
struct smin_
{
	typedef boost::mpl::int_<prp> prop;

#ifdef __NVCC__
	template<typename red_type> using op_red = mgpu::minimum_t<red_type>;
#endif

	template<typename red_type> __device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		return (r1 < r2)?r1:r2;
	}

	static bool is_special()
	{
		return false;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{}
};

#ifdef __NVCC__

template<typename type_t, unsigned int blockLength>
struct minimum_block_t  : public std::binary_function<type_t, type_t, type_t> {
  MGPU_HOST_DEVICE type_t operator()(type_t a, type_t b) const {
  	type_t res;
  	for (int i=0; i<blockLength; ++i)
  	{
  		res[i] = min(a[i], b[i]);
  	}
    return res;
  }
};

#endif

template<unsigned int prp, unsigned int blockLength>
struct smin_block_
{
	typedef boost::mpl::int_<prp> prop;

#ifdef __NVCC__
	template<typename red_type> using op_red = minimum_block_t<red_type, blockLength>;
#endif

	template<typename red_type>
	__device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		red_type res;
		for (int i=0; i<blockLength; ++i)
		{
			res[i] = (r1[i] < r2[i])?r1[i]:r2[i];
		}
		return res;
	}

	static bool is_special()
	{
		return false;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{}
};


#ifdef __NVCC__

template<typename type_t>
struct bitwiseOr_t  : public std::binary_function<type_t, type_t, type_t> {
  MGPU_HOST_DEVICE type_t operator()(type_t a, type_t b) const {
    return a|b;
  }
};

template<unsigned int prp>
struct sBitwiseOr_
{
	typedef boost::mpl::int_<prp> prop;

	template<typename red_type> using op_red = bitwiseOr_t<red_type>;

	template<typename red_type>
	__device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		return r1|r2;
	}

	static bool is_special()
	{
		return false;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{}
};


template<unsigned int prp>
struct sstart_
{
	typedef boost::mpl::int_<prp> prop;

	template<typename red_type> using op_red = mgpu::minimum_t<red_type>;

	template<typename red_type> __device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		return 0;
	}

	static bool is_special()
	{
		return true;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{
		output.template get<0>(i) = seg_prev;
	}
};

template<unsigned int prp>
struct sstop_
{
	typedef boost::mpl::int_<prp> prop;

	template<typename red_type> using op_red = mgpu::minimum_t<red_type>;

	template<typename red_type> __device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		return 0;
	}

	static bool is_special()
	{
		return true;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{
		output.template get<0>(i) = seg_next;
	}
};

template<unsigned int prp>
struct snum_
{
	typedef boost::mpl::int_<prp> prop;

	template<typename red_type> using op_red = mgpu::minimum_t<red_type>;

	template<typename red_type> __device__ __host__ static red_type red(red_type & r1, red_type & r2)
	{
		return 0;
	}

	static bool is_special()
	{
		return true;
	}

	//! is not special reduction so it does not need it
	template<typename seg_type, typename output_type>
	__device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
	{
		output.template get<0>(i) = seg_next - seg_prev;
	}
};


template<typename vector_index_type>
__global__ void construct_insert_list_key_only(vector_index_type vit_block_data,
								 vector_index_type vit_block_n,
								 vector_index_type vit_block_scan,
								 vector_index_type vit_list_0,
								 vector_index_type vit_list_1,
								 int nslot)
{
	int n_move = vit_block_n.template get<0>(blockIdx.x);
	int n_block_move = vit_block_n.template get<0>(blockIdx.x) / blockDim.x;
	int start = vit_block_scan.template get<0>(blockIdx.x);

	int i = 0;
	for ( ; i < n_block_move ; i++)
	{
		vit_list_0.template get<0>(start + i*blockDim.x + threadIdx.x) = vit_block_data.template get<0>(nslot*blockIdx.x + i*blockDim.x + threadIdx.x);
		vit_list_1.template get<0>(start + i*blockDim.x + threadIdx.x) = start + i*blockDim.x + threadIdx.x;
	}

	// move remaining
	if (threadIdx.x < n_move - i*blockDim.x )
	{
		vit_list_0.template get<0>(start + i*blockDim.x + threadIdx.x) = vit_block_data.template get<0>(nslot*blockIdx.x + i*blockDim.x + threadIdx.x);
		vit_list_1.template get<0>(start + i*blockDim.x + threadIdx.x) = start + i*blockDim.x + threadIdx.x;
	}
}

template<typename vector_index_type, typename vector_data_type>
__global__ void construct_insert_list(vector_index_type vit_block_data,
								 vector_index_type vit_block_n,
								 vector_index_type vit_block_scan,
								 vector_index_type vit_list_0,
								 vector_index_type vit_list_1,
								 vector_data_type vdata_block,
								 vector_data_type vdata,
								 int nslot)
{
	int n_move = vit_block_n.template get<0>(blockIdx.x);
	int n_block_move = vit_block_n.template get<0>(blockIdx.x) / blockDim.x;
	int start = vit_block_scan.template get<0>(blockIdx.x);

	int i = 0;
	for ( ; i < n_block_move ; i++)
	{
		vit_list_0.template get<0>(start + i*blockDim.x + threadIdx.x) = vit_block_data.template get<0>(nslot*blockIdx.x + i*blockDim.x + threadIdx.x);
		vit_list_1.template get<0>(start + i*blockDim.x + threadIdx.x) = start + i*blockDim.x + threadIdx.x;
		vdata.get(start + i*blockDim.x + threadIdx.x) = vdata_block.get(nslot*blockIdx.x + i*blockDim.x + threadIdx.x);
	}

	// move remaining
	if (threadIdx.x < n_move - i*blockDim.x )
	{
		vit_list_0.template get<0>(start + i*blockDim.x + threadIdx.x) = vit_block_data.template get<0>(nslot*blockIdx.x + i*blockDim.x + threadIdx.x);
		vit_list_1.template get<0>(start + i*blockDim.x + threadIdx.x) = start + i*blockDim.x + threadIdx.x;
		vdata.get(start + i*blockDim.x + threadIdx.x) = vdata_block.get(nslot*blockIdx.x + i*blockDim.x + threadIdx.x);;
	}
}

template<typename vector_index_type>
__global__ void construct_remove_list(vector_index_type vit_block_data,
								 vector_index_type vit_block_n,
								 vector_index_type vit_block_scan,
								 vector_index_type vit_list_0,
								 vector_index_type vit_list_1,
								 int nslot)
{
	int n_move = vit_block_n.template get<0>(blockIdx.x);
	int n_block_move = vit_block_n.template get<0>(blockIdx.x) / blockDim.x;
	int start = vit_block_scan.template get<0>(blockIdx.x);

	int i = 0;
	for ( ; i < n_block_move ; i++)
	{
		vit_list_0.template get<0>(start + i*blockDim.x + threadIdx.x) = vit_block_data.template get<0>(nslot*blockIdx.x + i*blockDim.x + threadIdx.x);
		vit_list_1.template get<0>(start + i*blockDim.x + threadIdx.x) = start + i*blockDim.x + threadIdx.x;
	}

	// move remaining
	if (threadIdx.x < n_move - i*blockDim.x )
	{
		vit_list_0.template get<0>(start + i*blockDim.x + threadIdx.x) = vit_block_data.template get<0>(nslot*blockIdx.x + i*blockDim.x + threadIdx.x);
		vit_list_1.template get<0>(start + i*blockDim.x + threadIdx.x) = start + i*blockDim.x + threadIdx.x;
	}
}

template<typename e_type, typename v_reduce>
struct data_merger
{
	e_type src1;
	e_type src2;

	e_type dst;

	__device__ __host__ inline data_merger(const e_type & src1, const e_type & src2, const e_type & dst)
	:src1(src1),src2(src2),dst(dst)
	{
	};


	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<v_reduce,T>::type red_type;

		dst.template get<red_type::prop::value>() = red_type::red(src1.template get<red_type::prop::value>(),src2.template get<red_type::prop::value>());
	}
};

template<typename vector_index_type, typename vector_data_type, typename vector_index_type2, unsigned int block_dim, typename ... v_reduce>
__global__ void solve_conflicts(vector_index_type vct_index, vector_data_type vct_data,
								vector_index_type merge_index, vector_data_type vct_add_data,
		 	 	 	 	 	 	vector_index_type vct_index_out, vector_data_type vct_data_out,
		 	 	 	 	 	 	vector_index_type2 vct_tot_out,
		 	 	 	 	 	 	int base)
{
	typedef typename std::remove_reference<decltype(vct_index.template get<0>(0))>::type index_type;

    // Specialize BlockScan for a 1D block of 256 threads on type int
    typedef cub::BlockScan<int, block_dim> BlockScan;
    // Allocate shared memory for BlockScan
    __shared__ typename BlockScan::TempStorage temp_storage;

	int p = blockIdx.x * blockDim.x + threadIdx.x;

	if (p >= vct_index.size())	return;

	index_type id_check = (p == vct_index.size() - 1)?(index_type)-1:vct_index.template get<0>(p+1);
	int predicate = vct_index.template get<0>(p) != id_check;

	int scan = predicate;

	// in shared memory scan
	BlockScan(temp_storage).ExclusiveSum(scan, scan);

	if (predicate == 1)
	{
		vct_index_out.template get<0>(blockIdx.x*block_dim + scan) = vct_index.template get<0>(p);

		int index1 = merge_index.template get<0>(p);

		auto e = vct_data_out.get(blockIdx.x*block_dim + scan);

		if (index1 < base)
		{
			e = vct_data.get(index1);
			vct_data_out.get(blockIdx.x*block_dim + scan) = e;
		}
		else
		{
			e = vct_add_data.get(index1 - base);
			vct_data_out.get(blockIdx.x*block_dim + scan) = e;
		}
	}

	__syncthreads();

	if (predicate == 0)
	{
		//! at the border of the block the index must be copied
		if (threadIdx.x == blockDim.x-1)
		{vct_index_out.template get<0>(blockIdx.x*block_dim + scan) = vct_index.template get<0>(p);}

		// we have to merge the data

		typedef boost::mpl::vector<v_reduce ...> v_reduce_;

		int index1 = merge_index.template get<0>(p);
		int index2 = merge_index.template get<0>(p+1) - base;

		data_merger<decltype(vct_data.get(p)),v_reduce_> dm(vct_data.get(index1),
															vct_add_data.get(index2),
															vct_data_out.get(blockIdx.x*block_dim + scan));

		// data_merge
		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(v_reduce)>>(dm);
	}

	if (threadIdx.x == blockDim.x - 1 || p == vct_index.size() - 1)
	{
		vct_tot_out.template get<0>(blockIdx.x) = scan + predicate;
		vct_tot_out.template get<2>(blockIdx.x) = predicate;
	}
}

template<typename vector_index_type, typename vector_index_type2, unsigned int block_dim>
__global__ void solve_conflicts_remove(vector_index_type vct_index,
								vector_index_type merge_index,
		 	 	 	 	 	 	vector_index_type vct_index_out,
		 	 	 	 	 	 	vector_index_type vct_index_out_ps,
		 	 	 	 	 	 	vector_index_type2 vct_tot_out,
		 	 	 	 	 	 	int base)
{
	typedef typename std::remove_reference<decltype(vct_index.template get<0>(0))>::type index_type;

    // Specialize BlockScan for a 1D block of 256 threads on type int
    typedef cub::BlockScan<int, block_dim> BlockScan;
    // Allocate shared memory for BlockScan
    __shared__ typename BlockScan::TempStorage temp_storage;

	int p = blockIdx.x * blockDim.x + threadIdx.x;

	if (p >= vct_index.size())	return;

	index_type id_check_n = (p == vct_index.size() - 1)?(index_type)-1:vct_index.template get<0>(p+1);
	index_type id_check_p = (p == 0)?(index_type)-1:vct_index.template get<0>(p-1);
	index_type id_check = vct_index.template get<0>(p);
	int predicate = id_check != id_check_p;
	predicate &= id_check != id_check_n;
	int mi = merge_index.template get<0>(p);
	predicate &= (mi < base);

	int scan = predicate;

	// in shared memory scan
	BlockScan(temp_storage).ExclusiveSum(scan, scan);

	if (predicate == 1)
	{
		vct_index_out.template get<0>(blockIdx.x*block_dim + scan) = vct_index.template get<0>(p);
		vct_index_out_ps.template get<0>(blockIdx.x*block_dim + scan) = merge_index.template get<0>(p);
	}

	__syncthreads();

	if (threadIdx.x == blockDim.x - 1 || p == vct_index.size() - 1)
	{
		vct_tot_out.template get<0>(blockIdx.x) = scan + predicate;
	}
}

template<typename vector_type, typename vector_type2, typename red_op>
__global__ void reduce_from_offset(vector_type segment_offset,
		                           vector_type2 output,
		                           typename std::remove_reference<decltype(segment_offset.template get<1>(0))>::type max_index)
{
	int p = blockIdx.x * blockDim.x + threadIdx.x;

	if (p >= segment_offset.size())	return;

	typename std::remove_reference<decltype(segment_offset.template get<1>(0))>::type v;
	if (p == segment_offset.size()-1)
	{v = max_index;}
	else
	{v = segment_offset.template get<1>(p+1);}

	red_op::set(v,segment_offset.template get<1>(p),output,p);
}

template<typename vector_index_type, typename vector_data_type, typename vector_index_type2>
__global__ void realign(vector_index_type vct_index, vector_data_type vct_data,
		 	 	 	 	 	 	vector_index_type vct_index_out, vector_data_type vct_data_out,
		 	 	 	 	 	 	vector_index_type2 vct_tot_out_scan)
{
	int p = blockIdx.x * blockDim.x + threadIdx.x;

	if (p >= vct_index.size())	return;

	int tot = vct_tot_out_scan.template get<0>(blockIdx.x);

	if (threadIdx.x > tot)
	{return;}

	if (threadIdx.x == tot && vct_tot_out_scan.template get<2>(blockIdx.x) == 1)
	{return;}

	// this branch exist if the previous block (last thread) had a predicate == 0 in resolve_conflict in that case
	// the thread 0 of the next block does not have to produce any data
	if (threadIdx.x == 0 && blockIdx.x != 0 && vct_tot_out_scan.template get<2>(blockIdx.x - 1) == 0)
	{return;}

	int ds = vct_tot_out_scan.template get<1>(blockIdx.x);

	vct_index_out.template get<0>(ds+threadIdx.x) = vct_index.template get<0>(p);

	auto src = vct_data.get(p);
	auto dst = vct_data_out.get(ds+threadIdx.x);

	dst = src;
}

template<typename vector_index_type, typename vct_data_type, typename vector_index_type2>
__global__ void realign_remove(vector_index_type vct_index, vector_index_type vct_m_index, vct_data_type vct_data,
		 	 	 	 	 	 	vector_index_type vct_index_out, vct_data_type vct_data_out,
		 	 	 	 	 	 	vector_index_type2 vct_tot_out_scan)
{
	int p = blockIdx.x * blockDim.x + threadIdx.x;

	if (p >= vct_index.size())	return;

	int tot = vct_tot_out_scan.template get<0>(blockIdx.x);

	if (threadIdx.x >= tot)
	{return;}

	int ds = vct_tot_out_scan.template get<1>(blockIdx.x);

	vct_index_out.template get<0>(ds+threadIdx.x) = vct_index.template get<0>(p);

	int oi = vct_m_index.template get<0>(p);

	auto src = vct_data.get(oi);
	auto dst = vct_data_out.get(ds+threadIdx.x);

	dst = src;
}

template<typename vector_index_type, typename vector_data_type>
__global__ void reorder_vector_data(vector_index_type vi, vector_data_type v_data, vector_data_type v_data_ord)
{
	int p = blockIdx.x * blockDim.x + threadIdx.x;

	if (p >= vi.size())	return;

	// reorder based on v_ord the data vector

	v_data_ord.get_o(p) = v_data.get_o(vi.template get<0>(p));
}

template<unsigned int prp, typename vector_type>
__global__ void set_indexes(vector_type vd, int off)
{
	int p = blockIdx.x * blockDim.x + threadIdx.x;

	vd.template get<prp>(p) = p + off;
}

#endif

#endif /* MAP_VECTOR_SPARSE_CUDA_KERNELS_CUH_ */
