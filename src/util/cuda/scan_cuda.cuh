/*
 * scan_cuda.cuh
 *
 *  Created on: Jun 28, 2018
 *      Author: i-bird
 */

#ifndef SCAN_CUDA_CUH_
#define SCAN_CUDA_CUH_

#include "Vector/map_vector.hpp"

#ifdef __NVCC__

////////////////// Ratio reduction ///////////////////////////////////////////////////////////////

/*template<typename cnt_type, typename ids_type, unsigned int rt = sizeof(cnt_type)/sizeof(ids_type)>
struct ratio_reduction
{
	static inline __device__ cnt_type reduction(cnt_type * rd)
	{
		return -1;
	}
};

template<typename cnt_type, typename ids_type>
struct ratio_reduction<cnt_type,ids_type,1>
{
	static inline __device__ unsigned int reduction(cnt_type * val)
	{
		return val[0] + val[1] + val[2] + val[3];
	}
};

template<typename cnt_type, typename ids_type>
struct ratio_reduction<cnt_type,ids_type,2>
{
	static inline __device__ unsigned int reduction(cnt_type * val)
	{
	    val[0] = ((val[0] & 0xFFFF0000) >> 16) + (val[0] & 0xFFFF);
	    val[1] = ((val[1] & 0xFFFF0000) >> 16) + (val[1] & 0xFFFF);
	    val[2] = ((val[2] & 0xFFFF0000) >> 16) + (val[2] & 0xFFFF);
	    val[3] = ((val[3] & 0xFFFF0000) >> 16) + (val[3] & 0xFFFF);
	    return val[0] + val[1] + val[2] + val[3];
	}
};

template<typename cnt_type, typename ids_type>
struct ratio_reduction<cnt_type,ids_type,4>
{
	static inline __device__ unsigned int reduction(cnt_type * val)
	{
	    val[0] = ((val[0] & 0xFF000000) >> 24) + ((val[0] & 0xFF0000) >> 16) + ((val[0] & 0xFF00) >> 8) + (val[0] & 0xFF);
	    val[1] = ((val[1] & 0xFF000000) >> 24) + ((val[1] & 0xFF0000) >> 16) + ((val[1] & 0xFF00) >> 8) + (val[1] & 0xFF);
	    val[2] = ((val[2] & 0xFF000000) >> 24) + ((val[2] & 0xFF0000) >> 16) + ((val[2] & 0xFF00) >> 8) + (val[2] & 0xFF);
	    val[3] = ((val[3] & 0xFF000000) >> 24) + ((val[3] & 0xFF0000) >> 16) + ((val[3] & 0xFF00) >> 8) + (val[3] & 0xFF);
	    return val[0] + val[1] + val[2] + val[3];
	}
};


///////////////////////// Ratio extends ////////////////////////////////////////////////////////////////////


template<typename T>
struct to_four_type
{
	typedef T type;
};

template<>
struct to_four_type<unsigned int>
{
	typedef uint4 type;
};

template<>
struct to_four_type<int>
{
	typedef int4 type;
};

template<typename T>
struct make4
{
	__device__ inline static T make()
	{
		return 0;
	}
};

template<>
struct make4<uint4>
{
	__device__ inline static uint4 make()
	{
		return make_uint4(0,0,0,0);
	}
};

template<>
struct make4<int4>
{
	__device__ inline static int4 make()
	{
		return make_int4(0,0,0,0);
	}
};

template<typename cnt_type, typename ids_type, unsigned int rt = sizeof(cnt_type)/sizeof(ids_type)>
struct ratio_extend
{
	static inline __device__ cnt_type reduction(ids_type * rd)
	{
		return -1;
	}
};

template<typename cnt_type, typename ids_type>
struct ratio_extend<cnt_type,ids_type,1>
{
	typedef boost::mpl::int_<1> size;

	typedef cnt_type cnt_type_;
	typedef typename to_four_type<cnt_type>::type cnt_type4;

	static inline __device__ void extend(cnt_type4 * val)
	{
		unsigned int tmp = val[0].w;
	    val[0].w = val[0].z;
	    val[0].z = val[0].y;
	    val[0].y = val[0].x;
	    val[0].x = tmp;
	}

	static inline __device__ unsigned int reduce(cnt_type4 * tu4, cnt_type4 * val)
	{
		tu4->x = val[0].x + val[0].y + val[0].z + val[0].w;

	    return tu4->x;
	}

	static inline __device__ void scan(unsigned int woff_w, cnt_type4 tu4, cnt_type4 * val)
	{
		unsigned int lps;
	    lps = woff_w + tu4.w;

	    val[0].x = lps;

	    val[0].y += val[0].x;

	    val[0].z += val[0].y;

	    val[0].w += val[0].z;
	}

	static inline __device__ cnt_type4 zero()
	{
		return make4<cnt_type4>::make();
	}
};



template<typename cnt_type, typename ids_type>
struct ratio_extend<cnt_type,ids_type,2>
{
	typedef boost::mpl::int_<2> size;

	typedef cnt_type cnt_type_;
	typedef typename to_four_type<cnt_type>::type cnt_type4;

	static inline __device__ void extend(cnt_type4 * val)
	{
	    val[1].w = (val[0].w & 0x0000FFFF);
	    val[1].z = (val[0].z & 0xFFFF0000) >>  16;

	    val[1].y = (val[0].z & 0x0000FFFF);
	    val[1].x = (val[0].w & 0xFFFF0000) >>  16;

	    val[0].w = (val[0].y & 0x0000FFFF);
	    val[0].z = (val[0].x & 0xFFFF0000) >>  16;

	    short int tmp = (val[0].y & 0xFFFF0000) >>  16;
	    val[0].y = (val[0].x & 0x0000FFFF);
	    val[0].x = tmp;
	}

	static inline __device__ unsigned int reduce(cnt_type4 * tu4, cnt_type4 * val)
	{
	    tu4->x = val[0].x + val[0].y + val[0].z + val[0].w;
	    tu4->y = val[1].x + val[1].y + val[1].z + val[1].w;

	    return tu4->x + tu4->y;
	}

	static inline __device__ void scan(unsigned int woff_w, cnt_type4 tu4, cnt_type4 * val)
	{
		uint2 lps;
	    lps.x = woff_w + tu4.w;
	    lps.y = lps.x + tu4.x;

	    val[0].x = lps.x;
	    val[1].x = lps.y;

	    val[0].y += val[0].x;
	    val[1].y += val[1].x;

	    val[0].z += val[0].y;
	    val[1].z += val[1].y;

	    val[0].w += val[0].z;
	    val[1].w += val[1].z;
	}

	static inline __device__ cnt_type4 zero()
	{
		return make4<cnt_type4>::make();
	}
};

template<typename cnt_type, typename ids_type>
struct ratio_extend<cnt_type,ids_type,4>
{
	typedef boost::mpl::int_<4> size;

	typedef cnt_type cnt_type_;
	typedef typename to_four_type<cnt_type>::type cnt_type4;

	static inline __device__ void extend(cnt_type4 * val)
	{
	    val[3].y = (val[0].w & 0x000000FF);
	    val[3].z = (val[0].w & 0x0000FF00) >>  8;
	    val[3].w = (val[0].w & 0x00FF0000) >> 16;
	    val[3].x = (val[0].w & 0xFF000000) >> 24;

	    val[2].y = (val[0].z & 0x000000FF);
	    val[2].z = (val[0].z & 0x0000FF00) >>  8;
	    val[2].w = (val[0].z & 0x00FF0000) >> 16;
	    val[2].x = (val[0].z & 0xFF000000) >> 24;

	    val[1].y = (val[0].y & 0x000000FF);
	    val[1].z = (val[0].y & 0x0000FF00) >>  8;
	    val[1].w = (val[0].y & 0x00FF0000) >> 16;
	    val[1].x = (val[0].y & 0xFF000000) >> 24;

	    val[0].y = (val[0].x & 0x000000FF);
	    val[0].z = (val[0].x & 0x0000FF00) >>  8;
	    val[0].w = (val[0].x & 0x00FF0000) >> 16;
	    val[0].x = (val[0].x & 0xFF000000) >> 24;
	}

	static inline __device__ unsigned int reduce(cnt_type4 * tu4, cnt_type4 * val)
	{
	    tu4->x = val[0].x + val[0].y + val[0].z + val[0].w;
	    tu4->y = val[1].x + val[1].y + val[1].z + val[1].w;
	    tu4->z = val[2].x + val[2].y + val[2].z + val[2].w;
	    tu4->w = val[3].x + val[3].y + val[3].z + val[3].w;

	    return tu4->x + tu4->y + tu4->z + tu4->w;
	}

	static inline __device__ void scan(unsigned int woff_w, cnt_type4 tu4, cnt_type4 * val)
	{
		uint4 lps;
	    lps.x = woff_w + tu4.w;
	    lps.y = lps.x + tu4.x;
	    lps.z = lps.y + tu4.y;
	    lps.w = lps.z + tu4.z;

	    val[0].x = lps.x;
	    val[1].x = lps.y;
	    val[2].x = lps.z;
	    val[3].x = lps.w;

	    val[0].y += val[0].x;
	    val[1].y += val[1].x;
	    val[2].y += val[2].x;
	    val[3].y += val[3].x;

	    val[0].z += val[0].y;
	    val[1].z += val[1].y;
	    val[2].z += val[2].y;
	    val[3].z += val[3].y;

	    val[0].w += val[0].z;
	    val[1].w += val[1].z;
	    val[2].w += val[2].z;
	    val[3].w += val[3].z;
	}

	static inline __device__ cnt_type4 zero()
	{
		return make4<cnt_type4>::make();
	}
};*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////

/*template<typename cnt_type, typename ids_type>
__global__ void compress4(const cnt_type n_cell, const cnt_type *const count, ids_type *const compressed)
{
    const int gid = 4*(threadIdx.x + blockDim.x * blockIdx.x);

    if (gid >= n_cell)
    {return;}

    compressed[gid] = count[gid];
    compressed[gid+1] = count[gid+1];
    compressed[gid+2] = count[gid+2];
    compressed[gid+3] = count[gid+3];
}

template <int NWARP, typename cnt_type, typename ids_type, typename reduction>
__global__ void breduce(int n, const cnt_type *vin, cnt_type *vout)
{
    const int wid = threadIdx.x/32;
    const int lid = threadIdx.x%32;
    const int tid = 4*(blockDim.x*blockIdx.x + threadIdx.x);
    cnt_type val[4] = {0,0,0,0};

    __shared__ unsigned int shtmp[NWARP];

    if (tid < n)
    {
    	val[0] = vin[0+tid];
    	val[1] = vin[1+tid];
    	val[2] = vin[2+tid];
    	val[3] = vin[3+tid];
    }

    val[0] = reduction::reduction(val);

#pragma unroll
    for(int i = 16; i > 0; i >>= 1)
    {val[0] += __shfl_down_sync(0xFFFFFFFF,(int)val[0],i);}

    if (0 == lid)
    {shtmp[wid] = val[0];}

    __syncthreads();
    if (0 == wid)
    {
        val[0] = (lid < NWARP) ? shtmp[lid] : 0;

#pragma unroll
        for(int i = 16; i > 0; i >>= 1)
        {val[0] += __shfl_down_sync(0xFFFFFFFF,(int)val[0], i);}
    }

    if (0 == threadIdx.x)
    {
    	vout[blockIdx.x] = val[0];
    }
}

template <int BDIM, typename cnt_type>
__global__ void bexscan(int n, cnt_type *v)
{
    extern __shared__ unsigned int shtmp[];

    for(int i = threadIdx.x; i < n; i += BDIM)
    {shtmp[i] = v[i];}

    int i = threadIdx.x;
    __syncthreads();
    for(; n >= BDIM; i += BDIM, n -= BDIM)
    {
        __syncthreads();
        if (i > 0 && 0 == threadIdx.x)
        {shtmp[i] += shtmp[i-1];}

        unsigned int a = 0;

#pragma unroll
        for(int j = 1; j < BDIM; j <<= 1)
        {
            a = 0;

            __syncthreads();
            if (threadIdx.x >= j) {a = shtmp[i] + shtmp[i-j];}

            __syncthreads();
            if (threadIdx.x >= j) {shtmp[i] = a;}
        }
        v[i] = shtmp[i];
    }
    if (threadIdx.x < n)
    {
        __syncthreads();
        if (i > 0 && 0 == threadIdx.x)
        {shtmp[i] += shtmp[i-1];}

        unsigned int a = 0;
        for(int j = 1; j < BDIM; j <<= 1)
        {
            a = 0;

            __syncthreads();
            if (threadIdx.x >= j) a = shtmp[i] + shtmp[i-j];

            __syncthreads();
            if (threadIdx.x >= j) shtmp[i] = a;
        }
        v[i] = shtmp[i];
    }
}

template <int NWARP, typename extend>
__global__ void gexscan(int n,
		                const typename extend::cnt_type4 *vin,
		                typename extend::cnt_type_ *offs,
						typename extend::cnt_type4 *vout)
{
    const int wid = threadIdx.x/32;
    const int lid = threadIdx.x%32;
    const int tid = blockDim.x*blockIdx.x + threadIdx.x;
    typename extend::cnt_type4 val[extend::size::value];

    __shared__ unsigned int woff[NWARP];

    if (tid < n) val[0] = vin[tid];
    else         val[0] = extend::zero();

    extend::extend(val);

    typename extend::cnt_type4 tu4;
    typename extend::cnt_type_ tmp = extend::reduce(&tu4,val);

    tu4.w = tmp;
#pragma unroll
    for(int i = 1; i < 32; i <<= 1)
    {tu4.w += (lid >= i)*__shfl_up_sync(0xFFFFFFFF,(int)tu4.w, i);}

    if (lid == 31)
    {
        if (wid < NWARP-1) woff[wid+1] = tu4.w;
        else               woff[0]     = (blockIdx.x > 0) ? offs[blockIdx.x-1] : 0;
    }

    tu4.w -= tmp;

    __syncthreads();

    if (0 == wid)
    {
        tmp = (lid < NWARP) ? woff[lid] : 0;
#pragma unroll
        for(int i = 1; i < NWARP; i <<= 1)
        {tmp += (lid >= i)*__shfl_up_sync(0xFFFFFFFF,(int)tmp, i);}

        if (lid < NWARP) woff[lid] = tmp;
    }
    __syncthreads();

    if (tid >= n) return;

    extend::scan(woff[wid],tu4,val);

    vout += tid*extend::size::value;

#pragma unroll
    for(int i = 0; i < extend::size::value; i++)
    {vout[i] = val[i];}
}*/

/*! \brief Scan is a class because it use internally temporary buffers that are heavy to reconstruct
 *
 *
 *
 */
/*template<typename cnt_type, typename ids_type>
class scan
{
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> red;
	openfpm::vector<aggregate<ids_type>,CudaMemory,typename memory_traits_inte<aggregate<ids_type>>::type,memory_traits_inte> compressed;

public:

	void scan_(openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> & cl_n,
			  openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> & cl_n_scan)
	{
		constexpr int THREADS = 128;
		constexpr int ratio = 4*sizeof(cnt_type)/sizeof(ids_type);

		int nblocks = (cl_n.size() + THREADS * ratio - 1 ) / (THREADS * ratio);
		red.resize(nblocks);

		auto ite = cl_n.getGPUIterator();

		compressed.resize(cl_n.size());
		compress4<cnt_type,ids_type><<<ite.wthr,ite.thr>>>((cnt_type)cl_n.size(),
														   static_cast<cnt_type *>(cl_n.template getDeviceBuffer<0>()),
														   static_cast<ids_type *>(compressed.template getDeviceBuffer<0>()));

		breduce<THREADS/32,cnt_type,ids_type,ratio_reduction<cnt_type,ids_type>><<<nblocks, THREADS                >>>(cl_n.size() / ratio * 4,
																  (cnt_type *)compressed.template getDeviceBuffer<0>(),
																  static_cast<cnt_type *>(red.template getDeviceBuffer<0>()));


		bexscan<THREADS,cnt_type><<<1, THREADS, nblocks*sizeof(uint)>>>(nblocks,
															   static_cast<cnt_type *>(red.template getDeviceBuffer<0>()));

		size_t raw_size = cl_n.size();

		// resize to be multiple of 16

		size_t ml = ((raw_size + ratio - 1) / ratio) *ratio;
		cl_n_scan.resize(ml);

		gexscan<THREADS/32,ratio_extend<cnt_type,ids_type>> <<< nblocks, THREADS >>>((cl_n.size() + ratio - 1 ) / ratio,
																							  static_cast<typename ratio_extend<cnt_type,ids_type>::cnt_type4 *>(compressed.template getDeviceBuffer<0>()),
																							  static_cast<cnt_type *>(red.template getDeviceBuffer<0>()),
																							  static_cast<typename ratio_extend<cnt_type,ids_type>::cnt_type4 *>(cl_n_scan.template getDeviceBuffer<0>()));

		cl_n_scan.resize(raw_size);
	}

};*/

#else

// In case we do not have NVCC we create a stub

template<typename cnt_type, typename ids_type>
class scan
{
};

#endif

#endif /* SCAN_CUDA_CUH_ */
