/*
 * Cuda_cell_list_util_func.hpp
 *
 *  Created on: Jun 17, 2018
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CUDA_CUDA_CELL_LIST_UTIL_FUNC_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CUDA_CUDA_CELL_LIST_UTIL_FUNC_HPP_

#include <boost/integer/integer_mask.hpp>
#include <Vector/map_vector_sparse.hpp>

template<unsigned int dim, typename cnt_type, typename ids_type, typename transform>
struct cid_
{
	static inline __device__ __host__ cnt_type get_cid(openfpm::array<ids_type,1,cnt_type> & div_c , ids_type * e)
	{
		cnt_type id = e[dim-1];

#pragma unroll
		for (int i = 1; i >= 0 ; i-- )
		{id = e[i] + div_c[i]*id;}

		return id;
	}

	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,dim,cnt_type> & div_c , const grid_key_dx<1,cnt_type> & e)
	{
		cnt_type id = e.get(dim-1);

#pragma unroll
		for (int i = 1; i >= 0 ; i-- )
		{id = e.get(i) + div_c[i]*id;}

		return id;
	}

	template<typename T> static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,dim,cnt_type> & div_c,
			                                                       const openfpm::array<T,dim,cnt_type> & spacing,
			                                                       const transform & t,
			                                                       const Point<dim,T> & p)
	{
		cnt_type id = p.get(dim-1) / spacing[dim-1];

#pragma unroll
		for (int i = 1; i >= 0 ; i-- )
		{id = t.transform(p.get(i),i) / spacing[i] + div_c[i]*id;}

		return id;
	}
};

template<typename cnt_type, typename ids_type, typename transform>
struct cid_<1,cnt_type,ids_type, transform>
{
	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,1,cnt_type> & div_c, ids_type * e)
	{
		return e[0];
	}

	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,1,cnt_type> & div_c, const grid_key_dx<1,ids_type> & e)
	{
		return e.get(0);
	}

	template<typename T> static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,1,cnt_type> & div_c,
			                                                       const openfpm::array<T,1,cnt_type> & spacing,
			                                                       const transform & t,
			                                                       const Point<1,T> & p)
	{
		return t.transform(p.get(0),0) / spacing[0];
	}
};

template<typename cnt_type, typename ids_type, typename transform>
struct cid_<2,cnt_type,ids_type,transform>
{
	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,2,cnt_type> & div_c, ids_type * e)
	{
		return e[0] + div_c[0] * e[1];
	}

	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,2,cnt_type> & div_c, const grid_key_dx<2,ids_type> & e)
	{
		return e.get(0) + div_c[0] * e.get(1);
	}

	template<typename T> static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,2,cnt_type> & div_c,
			                                                       const openfpm::array<T,2,cnt_type> & spacing,
			                                                       const openfpm::array<ids_type,2,cnt_type> & off,
			                                                       const transform & t,
			                                                       const Point<2,T> & p)
	{

		return openfpm::math::uint_floor(t.transform(p,0)/spacing[0]) + off[0] +
				  (openfpm::math::uint_floor(t.transform(p,1)/spacing[1]) + off[1])*div_c[0];
	}

	template<typename T> static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,2,cnt_type> & div_c,
			                                                       const openfpm::array<T,2,cnt_type> & spacing,
			                                                       const openfpm::array<ids_type,2,cnt_type> & off,
			                                                       const transform & t,
			                                                       const T * p,
			                                                       ids_type * e)
	{
		e[0] = openfpm::math::uint_floor(t.transform(p,0)/spacing[0]) + off[0];
		e[1] = openfpm::math::uint_floor(t.transform(p,1)/spacing[1]) + off[1];

		return e[0] + e[1]*div_c[0];
	}

	template<typename T>
	static inline __device__ __host__ grid_key_dx<2,ids_type> get_cid_key(const openfpm::array<T,2,cnt_type> & spacing,
			                                                       const openfpm::array<ids_type,2,cnt_type> & off,
			                                                       const transform & t,
			                                                       const Point<2,T> & p)
	{
		grid_key_dx<2,ids_type> e;

		e.set_d(0,openfpm::math::uint_floor(t.transform(p,0)/spacing[0]) + off[0]);
		e.set_d(1,openfpm::math::uint_floor(t.transform(p,1)/spacing[1]) + off[1]);

		return e;
	}

	template <typename U = cnt_type, typename sfinae=typename std::enable_if<std::is_same<ids_type,U>::value >::type >
	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,2,cnt_type> & div_c,
			                                  const grid_key_dx<2,cnt_type> & e)
	{
		return e.get(0) + e.get(1)*div_c[0];
	}
};


template<typename cnt_type, typename ids_type,typename transform>
struct cid_<3,cnt_type,ids_type,transform>
{

	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,3,cnt_type> & div_c,
			                                  const ids_type * e)
	{
		return e[0] + (e[1] + e[2]*div_c[1])*div_c[0];
	}

	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,3,cnt_type> & div_c,
			                                  const grid_key_dx<3,ids_type> & e)
	{
		return e.get(0) + (e.get(1) + e.get(2)*div_c[1])*div_c[0];
	}

	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,3,cnt_type> & div_c,
													   const openfpm::array<ids_type,3,cnt_type> & off,
			                                  const grid_key_dx<3,ids_type> & e)
	{
		return (e.get(0) + off[0]) + ((e.get(1) + off[1]) + (e.get(2) + off[2])*div_c[1])*div_c[0];
	}

	template<typename T> static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,3,cnt_type> & div_c,
			                                                       const openfpm::array<T,3,cnt_type> & spacing,
			                                                       const openfpm::array<ids_type,3,cnt_type> & off,
			                                                       const transform & t,
			                                                       const Point<3,T> & p)
	{

		return openfpm::math::uint_floor(t.transform(p,0)/spacing[0]) + off[0] +
				  (openfpm::math::uint_floor(t.transform(p,1)/spacing[1]) + off[1] +
				  (openfpm::math::uint_floor(t.transform(p,2)/spacing[2]) + off[2])*div_c[1])*div_c[0];
	}

	template<typename T> static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,3,cnt_type> & div_c,
			                                                       const openfpm::array<T,3,cnt_type> & spacing,
			                                                       const openfpm::array<ids_type,3,cnt_type> & off,
			                                                       const transform & t,
			                                                       const T * p,
			                                                       ids_type * e)
	{
		e[0] = openfpm::math::uint_floor(t.transform(p,0)/spacing[0]) + off[0];
		e[1] = openfpm::math::uint_floor(t.transform(p,1)/spacing[1]) + off[1];
		e[2] = openfpm::math::uint_floor(t.transform(p,2)/spacing[2]) + off[2];

		return e[0] + (e[1] + e[2]*div_c[1])*div_c[0];
	}

	template<typename T>
	static inline __device__ __host__ grid_key_dx<3,ids_type> get_cid_key(const openfpm::array<T,3,cnt_type> & spacing,
			                                                       const openfpm::array<ids_type,3,cnt_type> & off,
			                                                       const transform & t,
			                                                       const Point<3,T> & p)
	{
		grid_key_dx<3,ids_type> e;

		e.set_d(0,openfpm::math::uint_floor(t.transform(p,0)/spacing[0]) + off[0]);
		e.set_d(1,openfpm::math::uint_floor(t.transform(p,1)/spacing[1]) + off[1]);
		e.set_d(2,openfpm::math::uint_floor(t.transform(p,2)/spacing[2]) + off[2]);

		return e;
	}

	template <typename U = cnt_type, typename sfinae=typename std::enable_if<std::is_same<ids_type,U>::value >::type >
	static inline __device__ __host__ cnt_type get_cid(const openfpm::array<ids_type,3,cnt_type> & div_c,
			                                  const grid_key_dx<3,cnt_type> & e)
	{
		return e.get(0) + (e.get(1) + e.get(2)*div_c[1])*div_c[0];
	}
};

template<unsigned int bit_phases, typename cnt_type>
struct shift_ph
{
	typedef boost::mpl::int_<sizeof(cnt_type) - bit_phases> type;
	typedef boost::low_bits_mask_t<sizeof(size_t)*8-bit_phases>  mask_low;
};


template<typename cnt_type, typename ph>
__device__ __host__ cnt_type encode_phase_id(cnt_type ph_id,cnt_type pid)
{
    return pid + (ph_id << ph::type::value);
}


/////////////////////////// THIS ONLY WORK IF NVCC IS COMPILING THIS //////////////////////////////////////////////////////

#ifdef __NVCC__

template<bool is_sparse,unsigned int dim, typename pos_type,
         typename cnt_type, typename ids_type, typename transform>
__global__ void subindex(openfpm::array<ids_type,dim,cnt_type> div,
						 openfpm::array<pos_type,dim,cnt_type> spacing,
						 openfpm::array<ids_type,dim,cnt_type> off,
						 transform t,
						 int n_cap,
						 int n_part,
						 int n_cap2,
						 cnt_type start,
						 pos_type * p_pos,
						 cnt_type *counts,
						 ids_type * p_ids)
{
    cnt_type i, cid, ins;
    ids_type e[dim+1];

    i = threadIdx.x + blockIdx.x * blockDim.x + start;
    ins = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= n_part) return;

    pos_type p[dim];

    for (size_t k = 0 ; k < dim ; k++)
    {p[k] = p_pos[i+k*n_cap];}

    cid = cid_<dim,cnt_type,ids_type,transform>::get_cid(div,spacing,off,t,p,e);

    if (is_sparse == false)
    {
    	e[dim] = atomicAdd(counts + cid, 1);

        for (size_t k = 0 ; k <= dim ; k++)
        {p_ids[ins+k*(n_cap2)] = e[k];}
    }
    else
    {
        for (size_t k = 0 ; k <= dim ; k++)
        {p_ids[ins+k*(n_cap2)] = e[k];}

        counts[ins] = cid;
    }
}

template<typename vector_sparse, typename vector_cell>
__global__ void fill_cells_sparse(vector_sparse vs, vector_cell vc)
{
	vs.init();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	if (p >= vc.size())	{return;}

	int c = vc.template get<0>(p);

	vs.template insert<0>(c) = p;
	vs.flush_block_insert();
}


template<unsigned int dim, typename pos_type, typename cnt_type, typename ids_type, typename transform>
__global__ void subindex_without_count(openfpm::array<ids_type,dim,cnt_type> div,
						 openfpm::array<pos_type,dim,cnt_type> spacing,
						 openfpm::array<ids_type,dim,cnt_type> off,
						 transform t,
						 int n_cap,
						 int n_part,
						 int n_cap2,
						 pos_type * p_pos,
						 cnt_type *counts,
						 ids_type * p_ids)
{
    cnt_type i, cid;
    ids_type e[dim+1];

    i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= n_part) return;

    pos_type p[dim];

    for (size_t k = 0 ; k < dim ; k++)
    {p[k] = p_pos[i+k*n_cap];}

    cid = cid_<dim,cnt_type,ids_type,transform>::get_cid(div,spacing,off,t,p,e);

    e[dim] = atomicAdd(counts + cid, 1);

    for (size_t k = 0 ; k <= dim ; k++)
    {p_ids[i+k*(n_cap2)] = e[k];}
}

template<unsigned int dim, typename cnt_type, typename ids_type, typename ph>
__global__ void fill_cells(cnt_type phase_id ,
		                   openfpm::array<ids_type,dim,cnt_type> div_c,
		                   openfpm::array<ids_type,dim,cnt_type> off,
		                   cnt_type n,
		                   cnt_type n_cap,
		                   cnt_type start_p,
		                   const cnt_type *starts,
		                   const ids_type * p_ids,
		                   cnt_type *cells)
{
    cnt_type i, cid, id, start;
    ids_type e[dim+1];

    i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= n) return;

#pragma unroll
    for (int j = 0 ; j < dim+1 ; j++)
    {e[j] = p_ids[j*n_cap+i];}

    cid = cid_<dim,cnt_type,ids_type,int>::get_cid(div_c, e);

    start = starts[cid];
    id = start + e[dim];

    cells[id] = encode_phase_id<cnt_type,ph>(phase_id,i + start_p);
}






template<typename cnt_type, typename ph>
__device__ inline void phase_and_id(cnt_type c, cnt_type *ph_id, cnt_type *pid)
{
    *pid = c & ph::mask_low::value;
    *ph_id   = c << ph::type::value;
}


/*template <unsigned int n_prp , typename cnt_type, typename T>
__device__ inline void reorderMP(const boost_array_openfpm<T*,n_prp> input,
		              boost_array_openfpm<T*,dim> output,
		              cnt_type i)
{
    int src_id, ph_id;
    phase_and_id(i, &ph_id, &src_id);
    output[ph_id][i] = src.d[ph_id][src_id];
}*/

template <typename vector_prp , typename cnt_type>
__device__ inline void reorder(const vector_prp & input,
		                       vector_prp & output,
		                       cnt_type src_id,
		                       cnt_type dst_id)
{
	output.set(dst_id,input,src_id);
}

template <typename vector_prp, typename vector_pos, typename vector_ns, typename cnt_type, typename sh>
__global__ void reorder_parts(int n,
		                      const vector_prp input,
		                      vector_prp output,
		                      const vector_pos input_pos,
		                      vector_pos output_pos,
		                      vector_ns sorted_non_sorted,
		                      vector_ns non_sorted_to_sorted,
		                      const cnt_type * cells)
{
    cnt_type i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= n) return;

    cnt_type code = cells[i];
    reorder(input, output, code,i);
    reorder(input_pos,output_pos,code,i);

    sorted_non_sorted.template get<0>(i) = code;
    non_sorted_to_sorted.template get<0>(code) = i;
}

template<typename vector_sort_index, typename vector_out_type>
__global__ void mark_domain_particles(vector_sort_index vsi, vector_out_type vout_dg, int g_m)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i >= vsi.size()) return;

	vout_dg.template get<0>(i) = (vsi.template get<0>(i) < g_m)?1:0;
}

template<typename scan_type, typename vector_out_type>
__global__ void collect_domain_ghost_ids(scan_type scan, vector_out_type vout_ids)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i >= scan.size()-1) return;

	auto pp = scan.template get<0>(i+1);
	auto p = scan.template get<0>(i);

	if (pp != p)
	{vout_ids.template get<0>(scan.template get<0>(i)) = i;}
}

template<typename cl_sparse_type, typename vector_type, typename vector_type2>
__global__ void count_nn_cells(cl_sparse_type cl_sparse, vector_type output, vector_type2 nn_to_test)
{
	typedef typename cl_sparse_type::index_type index_type;

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	openfpm::sparse_index<index_type> id;
	id.id = tid;

	if (tid >= cl_sparse.size()) {return;}

	index_type cell = cl_sparse.get_index(id);

	for (int i = 0 ; i < nn_to_test.size() ; i++)
	{
		index_type cell_n = cell + nn_to_test.template get<0>(i);

		index_type start = cl_sparse.template get<0>(cell_n);

		if (start != (index_type)-1)
		{
			// Cell exist
			output.template get<0>(tid) += 1;
		}
	}
};

template<typename cl_sparse_type, typename vector_type, typename vector_type2, typename vector_type3>
__global__ void fill_nn_cells(cl_sparse_type cl_sparse, vector_type starts, vector_type2 nn_to_test, vector_type3 output, typename cl_sparse_type::index_type max_stop)
{
	typedef typename cl_sparse_type::index_type index_type;

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	openfpm::sparse_index<index_type> id;
	id.id = tid;

	if (tid >= cl_sparse.size())	{return;}

	index_type cell = cl_sparse.get_index(id);

	int cnt = 0;
	for (int i = 0 ; i < nn_to_test.size() ; i++)
	{
		index_type cell_n = cell + nn_to_test.template get<0>(i);

		auto sid = cl_sparse.get_sparse(cell_n);

		if (sid.id != (index_type)cl_sparse.size())
		{
			index_type start = cl_sparse.template get<0>(sid);
			// Cell exist
			output.template get<0>(starts.template get<0>(tid) + cnt) = start;

			if (sid.id == cl_sparse.size() - 1)
			{output.template get<1>(starts.template get<0>(tid) + cnt) = max_stop;}
			else
			{
				decltype(sid) sid_t;
				sid_t.id = sid.id+1;
				output.template get<1>(starts.template get<0>(tid) + cnt) = cl_sparse.template get<0>(sid_t);
			}
			++cnt;
		}
	}
};

template<typename T>
struct to_type4
{
	typedef void type;
};

template<>
struct to_type4<unsigned int>
{
	typedef uint4 type;
};


template<>
struct to_type4<int>
{
	typedef int4 type;
};

template<>
struct to_type4<float>
{
	typedef float4 type;
};

/////////////////////////// THIS ONLY WORK IF NVCC IS COMPILING THIS //////////////////////////////////////////////////////

#endif

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CUDA_CUDA_CELL_LIST_UTIL_FUNC_HPP_ */
