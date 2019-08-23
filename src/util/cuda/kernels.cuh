/*
 * kernels.cuh
 *
 *  Created on: Jan 26, 2019
 *      Author: i-bird
 */

#ifndef KERNELS_CUH_
#define KERNELS_CUH_


template<unsigned int prp_off, typename vector_type,typename vector_type_offs>
__global__  void find_buffer_offsets_zero(vector_type vd, int * cnt, vector_type_offs offs)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= (int)vd.size()) return;

    if (p == 0)
    {
    	int i = atomicAdd(cnt, 1);
    	offs.template get<1>(i) = 0;
    	offs.template get<0>(i) = vd.template get<prp_off>(0);
    	return;
    }

    if (vd.template get<prp_off>(p-1) != vd.template get<prp_off>(p))
	{
    	int i = atomicAdd(cnt, 1);
    	offs.template get<1>(i) = p;
    	offs.template get<0>(i) = vd.template get<prp_off>(p);
	}
}

template<unsigned int prp_off, typename vector_type2, typename vector_type,typename vector_type_offs>
__global__  void construct_index_unique(vector_type2 vd_input, vector_type vd, vector_type_offs offs)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= (int)vd.size() - 1) return;

    unsigned int id = vd.template get<0>(p);
    unsigned int id_p1 = vd.template get<0>(p+1);

    if (p == 0)
    {
        offs.template get<0>(id) = vd_input.template get<prp_off>(0);
        offs.template get<1>(id) = p;
    }

    if (id != id_p1)
    {
    	offs.template get<0>(id) = vd_input.template get<prp_off>(p);
    	offs.template get<1>(id) = p;
    }
}

template<unsigned int prp_off, typename vector_type,typename vector_type_offs>
__global__  void find_buffer_offsets_for_scan(vector_type vd, vector_type_offs offs)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= (int)vd.size()) return;

    unsigned int pm1 = (p == 0)?p:p-1;

    bool predicate = vd.template get<prp_off>(pm1) != vd.template get<prp_off>(p) || (p == 0);

    offs.template get<0>(p) = predicate;
}


template<unsigned int prp_off, typename vector_type,typename vector_type_offs>
__global__  void find_buffer_offsets(vector_type vd, int * cnt, vector_type_offs offs)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= (int)vd.size() - 1) return;

    if (vd.template get<prp_off>(p) != vd.template get<prp_off>(p+1))
	{
    	int i = atomicAdd(cnt, 1);

    	offs.template get<0>(i) = p+1;
    	offs.template get<1>(i) = vd.template get<prp_off>(p);
	}
}

template<unsigned int prp_off, typename vector_type,typename vector_type_offs>
__global__  void find_buffer_offsets_no_prc(vector_type vd, int * cnt, vector_type_offs offs, int g_m)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= (int)g_m - 1) return;

    if (vd.template get<prp_off>(p) != vd.template get<prp_off>(p+1))
	{
    	int i = atomicAdd(cnt, 1);
    	offs.template get<0>(i) = p+1;
	}
}


#endif /* KERNELS_CUH_ */
