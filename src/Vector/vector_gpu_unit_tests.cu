/*
 * map_vector_cuda_funcs_tests.cu
 *
 *  Created on: Aug 17, 2018
 *      Author: i-bird
 */


#include "util/cuda_util.hpp"
#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Vector/map_vector.hpp"

template<typename vector_vector_type, typename vector_out_type>
__global__ void vv_test_size(vector_vector_type vvt, vector_out_type out)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vvt.size()) return;

	out.template get<0>(p) = vvt.template get<0>(p).size();
}

template<typename vector_vector_type, typename vector_out_type>
__global__ void vv_test_pointer(vector_vector_type vvt, vector_out_type out)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vvt.size()) return;

	out.template get<0>(p) = (size_t)vvt.template get<0>(p).template getPointer<0>();
	out.template get<1>(p) = (size_t)vvt.template get<0>(p).template getPointer<1>();
}

template<typename vector_vector_type, typename vector_out_type>
__global__ void vv_test_data_get(vector_vector_type vvt, vector_out_type out, int i_sz)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= out.size()) return;

	int id1 = p/i_sz;
	int id2 = p%i_sz;

	out.template get<0>(p)[0] = (size_t)vvt.template get<0>(id1).template get<0>(id2)[0];
	out.template get<0>(p)[1] = (size_t)vvt.template get<0>(id1).template get<0>(id2)[1];
	out.template get<0>(p)[2] = (size_t)vvt.template get<0>(id1).template get<0>(id2)[2];

	out.template get<1>(p)[0] = (size_t)vvt.template get<0>(id1).template get<1>(id2)[0];
	out.template get<1>(p)[1] = (size_t)vvt.template get<0>(id1).template get<1>(id2)[1];
	out.template get<1>(p)[2] = (size_t)vvt.template get<0>(id1).template get<1>(id2)[2];
}

BOOST_AUTO_TEST_SUITE( vector_cuda_tests )



BOOST_AUTO_TEST_CASE ( test_vector_of_vector_gpu )
{
	typedef openfpm::vector<Box<3,float>,CudaMemory,memory_traits_inte> proc_boxes;

	openfpm::vector<aggregate<proc_boxes>,CudaMemory,memory_traits_inte> vb_int_proc;

	vb_int_proc.resize_no_device(5);

	openfpm::vector<std::pair<void *,void *>> ptr_dev;
	ptr_dev.resize(5);

	for (size_t i = 0 ; i< vb_int_proc.size() ; i++)
	{
		vb_int_proc.template get<0>(i).resize(7);

		for (size_t j = 0 ; j < vb_int_proc.template get<0>(i).size() ; j++)
		{
			for (size_t k = 0 ; k < 3 ; k++)
			{
				vb_int_proc.template get<0>(i).template get<0>(j)[k] = i+j+k;
				vb_int_proc.template get<0>(i).template get<1>(j)[k] = 100+i+j+k;
			}
		}

		vb_int_proc.template get<0>(i).template hostToDevice<0,1>();
		ptr_dev.get(i).first = vb_int_proc.template get<0>(i).template getDeviceBuffer<0>();
		ptr_dev.get(i).second = vb_int_proc.template get<0>(i).template getDeviceBuffer<1>();
	}

	vb_int_proc.template hostToDevice<0>();

	openfpm::vector_gpu<aggregate<unsigned int>> out;
	out.resize(vb_int_proc.size());

	auto ite = vb_int_proc.getGPUIterator();

	CUDA_LAUNCH_DIM3((vv_test_size<decltype(vb_int_proc.toKernel()),decltype(out.toKernel())>),ite.wthr,ite.thr,vb_int_proc.toKernel(),out.toKernel());

	out.deviceToHost<0>();

	for (size_t i = 0 ; i < out.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(out.template get<0>(i),7);
	}

	openfpm::vector_gpu<aggregate<size_t,size_t>> out_pointer;
	out_pointer.resize(vb_int_proc.size());

	CUDA_LAUNCH_DIM3((vv_test_pointer<decltype(vb_int_proc.toKernel()),decltype(out_pointer.toKernel())>),ite.wthr,ite.thr,vb_int_proc.toKernel(),out_pointer.toKernel());

	out_pointer.deviceToHost<0,1>();

	for (size_t i = 0 ; i < out_pointer.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL((size_t)out_pointer.template get<0>(i),(size_t)ptr_dev.get(i).first);
		BOOST_REQUIRE_EQUAL((size_t)out_pointer.template get<1>(i),(size_t)ptr_dev.get(i).second);
	}

	openfpm::vector_gpu<aggregate<float[3],float[3]>> out_data;
	out_data.resize(vb_int_proc.size()*7);

	auto ite2 = out_data.getGPUIterator();

	CUDA_LAUNCH_DIM3((vv_test_data_get<decltype(vb_int_proc.toKernel()),decltype(out_data.toKernel())>),ite2.wthr,ite2.thr,vb_int_proc.toKernel(),out_data.toKernel(),7);

	out_data.template deviceToHost<0,1>();

	size_t i_sz = 7;

	for (size_t p = 0 ; p < out_data.size() ; p++)
	{
		int id1 = p/i_sz;
		int id2 = p%i_sz;

		BOOST_REQUIRE_EQUAL(out_data.template get<0>(p)[0],vb_int_proc.template get<0>(id1).template get<0>(id2)[0] );
		BOOST_REQUIRE_EQUAL(out_data.template get<0>(p)[1],vb_int_proc.template get<0>(id1).template get<0>(id2)[1] );
		BOOST_REQUIRE_EQUAL(out_data.template get<0>(p)[2],vb_int_proc.template get<0>(id1).template get<0>(id2)[2] );

		BOOST_REQUIRE_EQUAL(out_data.template get<1>(p)[0],vb_int_proc.template get<0>(id1).template get<1>(id2)[0] );
		BOOST_REQUIRE_EQUAL(out_data.template get<1>(p)[1],vb_int_proc.template get<0>(id1).template get<1>(id2)[1] );
		BOOST_REQUIRE_EQUAL(out_data.template get<1>(p)[2],vb_int_proc.template get<0>(id1).template get<1>(id2)[2] );
	}
}


BOOST_AUTO_TEST_SUITE_END()
