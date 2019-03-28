/*
 * map_vector_sparse_cuda_ker_unit_tests.cuh
 *
 *  Created on: Jan 24, 2019
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_SPARSE_CUDA_KER_UNIT_TESTS_CUH_
#define MAP_VECTOR_SPARSE_CUDA_KER_UNIT_TESTS_CUH_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Vector/map_vector_sparse.hpp"
#include "Vector/map_vector.hpp"

template<typename vd_type>
__global__ void test_insert_sparse(vd_type vd_insert)
{
	vd_insert.init();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	auto ie = vd_insert.insert(10000 - p);
	ie.template get<0>() = p + 100;
	ie.template get<1>() = p + 10100;
	ie.template get<2>() = p + 20100;
	vd_insert.flush_block();
}

template<typename vd_type>
__global__ void test_insert_sparse2(vd_type vd_insert)
{
	vd_insert.init();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	auto ie = vd_insert.insert(9000 - p);
	ie.template get<0>() = p + 3000;
	ie.template get<1>() = p + 13000;
	ie.template get<2>() = p + 23000;

	vd_insert.flush_block();
}

template<typename vd_type>
__global__ void test_insert_sparse2_inc(vd_type vd_insert)
{
	vd_insert.init_inc();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	auto ie = vd_insert.insert(9000 - p);
	ie.template get<0>() = p + 3000;
	ie.template get<1>() = p + 13000;
	ie.template get<2>() = p + 23000;

	vd_insert.flush_block();
}

template<typename vd_type>
__global__ void test_insert_sparse3(vd_type vd_insert)
{
	vd_insert.init();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	auto ie = vd_insert.insert(p);
	ie.template get<0>() = 5;
	ie.template get<1>() = 1;
	ie.template get<2>() = 1;

	vd_insert.flush_block();
}

template<typename vd_sparse_type, typename vector_out_type>
__global__ void test_sparse_get_test(vd_sparse_type vd_test, vector_out_type output)
{
	int p = blockIdx.x*blockDim.x + threadIdx.x;
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	int v;

	output.template get<0>(i) = vd_test.template get<0>(10000 - p,v);
	output.template get<1>(i) = vd_test.template get_ele<1>(v);
	output.template get<2>(i) = vd_test.template get_ele<2>(v);
}

BOOST_AUTO_TEST_SUITE( vector_cuda_sparse )

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_gpu )
{
	openfpm::vector_sparse_gpu<aggregate<size_t,float,double>> vs;

	vs.template getBackground<0>() = 17;

	vs.template getBackground<1>() = 18;

	vs.template getBackground<2>() = 19;

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());


	mgpu::ofp_context_t ctx;
	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flust_type::FLUSH_ON_DEVICE);

	vs.setGPUInsertBuffer(10,1024);
	test_insert_sparse2<<<10,100>>>(vs.toKernel());

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flust_type::FLUSH_ON_DEVICE);

	vs.setGPUInsertBuffer(4000,512);
	test_insert_sparse3<<<4000,256>>>(vs.toKernel());

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flust_type::FLUSH_ON_DEVICE,1);

	openfpm::vector_gpu<aggregate<size_t,float,double>> output;

	output.resize(1500);

	test_sparse_get_test<<<10,150>>>(vs.toKernel(),output.toKernel());

	output.template deviceToHost<0,1,2>();
	vs.template deviceToHost<0,1,2>();

	bool match = true;
	for (size_t i = 0 ; i < output.size()  ; i++)
	{
		match &= output.template get<0>(i) == vs.template get<0>(10000 - 2*i);
		match &= output.template get<1>(i) == vs.template get<1>(10000 - 2*i);
		match &= output.template get<2>(i) == vs.template get<2>(10000 - 2*i);
	}

	BOOST_REQUIRE_EQUAL(match,true);

	vs.clear();

	test_sparse_get_test<<<10,150>>>(vs.toKernel(),output.toKernel());

	output.template deviceToHost<0,1,2>();
	vs.template deviceToHost<0,1,2>();

	match = true;
	for (size_t i = 0 ; i < output.size()  ; i++)
	{
		match &= output.template get<0>(i) == 17;
		match &= output.template get<1>(i) == 18;
		match &= output.template get<2>(i) == 19;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}


BOOST_AUTO_TEST_CASE( vector_sparse_cuda_gpu_incremental_add )
{
	openfpm::vector_sparse_gpu<aggregate<size_t,float,double>> vs;

	vs.template getBackground<0>() = 17;

	vs.template getBackground<1>() = 18;

	vs.template getBackground<2>() = 19;

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());

	mgpu::ofp_context_t ctx;

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flust_type::FLUSH_ON_DEVICE);

	vs.template deviceToHost<0,1,2>();

	BOOST_REQUIRE_EQUAL(vs.size(),1500);

	bool match = true;
	for (size_t i = 500 ; i < 1000 ; i++)
	{
		match &= vs.template get<0>(9000 - 2*i) == 3*(2*i + 3000);
		match &= vs.template get<1>(9000 - 2*i) == 2*i + 13000;
		match &= vs.template get<2>(9000 - 2*i) == 2*i + 23000;
	}

	for (size_t i = 0 ; i < 500 ; i++)
	{
		match &= vs.template get<0>(9000 - 2*i) == 3*(2*i + 3000) + 2*i + 1100;
		match &= vs.template get<1>(9000 - 2*i) == 2*i + 11100;
		match &= vs.template get<2>(9000 - 2*i) == 2*i + 23000;
	}

	for (size_t i = 0 ; i < 500 ; i++)
	{
		match &= vs.template get<0>(10000 - 2*i) == 2*i + 100;
		match &= vs.template get<1>(10000 - 2*i) == 2*i + 10100;
		match &= vs.template get<2>(10000 - 2*i) == 2*i + 20100;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_gpu_get )
{
	openfpm::vector_sparse_gpu<aggregate<size_t,float,double>> vs;

	vs.template getBackground<0>() = 0;

	vs.template getBackground<1>() = 0;

	vs.template getBackground<2>() = 0;

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());


	mgpu::ofp_context_t ctx;
	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flust_type::FLUSH_ON_DEVICE);

	vs.template deviceToHost<0,1,2>();

	bool match = true;
	for (size_t i = 0 ; i < 1000 ; i++)
	{
		match &= vs.template get<0>(10000 - 2*i) == 2*i + 100;
		match &= vs.template get<1>(10000 - 2*i) == 2*i + 10100;
		match &= vs.template get<2>(10000 - 2*i) == 2*i + 20100;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	vs.setGPUInsertBuffer(10,1024);
	test_insert_sparse2<<<10,100>>>(vs.toKernel());

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flust_type::FLUSH_ON_DEVICE);

	vs.template deviceToHost<0,1,2>();

	BOOST_REQUIRE_EQUAL(vs.size(),1500);

	match = true;
	for (size_t i = 500 ; i < 1000 ; i++)
	{
		match &= vs.template get<0>(9000 - 2*i) == 2*i + 3000;
		match &= vs.template get<1>(9000 - 2*i) == 2*i + 13000;
		match &= vs.template get<2>(9000 - 2*i) == 2*i + 23000;
	}

	for (size_t i = 0 ; i < 500 ; i++)
	{
		match &= vs.template get<0>(9000 - 2*i) == 2*i + 3000 + 2*i + 1100;
		match &= vs.template get<1>(9000 - 2*i) == 2*i + 11100;
		match &= vs.template get<2>(9000 - 2*i) == 2*i + 23000;

		if (match == false)
		{
			std::cout << vs.template get<2>(9000 - 2*i) << "   " << 2*i + 23000 << std::endl;
		}
	}

	for (size_t i = 0 ; i < 500 ; i++)
	{
		match &= vs.template get<0>(10000 - 2*i) == 2*i + 100;
		match &= vs.template get<1>(10000 - 2*i) == 2*i + 10100;
		match &= vs.template get<2>(10000 - 2*i) == 2*i + 20100;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	vs.setGPUInsertBuffer(4000,512);
	test_insert_sparse3<<<4000,256>>>(vs.toKernel());

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flust_type::FLUSH_ON_DEVICE,1);
	vs.template deviceToHost<0,1,2>();

	for (size_t i = 0 ; i <= 3500 ; i++)
	{
		match &= vs.template get<0>(2*i) == 5;
		match &= vs.template get<1>(2*i) == 1;
		match &= vs.template get<2>(2*i) == 1;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	for (size_t i = 3501 ; i <= 4000 ; i++)
	{
		match &= vs.template get<0>(2*i) == 5 - 2*i + 3000 + 9000;
		match &= vs.template get<1>(2*i) == 1;
		match &= vs.template get<2>(2*i) == 23000 + 9000 - 2*i;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	for (size_t i = 4001 ; i <= 4500 ; i++)
	{
		match &= vs.template get<0>(2*i) == 5 - 2*i + 1100 - 2*i + 3000 + 18000;
		match &= vs.template get<1>(2*i) == 1;
		match &= vs.template get<2>(2*i) == 23000 + 9000 - 2*i;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	for (size_t i = 4501 ; i <= 5000 ; i++)
	{
		match &= vs.template get<0>(2*i) == 5 - 2*i + 1100 + 9000;
		match &= vs.template get<1>(2*i) == 1;
		match &= vs.template get<2>(2*i) == 21100 + 9000 - 2*i;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_gpu_special_function )
{
	openfpm::vector_sparse_gpu<aggregate<size_t,float,double>> vs;

	vs.template getBackground<0>() = 17;

	vs.template getBackground<1>() = 18;

	vs.template getBackground<2>() = 19;

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());

	mgpu::ofp_context_t ctx;

	vs.flush<sstart_<0>>(ctx,flust_type::FLUSH_ON_DEVICE);

	vs.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(vs.size(),1500);

	bool match = true;

	size_t count = 0;
	for (size_t i = 999 ; i >= 500 ; i--)
	{
		match &= vs.template get<0>(9000 - 2*i) == count;
		count += 3;
	}

	for (long int i = 499 ; i >= 0 ; i--)
	{
		match &= vs.template get<0>(9000 - 2*i) == count;
		count += 4;
	}

	for (long int i = 499 ; i >= 0 ; i--)
	{
		match &= vs.template get<0>(10000 - 2*i) == count;
		count += 1;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}


BOOST_AUTO_TEST_SUITE_END()

#endif /* MAP_VECTOR_SPARSE_CUDA_KER_UNIT_TESTS_CUH_ */
