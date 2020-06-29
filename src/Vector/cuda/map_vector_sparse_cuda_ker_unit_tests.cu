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
#include <map>

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
	vd_insert.flush_block_insert();
}

//template<unsigned int blockLength, typename vd_type>
//__global__ void test_insert_sparse_block(vd_type vd_insert)
//{
//	vd_insert.init();
//
//	int p = blockIdx.x*blockDim.x + threadIdx.x;
//
//	p *= 2;
//
//	auto ie = vd_insert.insert(10000 - p);
//	ie.template get<0>() = p + 100;
//	for (unsigned int i = 0 ; i < blockLength ; i++)
//	{
//		ie.template get<1>()[i] = p + 10100 + i;
//	}
//
//	vd_insert.flush_block_insert();
//}

template<typename vd_type>
__global__ void test_remove_sparse(vd_type vd_insert)
{
	vd_insert.init();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	vd_insert.remove(10000 - p);
	vd_insert.flush_block_remove();
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

	vd_insert.flush_block_insert();
}

template<typename vd_type>
__global__ void test_remove_sparse2(vd_type vd_insert)
{
	vd_insert.init();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	vd_insert.remove(9000 - p);

	vd_insert.flush_block_remove();
}

template<typename vd_type>
__global__ void test_insert_sparse2_inc(vd_type vd_insert)
{
	vd_insert.init_ins_inc();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	auto ie = vd_insert.insert(9000 - p);
	ie.template get<0>() = p + 3000;
	ie.template get<1>() = p + 13000;
	ie.template get<2>() = p + 23000;

	vd_insert.flush_block_insert();
}

template<typename vd_type>
__global__ void test_remove_sparse2_inc(vd_type vd_insert)
{
	vd_insert.init_rem_inc();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	vd_insert.remove(9000 - p);

	vd_insert.flush_block_remove();
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

	vd_insert.flush_block_insert();
}

template<typename vd_type>
__global__ void test_remove_sparse3(vd_type vd_insert)
{
	vd_insert.init();

	int p = blockIdx.x*blockDim.x + threadIdx.x;

	p *= 2;

	vd_insert.remove(p);

	vd_insert.flush_block_remove();
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

//template<unsigned int blockLength, typename vd_sparse_type, typename vector_out_type>
//__global__ void test_sparse_get_test_block(vd_sparse_type vd_test, vector_out_type output)
//{
//	int p = blockIdx.x*blockDim.x + threadIdx.x;
//	int i = blockIdx.x*blockDim.x + threadIdx.x;
//
//	p *= 2;
//
//	int v;
//
//	output.template get<0>(i) = vd_test.template get<0>(10000 - p,v);
//	for (int j=0; j<blockLength; ++j)
//	{
//		output.template get<1>(i)[j] = vd_test.template get_ele<1>(v)[j];
//	}
//}

BOOST_AUTO_TEST_SUITE( vector_cuda_sparse )

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_gpu )
{
	openfpm::vector_sparse_gpu<aggregate<size_t,float,double>> vs;

	vs.template setBackground<0>(17);

	vs.template setBackground<1>(18);

	vs.template setBackground<2>(19);

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());

	mgpu::ofp_context_t ctx;
	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);

	vs.setGPUInsertBuffer(10,1024);
	test_insert_sparse2<<<10,100>>>(vs.toKernel());

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);

	vs.setGPUInsertBuffer(4000,512);
	test_insert_sparse3<<<4000,256>>>(vs.toKernel());

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);

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


		if (match == false){break;}
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

//BOOST_AUTO_TEST_CASE( vector_sparse_cuda_gpu_block )
//{
//	openfpm::vector_sparse_gpu<aggregate<size_t, DataBlock>> vs;
//
//	// Set the background values
//	vs.template getBackground<0>() = 17;
//
//	for (int i=0; i<DataBlock::length; ++i)
//	{
//		vs.template getBackground<1>()[i] = 666;
//	}
//
//	const unsigned int gridSize = 500;
//	const unsigned int blockSizeInsert = 128;
//	const unsigned int blockSizeRead = 256;
//
//	// Prealloc insertion buffer
//	vs.setGPUInsertBuffer(gridSize,1024);
//
//	// Insert some data on the vector_sparse_gpu
//	test_insert_sparse_block<DataBlock::length><<<gridSize,blockSizeInsert>>>(vs.toKernel());
//
//
//	mgpu::ofp_context_t ctx;
//
//	// Flushing the inserts with some reduction operator
//	vs.flush<sadd_<0>, smax_block_<1, DataBlock::length>>(ctx,flush_type::FLUSH_ON_DEVICE);
//
//	openfpm::vector_gpu<aggregate<size_t,DataBlock>> output;
//
//	output.resize(gridSize*blockSizeRead);
//
//	// Copy results to an output vector
//	test_sparse_get_test_block<DataBlock::length><<<gridSize,blockSizeRead>>>(vs.toKernel(),output.toKernel());
//
//	output.template deviceToHost<0,1>();
//	vs.template deviceToHost<0,1>();
//
//	bool match = true;
//	for (size_t i = 0 ; i < output.size()  ; i++)
//	{
//		match &= output.template get<0>(i) == vs.template get<0>(10000 - 2*i);
//		// std::cout << "SCALAR " << output.template get<0>(i) << std::endl;
//		// std::cout << i;
//		for (int j=0; j<DataBlock::length; ++j)
//		{
//			// std::cout << " " <<  output.template get<1>(i)[j] << "  " << vs.template get<1>(10000 - 2*i)[j] << ", ";
//			match &= output.template get<1>(i)[j] == vs.template get<1>(10000 - 2*i)[j];
//		}
//		// std::cout << std::endl;
//	}
//
//	BOOST_REQUIRE_EQUAL(match,true);
//
//	vs.clear();
//}

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_gpu_incremental_add )
{
	openfpm::vector_sparse_gpu<aggregate<size_t,float,double>> vs;

	vs.template setBackground<0>(17);

	vs.template setBackground<1>(18);

	vs.template setBackground<2>(19);

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());

	mgpu::ofp_context_t ctx;

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);

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

	vs.template setBackground<0>(0);

	vs.template setBackground<1>(0);

	vs.template setBackground<2>(0);

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());


	mgpu::ofp_context_t ctx;
	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);

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

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);

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

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);
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

	vs.template setBackground<0>(17);

	vs.template setBackground<1>(18);

	vs.template setBackground<2>(19);

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());

	mgpu::ofp_context_t ctx;

	vs.flush<sstart_<0>>(ctx,flush_type::FLUSH_ON_DEVICE);

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

//////////////////////////////// REMOVE test section

void check_lines(openfpm::vector_sparse_gpu<aggregate<size_t,float,double>> & vs,
		         bool s1,
		         bool s2,
		         bool s3,
		         bool s4)
{
	bool match = true;

	for (size_t i = 0 ; i <= 3500 ; i++)
	{
		if (s1 == true)
		{
			match &= vs.template get<0>(2*i) == 5;
			match &= vs.template get<1>(2*i) == 1;
			match &= vs.template get<2>(2*i) == 1;
		}
		else
		{
			match &= vs.template get<0>(2*i) == 17;
			match &= vs.template get<1>(2*i) == 18;
			match &= vs.template get<2>(2*i) == 19;
		}
	}

	BOOST_REQUIRE_EQUAL(match,true);

	for (size_t i = 3501 ; i <= 4000 ; i++)
	{
		if (s2 == true)
		{
			match &= vs.template get<0>(2*i) == 5 - 2*i + 3000 + 9000;
			match &= vs.template get<1>(2*i) == 1;
			match &= vs.template get<2>(2*i) == 23000 + 9000 - 2*i;
		}
		else
		{
			match &= vs.template get<0>(2*i) == 17;
			match &= vs.template get<1>(2*i) == 18;
			match &= vs.template get<2>(2*i) == 19;
		}
	}

	BOOST_REQUIRE_EQUAL(match,true);

	for (size_t i = 4001 ; i <= 4500 ; i++)
	{
		if (s3 == true)
		{
			match &= vs.template get<0>(2*i) == 5 - 2*i + 1100 - 2*i + 3000 + 18000;
			match &= vs.template get<1>(2*i) == 1;
			match &= vs.template get<2>(2*i) == 23000 + 9000 - 2*i;
		}
		else
		{
			match &= vs.template get<0>(2*i) == 17;
			match &= vs.template get<1>(2*i) == 18;
			match &= vs.template get<2>(2*i) == 19;
		}
	}

	BOOST_REQUIRE_EQUAL(match,true);

	for (size_t i = 4501 ; i <= 5000 ; i++)
	{
		if (s4 == true)
		{
			match &= vs.template get<0>(2*i) == 5 - 2*i + 1100 + 9000;
			match &= vs.template get<1>(2*i) == 1;
			match &= vs.template get<2>(2*i) == 21100 + 9000 - 2*i;
		}
		else
		{
			match &= vs.template get<0>(2*i) == 17;
			match &= vs.template get<1>(2*i) == 18;
			match &= vs.template get<2>(2*i) == 19;
		}
	}
}


BOOST_AUTO_TEST_CASE( vector_sparse_cuda_gpu_remove )
{
	openfpm::vector_sparse_gpu<aggregate<size_t,float,double>> vs;

	vs.template setBackground<0>(17);

	vs.template setBackground<1>(18);

	vs.template setBackground<2>(19);

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());

	mgpu::ofp_context_t ctx;
	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);

	vs.setGPUInsertBuffer(10,1024);
	test_insert_sparse2<<<10,100>>>(vs.toKernel());

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);

	vs.setGPUInsertBuffer(4000,512);
	test_insert_sparse3<<<4000,256>>>(vs.toKernel());

	vs.flush<sadd_<0>,smin_<1>,smax_<2> >(ctx,flush_type::FLUSH_ON_DEVICE);

	// we launch a kernel to insert data
	vs.setGPURemoveBuffer(10,1024);
	test_remove_sparse<<<10,100>>>(vs.toKernel());

	size_t sz = vs.size();

	vs.flush_remove(ctx,flush_type::FLUSH_ON_DEVICE);
	vs.template deviceToHost<0,1,2>();

	BOOST_REQUIRE_EQUAL(vs.size(),sz - 1000);

	check_lines(vs,
			    true,
			    true,
			    false,
			    false);

	// we launch a kernel to insert data
	vs.setGPURemoveBuffer(10,1024);
	test_remove_sparse2<<<10,100>>>(vs.toKernel());

	vs.flush_remove(ctx,flush_type::FLUSH_ON_DEVICE);

	BOOST_REQUIRE_EQUAL(vs.size(),sz - 1500);

	vs.template deviceToHost<0,1,2>();

	check_lines(vs,
			    true,
			    false,
			    false,
			    false);

	vs.setGPURemoveBuffer(4000,512);
	test_remove_sparse3<<<4000,256>>>(vs.toKernel());

	vs.flush_remove(ctx,flush_type::FLUSH_ON_DEVICE);

	BOOST_REQUIRE_EQUAL(vs.size(),0);

	vs.template deviceToHost<0,1,2>();

	BOOST_REQUIRE_EQUAL(vs.size(),0);

	check_lines(vs,
			    false,
			    false,
			    false,
			    false);
}

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_gpu_remove_incremental )
{
	openfpm::vector_sparse_gpu<aggregate<size_t,float,double>> vs;

	vs.template setBackground<0>(17);

	vs.template setBackground<1>(18);

	vs.template setBackground<2>(19);

	vs.setGPUInsertBuffer(10,1024);

	// we launch a kernel to insert data
	test_insert_sparse<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());
	test_insert_sparse2_inc<<<10,100>>>(vs.toKernel());

	mgpu::ofp_context_t ctx;

	vs.flush<sadd_<0>,sadd_<1>,sadd_<2>>(ctx,flush_type::FLUSH_ON_DEVICE);

	// we launch a kernel to insert data
	vs.setGPURemoveBuffer(10,1024);
	test_remove_sparse<<<10,100>>>(vs.toKernel());
	test_remove_sparse2_inc<<<10,99>>>(vs.toKernel());
	test_remove_sparse2_inc<<<10,99>>>(vs.toKernel());
	test_remove_sparse2_inc<<<10,99>>>(vs.toKernel());

	vs.flush_remove(ctx,flush_type::FLUSH_ON_DEVICE);

	BOOST_REQUIRE_EQUAL(vs.size(),10);

	vs.template deviceToHost<0,1,2>();

	BOOST_REQUIRE_EQUAL(vs.template get<0>(7020),14940);
	BOOST_REQUIRE_EQUAL(vs.template get<0>(7018),14946);
	BOOST_REQUIRE_EQUAL(vs.template get<0>(7016),14952);
	BOOST_REQUIRE_EQUAL(vs.template get<0>(7014),14958);
	BOOST_REQUIRE_EQUAL(vs.template get<0>(7012),14964);
	BOOST_REQUIRE_EQUAL(vs.template get<0>(7010),14970);
	BOOST_REQUIRE_EQUAL(vs.template get<0>(7008),14976);
	BOOST_REQUIRE_EQUAL(vs.template get<0>(7006),14982);
	BOOST_REQUIRE_EQUAL(vs.template get<0>(7004),14988);
	BOOST_REQUIRE_EQUAL(vs.template get<0>(7002),14994);

	BOOST_REQUIRE_EQUAL(vs.template get<1>(7020),44940);
	BOOST_REQUIRE_EQUAL(vs.template get<1>(7018),44946);
	BOOST_REQUIRE_EQUAL(vs.template get<1>(7016),44952);
	BOOST_REQUIRE_EQUAL(vs.template get<1>(7014),44958);
	BOOST_REQUIRE_EQUAL(vs.template get<1>(7012),44964);
	BOOST_REQUIRE_EQUAL(vs.template get<1>(7010),44970);
	BOOST_REQUIRE_EQUAL(vs.template get<1>(7008),44976);
	BOOST_REQUIRE_EQUAL(vs.template get<1>(7006),44982);
	BOOST_REQUIRE_EQUAL(vs.template get<1>(7004),44988);
	BOOST_REQUIRE_EQUAL(vs.template get<1>(7002),44994);

	BOOST_REQUIRE_EQUAL(vs.template get<2>(7020),74940);
	BOOST_REQUIRE_EQUAL(vs.template get<2>(7018),74946);
	BOOST_REQUIRE_EQUAL(vs.template get<2>(7016),74952);
	BOOST_REQUIRE_EQUAL(vs.template get<2>(7014),74958);
	BOOST_REQUIRE_EQUAL(vs.template get<2>(7012),74964);
	BOOST_REQUIRE_EQUAL(vs.template get<2>(7010),74970);
	BOOST_REQUIRE_EQUAL(vs.template get<2>(7008),74976);
	BOOST_REQUIRE_EQUAL(vs.template get<2>(7006),74982);
	BOOST_REQUIRE_EQUAL(vs.template get<2>(7004),74988);
	BOOST_REQUIRE_EQUAL(vs.template get<2>(7002),74994);
}


BOOST_AUTO_TEST_SUITE_END()

#endif /* MAP_VECTOR_SPARSE_CUDA_KER_UNIT_TESTS_CUH_ */
