/*
 * map_vector_cuda_funcs_tests.cu
 *
 *  Created on: Aug 17, 2018
 *      Author: i-bird
 */


#define BOOST_GPU_ENABLED __host__ __device__

#include "util/cuda_launch.hpp"

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "util/cuda_util.hpp"
#include "Vector/map_vector.hpp"
#include "util/tokernel_transformation.hpp"

BOOST_AUTO_TEST_SUITE( vector_cuda_funcs_tests )


BOOST_AUTO_TEST_CASE( vector_cuda_funcs_add_prp_device )
{
	openfpm::vector_gpu<aggregate<float,float[3],float[3][3]>> vg_data;
	openfpm::vector_gpu<aggregate<float,float[3],float[3][3]>> vg_data2;

	vg_data.resize(100);
	vg_data2.resize(100);

	// we fill vg_data with something

	for (size_t i = 0 ; i < 100 ; i++)
	{
		vg_data.template get<0>(i) = 2.5 + i;

		vg_data.template get<1>(i)[0] = 4.6 + i;
		vg_data.template get<1>(i)[1] = 7.8 + i;
		vg_data.template get<1>(i)[2] = 9.0 + i;

		vg_data2.template get<0>(i) = 8.5 + i;

		vg_data2.template get<1>(i)[0] = 1.6 + i;
		vg_data2.template get<1>(i)[1] = 3.8 + i;
		vg_data2.template get<1>(i)[2] = 5.1 + i;
	}

	vg_data.hostToDevice<0,1>();
	vg_data2.hostToDevice<0,1>();

	vg_data.add_prp_device<aggregate<float,float[3],float[3][3]>,
						   CudaMemory,
						   openfpm::grow_policy_double,
						   OPENFPM_NATIVE,
						   memory_traits_inte,
	                       0,1>(vg_data2);

	vg_data.deviceToHost<0,1>();

	BOOST_REQUIRE_EQUAL(vg_data.size(),200);

	bool match = true;
	for (unsigned int i = 100 ; i < 200 ; i++)
	{
		match &= vg_data.template get<0>(i) == vg_data2.template get<0>(i-100);

		match &= vg_data.template get<1>(i)[0] == vg_data2.template get<1>(i-100)[0];
		match &= vg_data.template get<1>(i)[1] == vg_data2.template get<1>(i-100)[1];
		match &= vg_data.template get<1>(i)[2] == vg_data2.template get<1>(i-100)[2];
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( vector_cuda_to_kernel_recursive2 )
{
	typedef openfpm::vector_gpu<aggregate<int,long int>> test1_type;
	typedef openfpm::vector_gpu<aggregate<int,openfpm::vector_gpu<aggregate<long int>>>> test2_type;
	typedef openfpm::vector_gpu<aggregate<int,openfpm::vector_gpu<aggregate<Box<2,float>>>>> test3_type;
	typedef openfpm::vector<Box<3,float>,CudaMemory,memory_traits_inte> test4_type;

	typedef typename toKernel_transform<memory_traits_inte,test1_type>::type tker1;
	typedef typename toKernel_transform<memory_traits_inte,test2_type>::type tker2;
	typedef typename toKernel_transform<memory_traits_inte,test3_type>::type tker3;
	typedef typename toKernel_transform<memory_traits_inte,test4_type>::type tker4;

	bool test = std::is_same<tker1,openfpm::vector_gpu_ker<aggregate<int, long>, memory_traits_inte>>::value;

	BOOST_REQUIRE_EQUAL(test,true);

	test = std::is_same<tker2,openfpm::vector_gpu_ker<aggregate<int, openfpm::vector_gpu_ker<aggregate<long>, memory_traits_inte> >, memory_traits_inte>>::value;

	BOOST_REQUIRE_EQUAL(test,true);

	test = std::is_same<tker3,openfpm::vector_gpu_ker<aggregate<int, openfpm::vector_gpu_ker<aggregate<Box<2,float>>, memory_traits_inte> >, memory_traits_inte>>::value;

	BOOST_REQUIRE_EQUAL(test,true);

	test = std::is_same<tker4,openfpm::vector_gpu_ker<Box<3,float>,memory_traits_inte>>::value;

	BOOST_REQUIRE_EQUAL(test,true);
}

template<typename vv_rc,typename vector_output_type>
__global__ void kernel_recursive_check(vv_rc vvrc, vector_output_type vot)
{
	int k = 0;
	for (int i = 0 ; i < vvrc.size() ; i++)
	{
		for (int j = 0 ; j < vvrc.template get<1>(i).size() ; j++)
		{
			vot.template get<0>(k) = vvrc.template get<1>(i).template get<0>(j);
			k++;
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_cuda_to_kernel_recursive2_test_toKernel )
{
	typedef openfpm::vector_gpu<aggregate<int,openfpm::vector_gpu<aggregate<long int>>>> test2_type;
	typedef openfpm::vector_gpu<aggregate<int,openfpm::vector_gpu<aggregate<Box<2,float>>>>> test3_type;

	test2_type tt2;
	test3_type tt3;

	tt2.add_no_device();
	tt2.add_no_device();
	tt2.add_no_device();

/*	tt3.add();
	tt3.add();
	tt3.add();*/

	tt2.template get<0>(0) = 80;
	tt2.template get<1>(0).add();
	tt2.template get<1>(0).template get<0>(0) = 500;
	tt2.template get<0>(0) = 180;
	tt2.template get<1>(0).add();
	tt2.template get<1>(0).template get<0>(1) = 600;
	tt2.template get<0>(0) = 280;;
	tt2.template get<1>(0).add();
	tt2.template get<1>(0).template get<0>(2) = 700;
	tt2.template get<1>(0).template hostToDevice<0>();

	tt2.template get<0>(1) = 10080;
	tt2.template get<1>(1).add();
	tt2.template get<1>(1).template get<0>(0) = 1500;
	tt2.template get<0>(1) = 20080;
	tt2.template get<1>(1).add();
	tt2.template get<1>(1).template get<0>(1) = 1600;
	tt2.template get<0>(1) = 30080;
	tt2.template get<1>(1).add();
	tt2.template get<1>(1).template get<0>(2) = 1700;
	tt2.template get<1>(1).template hostToDevice<0>();

	tt2.template get<0>(2) = 40080;
	tt2.template get<1>(2).add();
	tt2.template get<1>(2).template get<0>(0) = 2500;
	tt2.template get<0>(2) = 50080;
	tt2.template get<1>(2).add();
	tt2.template get<1>(2).template get<0>(1) = 2600;
	tt2.template get<0>(2) = 60080;
	tt2.template get<1>(2).add();
	tt2.template get<1>(2).template get<0>(2) = 2700;
	tt2.template get<1>(2).template hostToDevice<0>();

	tt2.template hostToDevice<1>();
	openfpm::vector_gpu<aggregate<long int>> vg;
	vg.resize(9);

	CUDA_LAUNCH_DIM3(kernel_recursive_check,1,1,tt2.toKernel(),vg.toKernel());

	vg.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(vg.template get<0>(0),500);
	BOOST_REQUIRE_EQUAL(vg.template get<0>(1),600);
	BOOST_REQUIRE_EQUAL(vg.template get<0>(2),700);
	BOOST_REQUIRE_EQUAL(vg.template get<0>(3),1500);
	BOOST_REQUIRE_EQUAL(vg.template get<0>(4),1600);
	BOOST_REQUIRE_EQUAL(vg.template get<0>(5),1700);
	BOOST_REQUIRE_EQUAL(vg.template get<0>(6),2500);
	BOOST_REQUIRE_EQUAL(vg.template get<0>(7),2600);
	BOOST_REQUIRE_EQUAL(vg.template get<0>(8),2700);
}

BOOST_AUTO_TEST_CASE( vector_cuda_to_cpu_operator_equal )
{
	openfpm::vector_gpu<aggregate<int,int,double>> v1;
	openfpm::vector<aggregate<int,int,double>> v2;
	openfpm::vector<aggregate<int,int,double>,HeapMemory, memory_traits_inte > v3;
	openfpm::vector<aggregate<int,int,double>> v4;

	v2.resize(3000);

	for (size_t i = 0 ; i < 3000 ; i++)
	{
		v2.template get<0>(i) = i;
		v2.template get<1>(i) = i+300;
		v2.template get<2>(i) = i+6123.0;
	}

	v1 = v2;
	v3 = v2;
	v4 = v1;

	for (size_t i = 0 ; i < v2.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(v2.template get<0>(i),v1.template get<0>(i));
		BOOST_REQUIRE_EQUAL(v2.template get<0>(i),v3.template get<0>(i));
		BOOST_REQUIRE_EQUAL(v2.template get<0>(i),v4.template get<0>(i));

		BOOST_REQUIRE_EQUAL(v2.template get<1>(i),v1.template get<1>(i));
		BOOST_REQUIRE_EQUAL(v2.template get<1>(i),v3.template get<1>(i));
		BOOST_REQUIRE_EQUAL(v2.template get<1>(i),v4.template get<1>(i));

		BOOST_REQUIRE_EQUAL(v2.template get<2>(i),v1.template get<2>(i));
		BOOST_REQUIRE_EQUAL(v2.template get<2>(i),v3.template get<2>(i));
		BOOST_REQUIRE_EQUAL(v2.template get<2>(i),v4.template get<2>(i));
	}
}


BOOST_AUTO_TEST_CASE( vector_cuda_host_to_device_check )
{
	openfpm::vector_gpu<aggregate<int,int,double>> v1;

	v1.resize(3);

	for (size_t i = 0 ; i < v1.size() ; i++)
	{
		v1.template get<0>(i) = i;
		v1.template get<1>(i) = i+300;
		v1.template get<2>(i) = i+6123.0;
	}

	v1.hostToDevice<0,1,2>();

	// Now we reset the element 0, 1

	for (size_t i = 0 ; i < v1.size()-1 ; i++)
	{
		v1.template get<0>(i) = 0;
		v1.template get<1>(i) = 0;
		v1.template get<2>(i) = 0;
	}

	v1.hostToDevice<0,1,2>(v1.size()-1,v1.size()-1);

	v1.deviceToHost<0,1,2>();

	for (size_t i = 0 ; i < v1.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(v1.template get<0>(i),i);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i),i+300);
		BOOST_REQUIRE_EQUAL(v1.template get<2>(i),i+6123.0);
	}
}

BOOST_AUTO_TEST_CASE( vector_cuda_host_to_device_check_NUMA )
{
	openfpm::vector_gpu<aggregate<int,int,double>> v1;

	v1.resize(3);

	for (size_t i = 0 ; i < v1.size() ; i++)
	{
		v1.template get<0>(i) = i;
		v1.template get<1>(i) = i+300;
		v1.template get<2>(i) = i+6123.0;
	}

	v1.hostToDeviceNUMA<0,1,2>();

	// Now we reset the element 0, 1

	for (size_t i = 0 ; i < v1.size()-1 ; i++)
	{
		v1.template get<0>(i) = 0;
		v1.template get<1>(i) = 0;
		v1.template get<2>(i) = 0;
	}

	v1.hostToDeviceNUMA<0,1,2>(v1.size()-1,v1.size()-1);

	v1.deviceToHost<0,1,2>();

	for (size_t i = 0 ; i < v1.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(v1.template get<0>(i),i);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i),i+300);
		BOOST_REQUIRE_EQUAL(v1.template get<2>(i),i+6123.0);
	}
}

BOOST_AUTO_TEST_CASE( vector_cuda_host_to_device_vector_and_point_tensor )
{
	openfpm::vector_gpu<aggregate<float[3],float[3][3]>> v1;

	v1.resize(100);

	for (size_t i = 0 ; i < 50 ; i++)
	{
		v1.template get<0>(i)[0] = i+1500;
		v1.template get<0>(i)[1] = i+2200;
		v1.template get<0>(i)[2] = i+2600;

		v1.template get<1>(i)[0][0] = i+6000;
		v1.template get<1>(i)[0][1] = i+7200;
		v1.template get<1>(i)[0][2] = i+8600;
		v1.template get<1>(i)[1][0] = i+9000;
		v1.template get<1>(i)[1][1] = i+10200;
		v1.template get<1>(i)[1][2] = i+11600;
		v1.template get<1>(i)[2][0] = i+12800;
		v1.template get<1>(i)[2][1] = i+22200;
		v1.template get<1>(i)[2][2] = i+23600;
	}

	v1.hostToDevice<0,1>(0,50);

	for (size_t i = 50 ; i < 100 ; i++)
	{
		v1.template get<0>(i)[0] = i+1500;
		v1.template get<0>(i)[1] = i+2200;
		v1.template get<0>(i)[2] = i+2600;

		v1.template get<1>(i)[0][0] = i+6000;
		v1.template get<1>(i)[0][1] = i+7200;
		v1.template get<1>(i)[0][2] = i+8600;
		v1.template get<1>(i)[1][0] = i+9000;
		v1.template get<1>(i)[1][1] = i+10200;
		v1.template get<1>(i)[1][2] = i+11600;
		v1.template get<1>(i)[2][0] = i+12800;
		v1.template get<1>(i)[2][1] = i+22200;
		v1.template get<1>(i)[2][2] = i+23600;
	}

	v1.hostToDevice<0,1>(50,99);

	v1.deviceToHost<0,1>();

	for (size_t i = 0 ; i < 100 ; i++)
	{
		BOOST_REQUIRE_EQUAL(v1.template get<0>(i)[0],i+1500);
		BOOST_REQUIRE_EQUAL(v1.template get<0>(i)[1],i+2200);
		BOOST_REQUIRE_EQUAL(v1.template get<0>(i)[2],i+2600);

		BOOST_REQUIRE_EQUAL(v1.template get<1>(i)[0][0],i+6000);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i)[0][1],i+7200);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i)[0][2],i+8600);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i)[1][0],i+9000);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i)[1][1],i+10200);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i)[1][2],i+11600);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i)[2][0],i+12800);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i)[2][1],i+22200);
		BOOST_REQUIRE_EQUAL(v1.template get<1>(i)[2][2],i+23600);
	}
}

BOOST_AUTO_TEST_CASE( vector_cuda_copy )
{
	openfpm::vector_gpu<aggregate<float,float[3],float[3][3]>> v1;
	openfpm::vector_gpu<aggregate<float,float[3],float[3][3]>> v2;

	v1.resize(100);

	auto ite = v1.getIterator();

	while (ite.isNext())
	{
		auto p = ite.get();

		v1.template get<0>(p) = p + 100;

		v1.template get<0>(p) = p + 2000;
		v1.template get<0>(p) = p + 3000;
		v1.template get<0>(p) = p + 4000;

		v1.template get<1>(p)[0] = p + 5000;
		v1.template get<1>(p)[1] = p + 6000;
		v1.template get<1>(p)[2] = p + 7000;

		v1.template get<2>(p)[0][0] = p + 8000;
		v1.template get<2>(p)[0][1] = p + 9000;
		v1.template get<2>(p)[0][2] = p + 10000;

		v1.template get<2>(p)[1][0] = p + 11000;
		v1.template get<2>(p)[1][1] = p + 12000;
		v1.template get<2>(p)[2][2] = p + 13000;

		v1.template get<2>(p)[2][0] = p + 14000;
		v1.template get<2>(p)[2][1] = p + 15000;
		v1.template get<2>(p)[2][2] = p + 16000;

		++ite;
	}

	v1.hostToDevice<0,1,2>();

	ite = v1.getIterator();

	while (ite.isNext())
	{
		auto p = ite.get();

		v1.template get<0>(p) = p + 6100;

		v1.template get<0>(p) = p + 62000;
		v1.template get<0>(p) = p + 63000;
		v1.template get<0>(p) = p + 64000;

		v1.template get<1>(p)[0] = p + 65000;
		v1.template get<1>(p)[1] = p + 66000;
		v1.template get<1>(p)[2] = p + 67000;

		v1.template get<2>(p)[0][0] = p + 68000;
		v1.template get<2>(p)[0][1] = p + 69000;
		v1.template get<2>(p)[0][2] = p + 610000;

		v1.template get<2>(p)[1][0] = p + 611000;
		v1.template get<2>(p)[1][1] = p + 612000;
		v1.template get<2>(p)[2][2] = p + 613000;

		v1.template get<2>(p)[2][0] = p + 614000;
		v1.template get<2>(p)[2][1] = p + 615000;
		v1.template get<2>(p)[2][2] = p + 616000;

		++ite;
	}

	v2 = v1;

	// first check the CPU

	bool match = true;

	ite = v2.getIterator();

	while (ite.isNext())
	{
		auto p = ite.get();

		match = v2.template get<0>(p) == p + 6100;

		match = v2.template get<0>(p) == p + 62000;
		match = v2.template get<0>(p) == p + 63000;
		match = v2.template get<0>(p) == p + 64000;

		match = v2.template get<1>(p)[0] == p + 65000;
		match = v2.template get<1>(p)[1] == p + 66000;
		match = v2.template get<1>(p)[2] == p + 67000;

		match = v2.template get<2>(p)[0][0] == p + 68000;
		match = v2.template get<2>(p)[0][1] == p + 69000;
		match = v2.template get<2>(p)[0][2] == p + 610000;

		match = v2.template get<2>(p)[1][0] == p + 611000;
		match = v2.template get<2>(p)[1][1] == p + 612000;
		match = v2.template get<2>(p)[2][2] == p + 613000;

		match = v2.template get<2>(p)[2][0] == p + 614000;
		match = v2.template get<2>(p)[2][1] == p + 615000;
		match = v2.template get<2>(p)[2][2] == p + 616000;

		++ite;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	v2.deviceToHost<0,1,2>();

	ite = v2.getIterator();

	while (ite.isNext())
	{
		auto p = ite.get();

		match = v2.template get<0>(p) == p + 100;

		match = v2.template get<0>(p) == p + 2000;
		match = v2.template get<0>(p) == p + 3000;
		match = v2.template get<0>(p) == p + 4000;

		match = v2.template get<1>(p)[0] == p + 5000;
		match = v2.template get<1>(p)[1] == p + 6000;
		match = v2.template get<1>(p)[2] == p + 7000;

		match = v2.template get<2>(p)[0][0] == p + 8000;
		match = v2.template get<2>(p)[0][1] == p + 9000;
		match = v2.template get<2>(p)[0][2] == p + 10000;

		match = v2.template get<2>(p)[1][0] == p + 11000;
		match = v2.template get<2>(p)[1][1] == p + 12000;
		match = v2.template get<2>(p)[2][2] == p + 13000;

		match = v2.template get<2>(p)[2][0] == p + 14000;
		match = v2.template get<2>(p)[2][1] == p + 15000;
		match = v2.template get<2>(p)[2][2] == p + 16000;

		if (match == false)
		{
			std::cout << v2.template get<0>(p) << std::endl;
		}

		++ite;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_SUITE_END()

