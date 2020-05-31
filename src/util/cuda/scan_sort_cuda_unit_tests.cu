#define BOOST_GPU_ENABLED __host__ __device__


#include <hip/hip_runtime.h>
#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "util/cuda_util.hpp"
#include "util/boost/boost_array_openfpm.hpp"
#include "Vector/map_vector.hpp"
#include "scan_cuda.cuh"

#define SCAN_WITH_CUB
#define SORT_WITH_CUB

#include "sort_ofp.cuh"
#include "scan_ofp.cuh"

BOOST_AUTO_TEST_SUITE( scan_tests )

template<typename cnt_type, typename ids_type>
void test_compress()
{
	openfpm::vector<aggregate<cnt_type,cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type,cnt_type>>::type,memory_traits_inte> cl_n;

	openfpm::vector<aggregate<ids_type>,CudaMemory,typename memory_traits_inte<aggregate<ids_type>>::type,memory_traits_inte> compressed;

	// fill some counting

	cl_n.resize(12);

	cl_n.template get<0>(0) = 3;
	cl_n.template get<0>(1) = 5;
	cl_n.template get<0>(2) = 8;
	cl_n.template get<0>(3) = 1;
	cl_n.template get<0>(4) = 0;
	cl_n.template get<0>(5) = 0;
	cl_n.template get<0>(6) = 21;
	cl_n.template get<0>(7) = 4;
	cl_n.template get<0>(8) = 4;
	cl_n.template get<0>(9) = 6;
	cl_n.template get<0>(10) = 10;

	compressed.resize(cl_n.size());

	auto ite = cl_n.getGPUIterator();
	ite.thr.x /= 4;

	cl_n.template hostToDevice<0>();

	hipLaunchKernelGGL(HIP_KERNEL_NAME(compress4<cnt_type,ids_type>), dim3(ite.wthr), dim3(ite.thr), 0, 0, cl_n.size(),
														  static_cast<cnt_type *>(cl_n.template getDeviceBuffer<0>()),
														  static_cast<ids_type *>(compressed.template getDeviceBuffer<0>()));

	compressed.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(compressed.template get<0>(0),3);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(1),5);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(2),8);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(3),1);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(4),0);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(5),0);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(6),21);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(7),4);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(8),4);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(9),6);
	BOOST_REQUIRE_EQUAL(compressed.template get<0>(10),10);
}


template<typename cnt_type, typename ids_type>
void test_breduce()
{
	openfpm::vector<aggregate<ids_type>,CudaMemory,typename memory_traits_inte<aggregate<ids_type>>::type,memory_traits_inte> cl_n;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> red;

	cl_n.resize(8192);

	constexpr int THREADS = 128;
	constexpr int ratio = 4*sizeof(cnt_type)/sizeof(ids_type);
	constexpr int red_tot = THREADS * ratio;

	int nblocks = ((cl_n.size() / (ratio) ) + THREADS - 1 ) / THREADS;

	// fill with some data

	openfpm::vector<cnt_type> block_red;
	block_red.resize(nblocks);

	for (size_t i = 0 ; i < cl_n.size() ; i++)
	{
		if ((i % red_tot)/256  == 0)
		{
			cl_n.template get<0>(i) = i%128;
			block_red.get(i/red_tot) += cl_n.template get<0>(i);
		}
		else if ((i % red_tot)/256  == 1)
		{
			cl_n.template get<0>(i) = i%7;
			block_red.get(i/red_tot) += cl_n.template get<0>(i);
		}
		else if ((i % red_tot)/256  == 2)
		{
			cl_n.template get<0>(i) = i%13;
			block_red.get(i/red_tot) += cl_n.template get<0>(i);
		}
		else if ((i % red_tot)/256  == 3)
		{
			cl_n.template get<0>(i) = i%17;
			block_red.get(i/red_tot) += cl_n.template get<0>(i);
		}
		else if ((i % red_tot)/256  == 4)
		{
			cl_n.template get<0>(i) = i%128;
			block_red.get(i/red_tot) += cl_n.template get<0>(i);
		}
		else if ((i % red_tot)/256  == 5)
		{
			cl_n.template get<0>(i) = i%7;
			block_red.get(i/red_tot) += cl_n.template get<0>(i);
		}
		else if ((i % red_tot)/256  == 6)
		{
			cl_n.template get<0>(i) = i%13;
			block_red.get(i/red_tot) += cl_n.template get<0>(i);
		}
		else if ((i % red_tot)/256  == 7)
		{
			cl_n.template get<0>(i) = i%17;
			block_red.get(i/red_tot) += cl_n.template get<0>(i);
		}
	}

	red.resize(nblocks);

	cl_n.template hostToDevice<0>();

	hipLaunchKernelGGL(HIP_KERNEL_NAME(breduce<THREADS/32,cnt_type,ids_type,ratio_reduction<cnt_type,ids_type>>), dim3(nblocks), dim3(THREADS), 0, 0, cl_n.size()/ratio*4,
														  static_cast<cnt_type *>(cl_n.template getDeviceBuffer<0>()),
														  static_cast<cnt_type *>(red.template getDeviceBuffer<0>()));

	red.template deviceToHost<0>();

	for (size_t i = 0 ; i < red.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(red.template get<0>(i),block_red.get(i));
	}
}

template<typename cnt_type>
void test_bexscan()
{
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> base;

	base.resize(500);

	constexpr int THREADS = 128;

	// fill with some data

	for (size_t i = 0 ; i < base.size() ; i++)
	{base.template get<0>(i) = 1;}

	int nblocks = base.size();

	base.template hostToDevice<0>();

	hipLaunchKernelGGL(HIP_KERNEL_NAME(bexscan<THREADS,cnt_type>), dim3(1), dim3(THREADS), nblocks*sizeof(unsigned int), 0, nblocks,
														  	  	  	  	  static_cast<cnt_type *>(base.template getDeviceBuffer<0>()));

	base.template deviceToHost<0>();

	for (size_t i = 0 ; i < base.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(base.template get<0>(i),i+1);
	}
}

template<typename cnt_type, typename ids_type>
void test_gexscan()
{
	size_t nb = 16;

	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> base;
	openfpm::vector<aggregate<ids_type>,CudaMemory,typename memory_traits_inte<aggregate<ids_type>>::type,memory_traits_inte> cl_n;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> cl_n_scan;

	constexpr int ratio = sizeof(cnt_type)/sizeof(ids_type);
	constexpr int THREADS = 128;
	constexpr int expand = THREADS * ratio;

	base.resize(nb);
	cl_n.resize(expand*nb);
	cl_n_scan.resize(expand*nb);

	// fill with some data

	for (size_t i = 0 ; i < base.size() ; i++)
	{base.template get<0>(i) = (i+1)*120*THREADS/4*ratio;}

	for (size_t i = 0 ; i < cl_n.size() ; i++)
	{cl_n.template get<0>(i) = i%16;}

	int nblocks = cl_n.size() / ratio;

	cl_n.template hostToDevice<0>();
	base.template hostToDevice<0>();

	hipLaunchKernelGGL(HIP_KERNEL_NAME(gexscan<THREADS/32,ratio_extend<cnt_type,ids_type>>), dim3(cl_n.size() / 4 / ratio / THREADS), dim3(THREADS), 0, 0, nblocks,
																									  static_cast<typename ratio_extend<cnt_type,ids_type>::cnt_type4 *>(cl_n.template getDeviceBuffer<0>()),
																									  static_cast<cnt_type *>(base.template getDeviceBuffer<0>()),
												                                                      static_cast<typename ratio_extend<cnt_type,ids_type>::cnt_type4 *>(cl_n_scan.template getDeviceBuffer<0>()));

	cl_n_scan.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(0,cl_n_scan.template get<0>(0));

	size_t scan = 0;
	for (size_t i = 1 ; i < cl_n_scan.size() ; i++)
	{
		scan += cl_n.template get<0>(i-1);
		BOOST_REQUIRE_EQUAL(cl_n_scan.template get<0>(i),scan);
	}
}

BOOST_AUTO_TEST_CASE (test_breduce_func )
{
	test_breduce<unsigned int, unsigned char>();

	test_breduce<unsigned int, unsigned short>();

	test_breduce<unsigned int, unsigned int>();

	test_breduce<unsigned int, char>();

	test_breduce<unsigned int, short>();

	test_breduce<unsigned int, int>();

	test_breduce<int, unsigned char>();

	test_breduce<int, unsigned short>();

	test_breduce<int, unsigned int>();

	test_breduce<int, char>();

	test_breduce<int, short>();

	test_breduce<int, int>();
}

BOOST_AUTO_TEST_CASE(test_compress_functions)
{
	test_compress<unsigned int, unsigned char>();

	test_compress<unsigned int, unsigned short>();

	test_compress<unsigned int, unsigned int>();

	test_compress<unsigned int, char>();

	test_compress<unsigned int, short>();

	test_compress<unsigned int, int>();

	test_compress<int, unsigned char>();

	test_compress<int, unsigned short>();

	test_compress<int, unsigned int>();

	test_compress<int, char>();

	test_compress<int, short>();

	test_compress<int, int>();
}

BOOST_AUTO_TEST_CASE(test_bexscan_functions)
{
	test_bexscan<unsigned int>();

	test_bexscan<int>();
}

BOOST_AUTO_TEST_CASE( test_gexscan_funcs )
{
	std::cout << "Test gexscan GPU base func" << "\n";

	test_gexscan<unsigned int, unsigned char>();

	test_gexscan<unsigned int, unsigned short>();

	test_gexscan<unsigned int, unsigned int>();

	test_gexscan<int, unsigned char>();

	test_gexscan<int, unsigned short>();

	test_gexscan<int, unsigned int>();

	std::cout << "End gexscan GPU" << "\n";

	// Test the cell list
}

template<typename cnt_type, typename ids_type>
void test_scan(size_t num)
{
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> cl_n;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> cl_n_scan;

	cl_n.resize(num);

	// fill with some data

	for (size_t i = 0 ; i < cl_n.size() ; i++)
	{cl_n.template get<0>(i) = 255.0*rand()/RAND_MAX;}

	cl_n.template hostToDevice<0>();

	scan<cnt_type,ids_type>sc;
	sc.scan_(cl_n,cl_n_scan);

	cl_n_scan.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(0,cl_n_scan.template get<0>(0));

	size_t scan = 0;
	for (size_t i = 1 ; i < cl_n_scan.size() ; i++)
	{
		scan += cl_n.template get<0>(i-1);
		BOOST_REQUIRE_EQUAL(cl_n_scan.template get<0>(i),scan);
	}
}

BOOST_AUTO_TEST_CASE( test_scan_algo )
{
	std::cout << "Test GPU obsolete scan" << "\n";

	test_scan<unsigned int, unsigned char>(8192);

	test_scan<unsigned int, unsigned char>(25);

	test_scan<unsigned int, unsigned char>(139);

	test_scan<unsigned int, unsigned char>(1025);

	test_scan<unsigned int, unsigned short>(8192);

	test_scan<unsigned int, unsigned short>(25);

	test_scan<unsigned int, unsigned short>(139);

	test_scan<unsigned int, unsigned short>(1025);

	test_scan<unsigned int, unsigned int>(8192);

	test_scan<unsigned int, unsigned int>(25);

	test_scan<unsigned int, unsigned int>(139);

	test_scan<unsigned int, unsigned int>(1025);

	std::cout << "End GPU obsolete scan" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( test_scan_cub_wrapper )
{
	std::cout << "Test scan CUB" << "\n";

	openfpm::vector_gpu<aggregate<unsigned int>> input;
	openfpm::vector_gpu<aggregate<unsigned int>> output;

	openfpm::vector_gpu<aggregate<unsigned char>> temporal;

	input.resize(10000);
	output.resize(10000);

	// fill input

	for (size_t i = 0 ; i < 10000; i++)
	{
		input.template get<0>(i) = 10.0*(float)rand() / RAND_MAX;
	}

	input.template hostToDevice<0>();

	mgpu::ofp_context_t context;
	openfpm::scan((unsigned int *)input.template getDeviceBuffer<0>(),input.size(),(unsigned int *)output.template getDeviceBuffer<0>(),context);

    output.template deviceToHost<0>();

    size_t cnt = 0;
    for (size_t i = 0 ; i < input.size() ; i++)
    {
    	BOOST_REQUIRE_EQUAL(cnt,output.template get<0>(i));
    	cnt += input.template get<0>(i);
    }

	std::cout << "End scan CUB" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( test_sort_cub_wrapper )
{
	std::cout << "Test sort CUB" << "\n";

	openfpm::vector_gpu<aggregate<unsigned int>> input;
	openfpm::vector_gpu<aggregate<unsigned int>> input_id;

	openfpm::vector_gpu<aggregate<unsigned char>> temporal;

	input.resize(10000);
	input_id.resize(10000);


	// fill input

	for (size_t i = 0 ; i < 10000; i++)
	{
		input.template get<0>(i) = 10000.0*(float)rand() / RAND_MAX;
		input_id.template get<0>(i) = i;
	}

	input.template hostToDevice<0>();
	input_id.template hostToDevice<0>();

	mgpu::ofp_context_t context;

	openfpm::sort((unsigned int *)input.template getDeviceBuffer<0>(),
				  (unsigned int *)input_id.template getDeviceBuffer<0>(),
			      input.size(),mgpu::template less_t<unsigned int>(),context);

	input.template deviceToHost<0>();
	input_id.template deviceToHost<0>();

    for (size_t i = 0 ; i < input.size() - 1 ; i++)
    {
    	BOOST_REQUIRE(input.template get<0>(i) <= input.template get<0>(i+1));
    }

	openfpm::sort((unsigned int *)input.template getDeviceBuffer<0>(),
				  (unsigned int *)input_id.template getDeviceBuffer<0>(),
			      input.size(),mgpu::template greater_t<unsigned int>(),context);

	input.template deviceToHost<0>();
	input_id.template deviceToHost<0>();

    for (size_t i = 0 ; i < input.size() - 1 ; i++)
    {
    	BOOST_REQUIRE(input.template get<0>(i) >= input.template get<0>(i+1));
    }

	std::cout << "End sort CUB" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_SUITE_END()
