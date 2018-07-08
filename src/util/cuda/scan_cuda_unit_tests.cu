#define BOOST_GPU_ENABLED __host__ __device__

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "util/cuda_util.hpp"
#include "util/boost/boost_array_openfpm.hpp"
#include "Vector/map_vector.hpp"
#include "scan_cuda.cuh"

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

	compress4<cnt_type,ids_type><<<ite.wthr,ite.thr>>>(cl_n.size(),
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

	// fill with some data

	for (size_t i = 0 ; i < cl_n.size() ; i++)
	{cl_n.template get<0>(i) = i%16;}

	int nblocks = ((cl_n.size() / (ratio) ) + THREADS - 1 ) / THREADS;

	red.resize(nblocks);

	cl_n.template hostToDevice<0>();

	breduce<THREADS/32,cnt_type,ids_type,ratio_reduction<cnt_type,ids_type>><<<nblocks,THREADS>>>(cl_n.size()/ratio,
														  static_cast<cnt_type *>(cl_n.template getDeviceBuffer<0>()),
														  static_cast<cnt_type *>(red.template getDeviceBuffer<0>()));

	red.template deviceToHost<0>();

	for (size_t i = 0 ; i < red.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(red.template get<0>(i),120*128);
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

	bexscan<THREADS,cnt_type><<<1,THREADS,nblocks*sizeof(unsigned int)>>>(nblocks,
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
	size_t nb = 4;

	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> base;
	openfpm::vector<aggregate<ids_type>,CudaMemory,typename memory_traits_inte<aggregate<ids_type>>::type,memory_traits_inte> cl_n;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> cl_n_scan;

	base.resize(nb);
	cl_n.resize(2048*nb);
	cl_n_scan.resize(2048*nb);

	constexpr int THREADS = 128;

	// fill with some data

	for (size_t i = 0 ; i < base.size() ; i++)
	{base.template get<0>(i) = (i+1)*120*128;}

	for (size_t i = 0 ; i < cl_n.size() ; i++)
	{cl_n.template get<0>(i) = i%16;}

	int nblocks = cl_n.size() / 16;

	cl_n.template hostToDevice<0>();
	base.template hostToDevice<0>();

	gexscan<THREADS/32,ratio_extend<unsigned int,unsigned char>> <<< cl_n.size() / 16 / THREADS, THREADS >>>(nblocks,
																									  static_cast<ratio_extend<unsigned int,unsigned char>::cnt_type4 *>(cl_n.template getDeviceBuffer<0>()),
																									  static_cast<cnt_type *>(base.template getDeviceBuffer<0>()),
												                                                      static_cast<ratio_extend<unsigned int,unsigned char>::cnt_type4 *>(cl_n_scan.template getDeviceBuffer<0>()));

	cl_n_scan.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(0,cl_n_scan.template get<0>(0));

	size_t scan = 0;
	for (size_t i = 1 ; i < cl_n_scan.size() ; i++)
	{
		scan += cl_n.template get<0>(i-1);
		BOOST_REQUIRE_EQUAL(cl_n_scan.template get<0>(i),scan);
	}
}

BOOST_AUTO_TEST_CASE( test_base_funcs )
{
	std::cout << "Test cell list GPU base func" << "\n";

	test_compress<int,unsigned char>();

	test_breduce<unsigned int, unsigned char>();

	test_bexscan<unsigned int>();

	test_gexscan<unsigned int, unsigned char>();

	std::cout << "End cell list GPU" << "\n";

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
	{cl_n.template get<0>(i) = i%16;}

	cl_n.template hostToDevice<0>();

	scan<cnt_type,ids_type>(cl_n,cl_n_scan);

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
	std::cout << "Test cell list GPU scan" << "\n";

	test_scan<unsigned int, unsigned char>(8192);

	test_scan<unsigned int, unsigned char>(25);

	test_scan<unsigned int, unsigned char>(139);

	test_scan<unsigned int, unsigned char>(1025);

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_SUITE_END()
