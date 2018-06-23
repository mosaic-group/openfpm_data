/*
 * CellList_gpu_test.cpp
 *
 *  Created on: Jun 13, 2018
 *      Author: i-bird
 */

#define BOOST_GPU_ENABLED __host__ __device__

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "util/cuda_util.hpp"
#include "CellList_gpu.hpp"

BOOST_AUTO_TEST_SUITE( CellList_gpu_test )

template<unsigned int dim, typename T, typename cnt_type, typename ids_type>
void test_sub_index()
{
	openfpm::vector<aggregate<cnt_type,cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type,cnt_type>>::type,memory_traits_inte> cl_n;
	openfpm::vector<aggregate<ids_type[dim+1]>,CudaMemory,typename memory_traits_inte<aggregate<ids_type[dim+1]>>::type,memory_traits_inte> part_ids;

	// fill with some particles

	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> pl;

	// create 3 particles

	Point<dim,T> p1({0.2,0.2,0.2});
	Point<dim,T> p2({0.9,0.2,0.2});
	Point<dim,T> p3({0.2,0.9,0.2});
	Point<dim,T> p4({0.2,0.2,0.9});
	Point<dim,T> p5({0.9,0.9,0.2});
	Point<dim,T> p6({0.9,0.2,0.9});
	Point<dim,T> p7({0.2,0.9,0.9});
	Point<dim,T> p8({0.9,0.9,0.9});
	Point<dim,T> p9({0.0,0.0,0.0});
	Point<dim,T> p10({0.205,0.205,0.205});

	pl.add(p1);
	pl.add(p2);
	pl.add(p3);
	pl.add(p4);
	pl.add(p5);
	pl.add(p6);
	pl.add(p7);
	pl.add(p8);
	pl.add(p9);
	pl.add(p10);

	CudaMemory spacing;
	CudaMemory div;

	spacing.allocate(sizeof(T)*dim);
	div.allocate(dim*sizeof(ids_type));

	ids_type (& div_p)[dim] = *static_cast<ids_type (*)[dim]>(div.getPointer());
	T (& spacing_p)[dim] = *static_cast<T (*)[dim]>(spacing.getPointer());

	for (size_t i = 0 ; i < dim ; i++)
	{
		div_p[i] = 17;
		spacing_p[i] = 0.1;
	}

	cl_n.resize(17*17*17);
	CUDA_SAFE(cudaMemset(cl_n.template getDeviceBuffer<0>(),0,cl_n.size()*sizeof(cnt_type)));

	part_ids.resize(9);

	size_t sz[3] = {17,17,17};
	grid_sm<3,void> gr(sz);

	// Force to copy into device
	spacing.getDevicePointer();
	div.getDevicePointer();

	auto ite = pl.getGPUIterator();

	subindex<dim,T,cnt_type,ids_type><<<ite.wthr,ite.thr>>>(*static_cast<ids_type (*)[dim]>(div.getDevicePointer()),
																	*static_cast<T (*)[dim]>(spacing.getDevicePointer()),
																	pl.capacity(),
																	pl.size(),
																	static_cast<T *>(pl.template getDeviceBuffer<0>()),
																	static_cast<cnt_type *>(cl_n.template getDeviceBuffer<0>()),
																	static_cast<ids_type *>(part_ids.template getDeviceBuffer<0>()));

	cl_n.template deviceToHost<0>();
	part_ids.template deviceToHost<0>();

	// I have to find != 0 in the cell with particles

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(0)[0],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(0)[1],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(0)[2],2);

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(1)[0],9);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(1)[1],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(1)[2],2);

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(2)[0],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(2)[1],9);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(2)[2],2);

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(3)[0],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(3)[1],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(3)[2],9);

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(4)[0],9);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(4)[1],9);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(4)[2],2);

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(5)[0],9);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(5)[1],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(5)[2],9);

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(6)[0],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(6)[1],9);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(6)[2],9);

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(7)[0],9);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(7)[1],9);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(7)[2],9);

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(8)[0],0);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(8)[1],0);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(8)[2],0);

	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(9)[0],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(9)[1],2);
	BOOST_REQUIRE_EQUAL(part_ids.template get<0>(9)[2],2);

	BOOST_REQUIRE_EQUAL(cl_n.template get<0>(gr.LinId({2,2,2})),2);
	BOOST_REQUIRE_EQUAL(cl_n.template get<0>(gr.LinId({9,2,2})),1);
	BOOST_REQUIRE_EQUAL(cl_n.template get<0>(gr.LinId({2,9,2})),1);
	BOOST_REQUIRE_EQUAL(cl_n.template get<0>(gr.LinId({2,2,9})),1);
	BOOST_REQUIRE_EQUAL(cl_n.template get<0>(gr.LinId({9,9,2})),1);
	BOOST_REQUIRE_EQUAL(cl_n.template get<0>(gr.LinId({9,2,9})),1);
	BOOST_REQUIRE_EQUAL(cl_n.template get<0>(gr.LinId({2,9,9})),1);
	BOOST_REQUIRE_EQUAL(cl_n.template get<0>(gr.LinId({9,9,9})),1);
	BOOST_REQUIRE_EQUAL(cl_n.template get<0>(gr.LinId({0,0,0})),1);
}

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
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> base;
	openfpm::vector<aggregate<ids_type>,CudaMemory,typename memory_traits_inte<aggregate<ids_type>>::type,memory_traits_inte> cl_n;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> cl_n_scan;

	base.resize(500);
	cl_n.resize(8192);

	constexpr int THREADS = 128;

	// fill with some data

	for (size_t i = 0 ; i < base.size() ; i++)
	{base.template get<0>(i) = i*120*128;}

	int nblocks = base.size();

	gexscan<THREADS/32,ratio_extend<unsigned int,unsigned char>> <<< base.size() / THREADS, THREADS >>>(nblocks,
																									  static_cast<ratio_extend<unsigned int,unsigned char>::cnt_type4 *>(cl_n.template getDeviceBuffer<0>()),
																									  static_cast<cnt_type *>(base.template getDeviceBuffer<0>()),
												                                                      static_cast<ratio_extend<unsigned int,unsigned char>::cnt_type4 *>(cl_n_scan.template getDeviceBuffer<0>()));

	base.template deviceToHost<0>();

	for (size_t i = 0 ; i < base.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(base.template get<0>(i),i+1);
	}
}

BOOST_AUTO_TEST_CASE( test_base_funcs )
{
	std::cout << "Test cell list GPU base func" << "\n";

	test_sub_index<3,float,int,unsigned char>();

	test_compress<int,unsigned char>();

	test_breduce<unsigned int, unsigned char>();

	test_bexscan<unsigned int>();

	test_gexscan<unsigned int, unsigned char>();

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

template<unsigned int dim, typename T, typename CellS> void Test_cell_gpu(SpaceBox<dim,T> & box)
{
	//Space where is living the Cell list
	//SpaceBox<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Subdivisions
	size_t div[dim] = {16,16,16};

	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cl2(box,div);

	// vector of particles

	openfpm::vector<Point<dim,double>,CudaMemory,typename memory_traits_inte<Point<dim,double>>::type,memory_traits_inte> pl;

	// create 3 particles

	Point<dim,double> p1({0.2,0.2,0.2});
	Point<dim,double> p2({0.9,0.2,0.2});
	Point<dim,double> p3({0.2,0.9,0.2});
	Point<dim,double> p4({0.2,0.2,0.9});
	Point<dim,double> p5({0.9,0.9,0.2});
	Point<dim,double> p6({0.9,0.2,0.9});
	Point<dim,double> p7({0.2,0.9,0.9});
	Point<dim,double> p8({0.9,0.9,0.9});
	Point<dim,double> p9({0.0,0.0,0.0});

	pl.add(p1);
	pl.add(p2);
	pl.add(p3);
	pl.add(p4);
	pl.add(p5);
	pl.add(p6);
	pl.add(p7);
	pl.add(p8);
	pl.add(p9);

	cl2.construct(pl);
}


BOOST_AUTO_TEST_CASE( CellList_gpu_use)
{
	std::cout << "Test cell list GPU" << "\n";

	SpaceBox<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	SpaceBox<3,double> box2({-1.0f,-1.0f,-1.0f},{1.0f,1.0f,1.0f});

	Test_cell_gpu<3,double,CellList_gpu<3,double,CudaMemory>>(box);
	Test_cell_gpu<3,double,CellList_gpu<3,double,CudaMemory>>(box2);

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_SUITE_END()

