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
#include "cuda/CellList_gpu.hpp"
#include "CellList.hpp"
#include "util/boost/boost_array_openfpm.hpp"
#include  "Point_test.hpp"
#include "util/cuda/moderngpu/kernel_load_balance.hxx"
#include "util/cuda/scan_cuda.cuh"

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

	openfpm::array<T,dim,cnt_type> spacing;
	openfpm::array<ids_type,dim,cnt_type> div;
	openfpm::array<ids_type,dim,cnt_type> off;

	for (size_t i = 0 ; i < dim ; i++)
	{
		div[i] = 17;
		spacing[i] = 0.1;
		off[i] = 0;
	}

	cl_n.resize(17*17*17);
	CUDA_SAFE(cudaMemset(cl_n.template getDeviceBuffer<0>(),0,cl_n.size()*sizeof(cnt_type)));

	part_ids.resize(pl.size());

	size_t sz[3] = {17,17,17};
	grid_sm<3,void> gr(sz);

	auto ite = pl.getGPUIterator();

	pl.template hostToDevice<0>();

	Matrix<dim,T> mt;
	Point<dim,T> pt;

	no_transform_only<dim,T> t(mt,pt);

	subindex<dim,T,cnt_type,ids_type,no_transform_only<dim,T>><<<ite.wthr,ite.thr>>>(div,
																	spacing,
																	off,
																	t,
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


template<unsigned int dim, typename T, typename cnt_type, typename ids_type>
void test_sub_index2()
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

	Point<dim,T> pt({-0.3,-0.3,-0.3});

	p1 += pt;
	p2 += pt;
	p3 += pt;
	p4 += pt;
	p5 += pt;
	p6 += pt;
	p7 += pt;
	p8 += pt;
	p9 += pt;
	p10 += pt;


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

	openfpm::array<T,dim,cnt_type> spacing;
	openfpm::array<ids_type,dim,cnt_type> div;
	openfpm::array<ids_type,dim,cnt_type> off;

	for (size_t i = 0 ; i < dim ; i++)
	{
		div[i] = 17;
		spacing[i] = 0.1;
		off[i] = 0;
	}

	cl_n.resize(17*17*17);
	CUDA_SAFE(cudaMemset(cl_n.template getDeviceBuffer<0>(),0,cl_n.size()*sizeof(cnt_type)));

	part_ids.resize(pl.size());

	size_t sz[3] = {17,17,17};
	grid_sm<3,void> gr(sz);

	auto ite = pl.getGPUIterator();

	pl.template hostToDevice<0>();

	Matrix<dim,T> mt;

	shift_only<dim,T> t(mt,pt);

	subindex<dim,T,cnt_type,ids_type,shift_only<dim,T>><<<ite.wthr,ite.thr>>>(div,
																	spacing,
																	off,
																	t,
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

template<unsigned int dim, typename T>
void create_n_part(int n_part,
		           openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> & pl,
		           CellList<dim,T, Mem_fast<>> & cl)
{
	pl.resize(n_part);

	auto it = pl.getIterator();

	while(it.isNext())
	{
		auto p = it.get();

		pl.template get<0>(p)[0] = (double)rand()/RAND_MAX;
		pl.template get<0>(p)[1] = (double)rand()/RAND_MAX;
		pl.template get<0>(p)[2] = (double)rand()/RAND_MAX;

		Point<dim,T> xp;
		xp.get(0) = pl.template get<0>(p)[0];
		xp.get(1) = pl.template get<0>(p)[1];
		xp.get(2) = pl.template get<0>(p)[2];

		size_t c = cl.getCell(xp);
		cl.addCell(c,p);

		++it;
	}
}

template<unsigned int dim, typename T, typename cnt_type, typename ids_type>
void create_starts_and_parts_ids(CellList<dim,T, Mem_fast<>> & cl,
								 grid_sm<dim,void> & gr,
								 size_t n_part,
								 size_t n_cell,
								 openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> & starts,
								 openfpm::vector<aggregate<ids_type[dim+1]>,CudaMemory,typename memory_traits_inte<aggregate<ids_type[dim+1]>>::type,memory_traits_inte> & part_ids,
								 openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> & cells)
{
	// Construct starts and part_ids

	part_ids.resize(n_part);
	starts.resize(n_cell);
	cells.resize(n_part);

	grid_key_dx_iterator<dim> itg(gr);

	size_t start = 0;

	while (itg.isNext())
	{
		auto cell = itg.get();

		size_t clin = gr.LinId(cell);

		for (size_t j = 0 ; j < cl.getNelements(clin) ; j++)
		{
			size_t p_id = cl.get(clin,j);

			for (size_t k = 0 ; k < dim ; k++)
			{part_ids.template get<0>(p_id)[k] = cell.get(k);}

			part_ids.template get<0>(p_id)[dim] = j;

			cells.template get<0>(start+j) = p_id;
		}

		starts.template get<0>(clin) = start;
		start += cl.getNelements(clin);

		++itg;
	}
}

template<unsigned int dim, typename T, typename cnt_type, typename ids_type>
void test_fill_cell()
{
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> cells;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> cells_out;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> starts;
	openfpm::vector<aggregate<ids_type[dim+1]>,CudaMemory,typename memory_traits_inte<aggregate<ids_type[dim+1]>>::type,memory_traits_inte> part_ids;


	// CellList to check the result

	Box<dim,T> domain;

	typename openfpm::array<T,dim> spacing;
	typename openfpm::array<ids_type,dim,cnt_type> div_c;
	typename openfpm::array<ids_type,dim,cnt_type> off;

	size_t div_host[dim];

	size_t tot = 1;
	for (size_t i = 0 ; i < dim ; i++)
	{
		div_host[i] = 10;
		div_c[i] = 10;
		tot *= div_host[i];
		spacing[i] = 0.1;
		domain.setLow(i,0.0);
		domain.setHigh(i,1.0);
		off[i] = 0;
	}

	CellList<dim,T, Mem_fast<>> cl(domain,div_host,0);
	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> pl;

	create_n_part(5000,pl,cl);

	grid_sm<dim,void> gr(div_host);

	create_starts_and_parts_ids(cl,gr,pl.size(),tot,starts,part_ids,cells_out);

	bool check = true;
	cells.resize(pl.size());
	for (size_t i = 0 ; i < gr.size() - 1 ; i++)
	{
		size_t tot_p = starts.template get<0>(i+1) - starts.template get<0>(i);

		check &= (tot_p == cl.getNelements(i));

		grid_key_dx<dim> key = gr.InvLinId(i);

		// Now we check the ids

		for (size_t j = 0 ; j < cl.getNelements(i) ; j++)
		{
			size_t p_id = cl.get(i,j);

			for (size_t k = 0 ; k < dim ; k++)
			{check &= part_ids.template get<0>(p_id)[k] == key.get(k);}
		}
	}

	BOOST_REQUIRE(check == true);

	auto itgg = part_ids.getGPUIterator();

	starts.template hostToDevice<0>();
	part_ids.template hostToDevice<0>();

	// Here we test fill cell
	fill_cells<dim,cnt_type,ids_type,shift_ph<0,cnt_type>><<<itgg.wthr,itgg.thr>>>(0,
																				   div_c,
																				   off,
																				   part_ids.size(),
																				   part_ids.capacity(),
																				   static_cast<cnt_type *>(starts.template getDeviceBuffer<0>()),
																				   static_cast<ids_type *>(part_ids.template getDeviceBuffer<0>()),
																				   static_cast<cnt_type *>(cells.template getDeviceBuffer<0>()) );

	cells.template deviceToHost<0>();

	for (size_t i = 0 ; i < gr.size() - 1  ; i++)
	{
		size_t tot_p = starts.template get<0>(i+1) - starts.template get<0>(i);

		check &= (tot_p == cl.getNelements(i));

		grid_key_dx<dim> key = gr.InvLinId(i);

		// Now we check the ids

		for (size_t j = 0 ; j < cl.getNelements(i) ; j++)
		{
			size_t p_id = cl.get(i,j);

			size_t p_id2 = cells.template get<0>(starts.template get<0>(i) + j);

			check &= (p_id == p_id2);
		}
	}

	BOOST_REQUIRE(check == true);
}

BOOST_AUTO_TEST_CASE( test_subindex_funcs )
{
	std::cout << "Test cell list GPU base func" << "\n";

	test_sub_index<3,float,int,unsigned char>();
	test_sub_index2<3,float,int,unsigned char>();

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE ( test_cell_fill )
{
	std::cout << "Test GPU fill cells" << "\n";

	test_fill_cell<3,float,unsigned int, unsigned char>();

	std::cout << "End GPU fill cells" << "\n";
}

template<unsigned int dim, typename T, typename cnt_type, typename ids_type>
void test_reorder_parts(size_t n_part)
{
	// Create n_part
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> cells;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> cells_out;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> starts;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> sort_to_not_sort;
	openfpm::vector<aggregate<ids_type[dim+1]>,CudaMemory,typename memory_traits_inte<aggregate<ids_type[dim+1]>>::type,memory_traits_inte> part_ids;

	openfpm::vector<aggregate<float,float,float[3],float[3][3]>,CudaMemory,typename memory_traits_inte<aggregate<float,float,float[3],float[3][3]>>::type,memory_traits_inte> parts_prp;
	openfpm::vector<aggregate<float,float,float[3],float[3][3]>,CudaMemory,typename memory_traits_inte<aggregate<float,float,float[3],float[3][3]>>::type,memory_traits_inte> parts_prp_out;

	// CellList to check the result

	Box<dim,T> domain;

	typename boost::array_openfpm<T,dim> spacing;
	typename boost::array_openfpm<ids_type,dim,cnt_type> div_c;

	size_t div_host[dim];

	size_t tot = 1;
	for (size_t i = 0 ; i < dim ; i++)
	{
		div_host[i] = 10;
		div_c[i] = tot;
		tot *= div_host[i];
		spacing[i] = 0.1;
		domain.setLow(i,0.0);
		domain.setHigh(i,1.0);
	}

	CellList<dim,T, Mem_fast<>> cl(domain,div_host,0);
	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> pl;
	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> pl_out;

	create_n_part(n_part,pl,cl);
	parts_prp.resize(n_part);
	parts_prp_out.resize(n_part);
	pl_out.resize(n_part);
	sort_to_not_sort.resize(n_part);

	auto p_it = parts_prp.getIterator();
	while (p_it.isNext())
	{
		auto p = p_it.get();

		parts_prp.template get<0>(p) = 10000 + p;
		parts_prp.template get<1>(p) = 20000 + p;

		parts_prp.template get<2>(p)[0] = 30000 + p;
		parts_prp.template get<2>(p)[1] = 40000 + p;
		parts_prp.template get<2>(p)[2] = 50000 + p;

		parts_prp.template get<3>(p)[0][0] = 60000 + p;
		parts_prp.template get<3>(p)[0][1] = 70000 + p;
		parts_prp.template get<3>(p)[0][2] = 80000 + p;
		parts_prp.template get<3>(p)[1][0] = 90000 + p;
		parts_prp.template get<3>(p)[1][1] = 100000 + p;
		parts_prp.template get<3>(p)[1][2] = 110000 + p;
		parts_prp.template get<3>(p)[2][0] = 120000 + p;
		parts_prp.template get<3>(p)[2][1] = 130000 + p;
		parts_prp.template get<3>(p)[0][2] = 140000 + p;

		++p_it;
	}

	grid_sm<dim,void> gr(div_host);

	create_starts_and_parts_ids(cl,gr,pl.size(),tot,starts,part_ids,cells_out);


	auto itgg = pl.getGPUIterator();

	cells_out.template hostToDevice<0>();

	auto ite = pl.getGPUIterator();

	parts_prp.template hostToDevice<0,1,2,3>();

	// Here we test fill cell
	reorder_parts<decltype(parts_prp.toKernel()),
			      decltype(pl.toKernel()),
			      decltype(sort_to_not_sort.toKernel()),
			      cnt_type,
			      shift_ph<0,cnt_type>><<<ite.wthr,ite.thr>>>(pl.size(),
			                                                  parts_prp.toKernel(),
			                                                  parts_prp_out.toKernel(),
			                                                  pl.toKernel(),
			                                                  pl_out.toKernel(),
			                                                  sort_to_not_sort.toKernel(),
			                                                  static_cast<cnt_type *>(cells_out.template getDeviceBuffer<0>()));

	bool check = true;
	parts_prp_out.template deviceToHost<0>();
	sort_to_not_sort.template deviceToHost<0>();

	size_t st = 0;
	for (size_t i = 0 ; i < tot ; i++)
	{
		size_t n = cl.getNelements(i);

		for (size_t j = 0 ; j < n ; j++)
		{
			size_t p = cl.get(i,j);

			check &= parts_prp_out.template get<0>(st) == parts_prp.template get<0>(p);
			check &= sort_to_not_sort.template get<0>(st) == p;

			st++;
		}
	}


	BOOST_REQUIRE_EQUAL(check,true);
}

BOOST_AUTO_TEST_CASE ( test_reorder_particles )
{
	std::cout << "Test GPU reorder" << "\n";

	test_reorder_parts<3,float,unsigned int, unsigned char>(5000);

	std::cout << "End GPU reorder" << "\n";
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
	openfpm::vector<Point<dim,double>,CudaMemory,typename memory_traits_inte<Point<dim,double>>::type,memory_traits_inte> pl_out;

	openfpm::vector<aggregate<float,float[3],float[3][3]>,CudaMemory,typename memory_traits_inte<aggregate<float,float[3],float[3][3]>>::type,memory_traits_inte> pl_prp;
	openfpm::vector<aggregate<float,float[3],float[3][3]>,CudaMemory,typename memory_traits_inte<aggregate<float,float[3],float[3][3]>>::type,memory_traits_inte> pl_prp_out;


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

	pl_prp.resize(pl.size());
	pl_prp_out.resize(pl.size());
	pl_out.resize(pl.size());

	for (size_t i = 0 ; i < pl.size() ; i++)
	{
		pl_prp.template get<0>(i) = pl.template get<0>(i)[0];

		pl_prp.template get<1>(i)[0] = pl.template get<0>(i)[0]+100.0;
		pl_prp.template get<1>(i)[1] = pl.template get<0>(i)[1]+100.0;
		pl_prp.template get<1>(i)[2] = pl.template get<0>(i)[2]+100.0;

		pl_prp.template get<2>(i)[0][0] = pl.template get<0>(i)[0]+1000.0;
		pl_prp.template get<2>(i)[0][1] = pl.template get<0>(i)[1]+1000.0;
		pl_prp.template get<2>(i)[0][2] = pl.template get<0>(i)[2]+1000.0;

		pl_prp.template get<2>(i)[1][0] = pl.template get<0>(i)[0]+2000.0;
		pl_prp.template get<2>(i)[1][1] = pl.template get<0>(i)[1]+3000.0;
		pl_prp.template get<2>(i)[1][2] = pl.template get<0>(i)[2]+4000.0;

		pl_prp.template get<2>(i)[2][0] = pl.template get<0>(i)[0]+5000.0;
		pl_prp.template get<2>(i)[2][1] = pl.template get<0>(i)[1]+6000.0;
		pl_prp.template get<2>(i)[2][2] = pl.template get<0>(i)[2]+7000.0;
	}

	pl_prp.resize(pl.size());
	pl_prp_out.resize(pl.size());

	pl.template hostToDevice<0>();
	pl_prp.template hostToDevice<0,1,2>();

	cl2.template construct<decltype(pl),decltype(pl_prp)>(pl,pl_out,pl_prp,pl_prp_out);

	// Check

	pl_prp_out.deviceToHost<0>();
	pl_prp_out.deviceToHost<1>();
	pl_prp_out.deviceToHost<2>();

	////// Correct order ////////////

	openfpm::vector<Point<dim,double>,CudaMemory,typename memory_traits_inte<Point<dim,double>>::type,memory_traits_inte> pl_correct;

	pl_correct.add(p9);
	pl_correct.add(p1);
	pl_correct.add(p2);
	pl_correct.add(p3);
	pl_correct.add(p5);
	pl_correct.add(p4);
	pl_correct.add(p6);
	pl_correct.add(p7);
	pl_correct.add(p8);

	for (size_t i = 0 ; i < pl_correct.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<0>(i),(float)pl_correct.template get<0>(i)[0]);
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<1>(i)[0],(float)(pl_correct.template get<0>(i)[0]+100.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<1>(i)[1],(float)(pl_correct.template get<0>(i)[1]+100.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<1>(i)[2],(float)(pl_correct.template get<0>(i)[2]+100.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<2>(i)[0][0],(float)(pl_correct.template get<0>(i)[0] + 1000.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<2>(i)[0][1],(float)(pl_correct.template get<0>(i)[1] + 1000.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<2>(i)[0][2],(float)(pl_correct.template get<0>(i)[2] + 1000.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<2>(i)[1][0],(float)(pl_correct.template get<0>(i)[0] + 2000.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<2>(i)[1][1],(float)(pl_correct.template get<0>(i)[1] + 3000.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<2>(i)[1][2],(float)(pl_correct.template get<0>(i)[2] + 4000.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<2>(i)[2][0],(float)(pl_correct.template get<0>(i)[0] + 5000.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<2>(i)[2][1],(float)(pl_correct.template get<0>(i)[1] + 6000.0));
		BOOST_REQUIRE_EQUAL(pl_prp_out.template get<2>(i)[2][2],(float)(pl_correct.template get<0>(i)[2] + 7000.0));
	}

	// Check the sort to non sort buffer

	auto & vsrt = cl2.private_get_sort_to_not_sorted();
	vsrt.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(vsrt.size(),9);

	BOOST_REQUIRE_EQUAL(vsrt.template get<0>(0),8);
	BOOST_REQUIRE_EQUAL(vsrt.template get<0>(1),0);
	BOOST_REQUIRE_EQUAL(vsrt.template get<0>(2),1);
	BOOST_REQUIRE_EQUAL(vsrt.template get<0>(3),2);
	BOOST_REQUIRE_EQUAL(vsrt.template get<0>(4),4);
	BOOST_REQUIRE_EQUAL(vsrt.template get<0>(5),3);
	BOOST_REQUIRE_EQUAL(vsrt.template get<0>(6),5);
	BOOST_REQUIRE_EQUAL(vsrt.template get<0>(7),6);
	BOOST_REQUIRE_EQUAL(vsrt.template get<0>(8),7);

}


BOOST_AUTO_TEST_CASE( CellList_gpu_use)
{
	std::cout << "Test cell list GPU" << "\n";

	SpaceBox<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	SpaceBox<3,double> box2({-1.0f,-1.0f,-1.0f},{1.0f,1.0f,1.0f});

	Test_cell_gpu<3,double,CellList_gpu<3,double,CudaMemory>>(box);

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

template<unsigned int dim, typename vector_ps, typename vector_pr>
void fill_random_parts(Box<dim,float> & box, vector_ps & vd_pos, vector_pr & vd_prp, size_t n)
{
	for (size_t i = 0 ; i < n ; i++)
	{
		Point<dim,float> p;

		p.get(0) = ((box.getHigh(0) - box.getLow(0) - 0.0001)*(float)rand()/RAND_MAX) + box.getLow(0);
		p.get(1) = ((box.getHigh(1) - box.getLow(1) - 0.0001)*(float)rand()/RAND_MAX) + box.getLow(1);
		p.get(2) = ((box.getHigh(2) - box.getLow(2) - 0.0001)*(float)rand()/RAND_MAX) + box.getLow(2);

		vd_pos.add(p);
		vd_prp.add();
		vd_prp.last().template get<0>() = i % 3;
	}
}


template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_number(vector_pos pos, vector_ns s_t_ns, CellList_type cl, vector_n_type vn)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= pos.size()) return;

    vn.template get<0>(p) = 0;

    Point<3,float> xp = pos.template get<0>(p);

    auto it = cl.getNNIterator(cl.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get();

    	int s1 = s_t_ns.template get<0>(q);

    	atomicAdd(&vn.template get<0>(s1), 1);

    	++it;
    }
}


template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_list(vector_pos pos, vector_ns s_t_ns, CellList_type cl, vector_n_type v_nscan ,vector_n_type v_list)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);
    int start_list = v_nscan.template get<0>(p);

    auto it = cl.getNNIterator(cl.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get();

    	int s1 = s_t_ns.template get<0>(q);

    	v_list.template get<0>(start_list) = s1;

    	++start_list;
    	++it;
    }
}

template<unsigned int dim, typename T, typename CellS> void Test_cell_gpu_force(SpaceBox<dim,T> & box, size_t npart)
{
	// Subdivisions
	size_t div[dim] = {16,16,16};

	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cl2(box,div);

	// vector of particles

	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> pl;
	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> pl_out;

	openfpm::vector<aggregate<T,T[3]>,CudaMemory,typename memory_traits_inte<aggregate<T,T[3]>>::type,memory_traits_inte> pl_prp;
	openfpm::vector<aggregate<T,T[3]>,CudaMemory,typename memory_traits_inte<aggregate<T,T[3]>>::type,memory_traits_inte> pl_prp_out;

	openfpm::vector<aggregate<unsigned int>,CudaMemory,typename memory_traits_inte<aggregate<unsigned int>>::type,memory_traits_inte> n_out;

	// create random particles

	fill_random_parts<3>(box,pl,pl_prp,npart);

	pl_prp_out.resize(pl.size());
	pl_out.resize(pl.size());
	pl_out.resize(pl.size()+1);
	n_out.fill<0>(0);

	pl_prp.resize(pl.size());
	pl_prp_out.resize(pl.size());

	pl.template hostToDevice<0>();
	pl_prp.template hostToDevice<0,1>();

	// Construct an equivalent CPU cell-list

	CellList<dim,T,Mem_fast<>,shift<dim,T>> cl_cpu(box,div);

	// construct

	auto it2 = pl.getIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,T> xp = pl.get(p);

		cl_cpu.add(xp,p);

		++it2;
	}

	cl2.template construct<decltype(pl),decltype(pl_prp)>(pl,pl_out,pl_prp,pl_prp_out);
	auto & s_t_ns = cl2.getSortToNonSort();

	pl.template hostToDevice<0>();

	auto ite = pl.getGPUIterator();

	calc_force_number<decltype(pl.toKernel()),
			          decltype(s_t_ns.toKernel()),
			          decltype(cl2.toKernel()),
			          decltype(n_out.toKernel())>
	<<<ite.wthr,ite.thr>>>(pl.toKernel(),
						   s_t_ns.toKernel(),
						   cl2.toKernel(),
						   n_out.toKernel());

	// Check

	n_out.deviceToHost<0>();

	{
	bool check = true;
	auto it = pl.getIterator();

	while(it.isNext())
	{
		auto p = it.get();

		Point<dim,T> xp = pl.get(p);

		// Get NN iterator

		auto NN_it = cl_cpu.getNNIterator(cl_cpu.getCell(xp));

		size_t n_ele = 0;
		while (NN_it.isNext())
		{
			auto q = NN_it.get();

			n_ele++;

			++NN_it;
		}

		check &= n_ele == n_out.template get<0>(p);

		++it;
	}
	BOOST_REQUIRE_EQUAL(check,true);
	}

	// now we scan the buffer

	openfpm::vector<aggregate<unsigned int>,CudaMemory,typename memory_traits_inte<aggregate<unsigned int>>::type,memory_traits_inte> n_out_scan;
	openfpm::vector<aggregate<unsigned int>,CudaMemory,typename memory_traits_inte<aggregate<unsigned int>>::type,memory_traits_inte> nn_list;

	scan<unsigned int,unsigned char>(n_out,n_out_scan);
	n_out_scan.template deviceToHost<0>();

	if (n_out_scan.template get<0>(pl.size()) == 0)
	{return;}

	nn_list.resize(n_out_scan.template get<0>(pl.size()));

	// Now for each particle we construct the list

	pl.template hostToDevice<0>();

	calc_force_list<decltype(pl.toKernel()),
			          decltype(s_t_ns.toKernel()),
			          decltype(cl2.toKernel()),
			          decltype(nn_list.toKernel())>
	<<<ite.wthr,ite.thr>>>(pl.toKernel(),
						   s_t_ns.toKernel(),
						   cl2.toKernel(),
						   n_out_scan.toKernel(),
						   nn_list.toKernel());

	nn_list.template deviceToHost<0>();

	// Check

	n_out.deviceToHost<0>();

	{
	bool check = true;
	auto it = pl.getIterator();

	while(it.isNext())
	{
		auto p = it.get();

		Point<dim,T> xp = pl.get(p);

		// Get NN iterator

		openfpm::vector<int> cpu_list;

		auto NN_it = cl_cpu.getNNIterator(cl_cpu.getCell(xp));

		while (NN_it.isNext())
		{
			auto q = NN_it.get();

			cpu_list.add(q);

			++NN_it;
		}

		openfpm::vector<int> gpu_list;

		for (size_t i = n_out_scan.template get<0>(p) ; i < n_out_scan.template get<0>(p+1) ; i++)
		{
			gpu_list.add(nn_list.template get<0>(i));
		}

		// sort bost vector

		cpu_list.sort();
		gpu_list.sort();

		for (size_t j = 0 ; j < cpu_list.size() ; j++)
		{check &= cpu_list.get(j) == gpu_list.get(j);}

		++it;
	}

	BOOST_REQUIRE_EQUAL(check,true);

	}
}

BOOST_AUTO_TEST_CASE( CellList_gpu_use_calc_force)
{
	std::cout << "Test cell list GPU" << "\n";

	SpaceBox<3,float> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	SpaceBox<3,float> box2({-0.3f,-0.3f,-0.3f},{1.0f,1.0f,1.0f});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>>(box,1000);
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>>(box,10000);

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>>(box2,1000);
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>>(box2,10000);

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

template<typename CellList_type, typename Vector_type, typename Vector_out>
__global__ void cl_offload_gpu(CellList_type cl, Vector_type parts, Vector_out output)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= parts.size()) return;

    Point<3,float> xp = parts.template get<0>(p);

    output.template get<0>(p) = cl.getNelements(cl.getCell(xp));
}

template<typename CellList_type, typename Vector_type, typename Vector_scan_type, typename Vector_list_type>
__global__ void cl_offload_gpu_list(CellList_type cl, Vector_type parts, Vector_scan_type scan, Vector_list_type list)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= parts.size()) return;

    Point<3,float> xp = parts.template get<0>(p);

    int id = cl.getCell(xp);
    int n_ele = cl.getNelements(id);
    int start = scan.template get<0>(p);

    for (int j = 0 ; j < n_ele ; j++)
    {
    	list.template get<0>(start+j) = cl.get(id,j);
    }

}

BOOST_AUTO_TEST_CASE( CellList_use_cpu_offload_test )
{
	std::cout << "Test cell list offload gpu" << "\n";

	// Subdivisions
	size_t div[3] = {10,10,10};

	// grid info
	grid_sm<3,void> g_info(div);

	Box<3,float> box({-1.0,-1.0,-1.0},{1.0,1.0,1.0});

	// CellS = CellListM<dim,T,8>
	CellList<3,float,Mem_fast<CudaMemory,int>,shift<3,float>> cl1(box,div);

	openfpm::vector_gpu<Point<3,float>> v;
	openfpm::vector_gpu<aggregate<int>> os;
	v.resize(10000);
	os.resize(v.size());

	for (size_t i = 0 ; i < v.size() ; i++)
	{
		v.template get<0>(i)[0] = 2.0 * (float)rand() / RAND_MAX - 1.0;
		v.template get<0>(i)[1] = 2.0 * (float)rand() / RAND_MAX - 1.0;
		v.template get<0>(i)[2] = 2.0 * (float)rand() / RAND_MAX - 1.0;

		Point<3,float> xp = v.template get<0>(i);

		cl1.add(xp,i);
	}

	auto ite = v.getGPUIterator();

	cl1.hostToDevice();
	v.hostToDevice<0>();

	cl_offload_gpu<decltype(cl1.toKernel()),decltype(v.toKernel()),decltype(os.toKernel())><<<ite.wthr,ite.thr>>>(cl1.toKernel(),v.toKernel(),os.toKernel());

	os.deviceToHost<0>();

	bool match = true;
	for (size_t i = 0 ; i < os.size() ; i++)
	{
		Point<3,float> xp = v.template get<0>(i);

		match &= os.template get<0>(i) == cl1.getNelements(cl1.getCell(xp));
	}

	BOOST_REQUIRE_EQUAL(match,true);

	// now we scan the vector out

	openfpm::vector_gpu<aggregate<int>> os_scan;
	os_scan.resize(v.size());

	scan<int,int>(os,os_scan);

	os_scan.deviceToHost<0>();
	os.deviceToHost<0>(os.size()-1,os.size()-1);
	size_t size_list = os_scan.template get<0>(os_scan.size()-1) + os.template get<0>(os.size()-1);

	openfpm::vector_gpu<aggregate<int>> os_list;
	os_list.resize(size_list);

	cl_offload_gpu_list<decltype(cl1.toKernel()),decltype(v.toKernel()),
			            decltype(os_scan.toKernel()),decltype(os_list.toKernel())><<<ite.wthr,ite.thr>>>
			            (cl1.toKernel(),v.toKernel(),os_scan.toKernel(),os_list.toKernel());

	os_list.deviceToHost<0>();

	match = true;
	for (size_t i = 0 ; i < os.size() ; i++)
	{
		Point<3,float> xp = v.template get<0>(i);

		for (size_t j = 0 ; j < cl1.getNelements(cl1.getCell(xp)) ; j++)
		{
			match &= os_list.template get<0>(os_scan.template get<0>(i)+j) == cl1.get(cl1.getCell(xp),j);
		}
	}

	BOOST_REQUIRE_EQUAL(match,true);

	std::cout << "End cell list offload gpu" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_SUITE_END()

