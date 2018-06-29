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
#include "CellList.hpp"
#include "util/boost/boost_array_openfpm.hpp"
#include  "Point_test.hpp"

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

template<unsigned int dim, typename T>
void create_n_part(int n_part,
		           openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> & pl,
		           CellList<dim,T, Mem_fast> & cl)
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
void create_starts_and_parts_ids(CellList<dim,T, Mem_fast> & cl,
								 grid_sm<dim,void> & gr,
								 size_t n_part,
								 size_t n_cell,
								 openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> & starts,
								 openfpm::vector<aggregate<ids_type[dim+1]>,CudaMemory,typename memory_traits_inte<aggregate<ids_type[dim+1]>>::type,memory_traits_inte> & part_ids)
{
	// Construct starts and part_ids

	part_ids.resize(n_part);
	starts.resize(n_cell);

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
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> starts;
	openfpm::vector<aggregate<ids_type[dim+1]>,CudaMemory,typename memory_traits_inte<aggregate<ids_type[dim+1]>>::type,memory_traits_inte> part_ids;

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

	CellList<dim,T, Mem_fast> cl(domain,div_host,0);
	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> pl;

	create_n_part(5000,pl,cl);

	grid_sm<dim,void> gr(div_host);

	create_starts_and_parts_ids(cl,gr,pl.size(),tot,starts,part_ids);

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
	openfpm::vector<aggregate<cnt_type>,CudaMemory,typename memory_traits_inte<aggregate<cnt_type>>::type,memory_traits_inte> starts;
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

	CellList<dim,T, Mem_fast> cl(domain,div_host,0);
	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> pl;

	create_n_part(n_part,pl,cl);
	parts_prp.resize(n_part);
	parts_prp_out.resize(n_part);

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

	parts_prp_out.set(0,parts_prp,0);

	grid_sm<dim,void> gr(div_host);

	create_starts_and_parts_ids(cl,gr,pl.size(),tot,starts,part_ids);


	auto itgg = pl.getGPUIterator();

	starts.template hostToDevice<0>();

	// Here we test fill cell
	reorder_parts<decltype(parts_prp.template toGPU<0,1,2,3>()),cnt_type,shift_ph<0,cnt_type>><<<1,1>>>(pl.size(),
			                                                                                                               parts_prp.template toGPU<0,1,2,3>(),
			                                                                                                               parts_prp_out.template toGPU<>(),
			                                                                                                               static_cast<cnt_type *>(starts.template getDeviceBuffer<0>()));

	parts_prp_out.template deviceToHost<0>();

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

