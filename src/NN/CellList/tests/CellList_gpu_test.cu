/*
 * CellList_gpu_test.cpp
 *
 *  Created on: Jun 13, 2018
 *      Author: i-bird
 */

#include "util/cuda_util.hpp"
#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "NN/CellList/cuda/CellList_gpu.hpp"
#include "NN/CellList/CellList.hpp"
#include "util/boost/boost_array_openfpm.hpp"
#include  "Point_test.hpp"
#include "util/cuda_util.hpp"

BOOST_AUTO_TEST_SUITE( CellList_gpu_test )

template<unsigned int dim, typename T, typename cnt_type, typename ids_type>
void test_sub_index()
{
	openfpm::vector<aggregate<unsigned int,unsigned int>,CudaMemory,memory_traits_inte> cl_n;
	openfpm::vector<aggregate<unsigned int[2]>,CudaMemory,memory_traits_inte> cellIndex_LocalIndex;

	// fill with some particles

	openfpm::vector<Point<dim,T>,CudaMemory,memory_traits_inte> vPos;

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

	vPos.add(p1);
	vPos.add(p2);
	vPos.add(p3);
	vPos.add(p4);
	vPos.add(p5);
	vPos.add(p6);
	vPos.add(p7);
	vPos.add(p8);
	vPos.add(p9);
	vPos.add(p10);

	openfpm::array<T,dim> spacing;
	openfpm::array<ids_type,dim> div;
	openfpm::array<ids_type,dim> off;

	for (size_t i = 0 ; i < dim ; i++)
	{
		div[i] = 17;
		spacing[i] = 0.1;
		off[i] = 0;
	}

	cl_n.resize(17*17*17);
	cl_n.template fill<0>(0);
	//CUDA_SAFE(cudaMemset(cl_n.template getDeviceBuffer<0>(),0,cl_n.size()*sizeof(unsigned int)));

	cellIndex_LocalIndex.resize(vPos.size());

	size_t sz[3] = {17,17,17};
	grid_sm<3,void> gr(sz);

	auto ite = vPos.getGPUIterator();

	vPos.template hostToDevice<0>();

	Matrix<dim,T> mt;
	Point<dim,T> pt;

	no_transform_only<dim,T> t(mt,pt);


	CUDA_LAUNCH_DIM3((fill_cellIndex_LocalIndex<dim,T,ids_type,no_transform_only<dim,T>>),ite.wthr,ite.thr,div,
																	spacing,
																	off,
																	t,
																	vPos.size(),
																	(size_t)0,
																	vPos.toKernel(),
																	cl_n.toKernel(),
																	cellIndex_LocalIndex.toKernel());

	cl_n.template deviceToHost<0>();
	cellIndex_LocalIndex.template deviceToHost<0>();

	// I have to find != 0 in the cell with particles

	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(0)[0],gr.LinId({2,2,2}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(1)[0],gr.LinId({9,2,2}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(2)[0],gr.LinId({2,9,2}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(3)[0],gr.LinId({2,2,9}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(4)[0],gr.LinId({9,9,2}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(5)[0],gr.LinId({9,2,9}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(6)[0],gr.LinId({2,9,9}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(7)[0],gr.LinId({9,9,9}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(8)[0],gr.LinId({0,0,0}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(9)[0],gr.LinId({2,2,2}));

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
	openfpm::vector<aggregate<unsigned int,unsigned int>,CudaMemory,memory_traits_inte> cl_n;
	openfpm::vector<aggregate<unsigned int[2]>,CudaMemory,memory_traits_inte> cellIndex_LocalIndex;

	// fill with some particles

	openfpm::vector<Point<dim,T>,CudaMemory,memory_traits_inte> vPos;

	// Make the test more complicated the test to pass because make different capacity of vPos and cellIndex_LocalIndex
	vPos.resize(256);
	vPos.resize(0);

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


	vPos.add(p1);
	vPos.add(p2);
	vPos.add(p3);
	vPos.add(p4);
	vPos.add(p5);
	vPos.add(p6);
	vPos.add(p7);
	vPos.add(p8);
	vPos.add(p9);
	vPos.add(p10);

	openfpm::array<T,dim> spacing;
	openfpm::array<ids_type,dim> div;
	openfpm::array<ids_type,dim> off;

	for (size_t i = 0 ; i < dim ; i++)
	{
		div[i] = 17;
		spacing[i] = 0.1;
		off[i] = 0;
	}

	cl_n.resize(17*17*17);
	cl_n.template fill<0>(0);
//	CUDA_SAFE(cudaMemset(cl_n.template getDeviceBuffer<0>(),0,cl_n.size()*sizeof(unsigned int)));

	cellIndex_LocalIndex.resize(vPos.size());

	size_t sz[3] = {17,17,17};
	grid_sm<3,void> gr(sz);

	auto ite = vPos.getGPUIterator();

	vPos.template hostToDevice<0>();

	Matrix<dim,T> mt;

	shift_only<dim,T> t(mt,pt);


	CUDA_LAUNCH_DIM3((fill_cellIndex_LocalIndex<dim,T,ids_type,shift_only<dim,T>>),ite.wthr,ite.thr,div,
																	spacing,
																	off,
																	t,
																	vPos.size(),
																	(size_t)0,
																	vPos.toKernel(),
																	cl_n.toKernel(),
																	cellIndex_LocalIndex.toKernel());

	cl_n.template deviceToHost<0>();
	cellIndex_LocalIndex.template deviceToHost<0>();

	// I have to find != 0 in the cell with particles

	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(0)[0],gr.LinId({2,2,2}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(1)[0],gr.LinId({9,2,2}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(2)[0],gr.LinId({2,9,2}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(3)[0],gr.LinId({2,2,9}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(4)[0],gr.LinId({9,9,2}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(5)[0],gr.LinId({9,2,9}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(6)[0],gr.LinId({2,9,9}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(7)[0],gr.LinId({9,9,9}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(8)[0],gr.LinId({0,0,0}));
	BOOST_REQUIRE_EQUAL(cellIndex_LocalIndex.template get<0>(9)[0],gr.LinId({2,2,2}));

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
		           openfpm::vector<Point<dim,T>,CudaMemory,memory_traits_inte> & vPos,
		           CellList<dim,T, Mem_fast<>> & cellList)
{
	vPos.resize(n_part);

	auto it = vPos.getIterator();

	while(it.isNext())
	{
		auto p = it.get();

		vPos.template get<0>(p)[0] = (double)rand()/RAND_MAX;
		vPos.template get<0>(p)[1] = (double)rand()/RAND_MAX;
		vPos.template get<0>(p)[2] = (double)rand()/RAND_MAX;

		Point<dim,T> xp;
		xp.get(0) = vPos.template get<0>(p)[0];
		xp.get(1) = vPos.template get<0>(p)[1];
		xp.get(2) = vPos.template get<0>(p)[2];

		size_t c = cellList.getCell(xp);
		cellList.addCell(c,p);

		++it;
	}
}

template<unsigned int dim, typename T, typename cnt_type, typename ids_type>
void create_starts_and_parts_ids(CellList<dim,T, Mem_fast<>> & cellList,
								 grid_sm<dim,void> & gr,
								 size_t n_part,
								 size_t n_cell,
								 openfpm::vector<aggregate<cnt_type>,CudaMemory,memory_traits_inte> & starts,
								 openfpm::vector<aggregate<ids_type[2]>,CudaMemory,memory_traits_inte> & cellIndex_LocalIndex,
								 openfpm::vector<aggregate<cnt_type>,CudaMemory,memory_traits_inte> & cells)
{
	// Construct starts and cellIndex_LocalIndex

	cellIndex_LocalIndex.resize(n_part);
	starts.resize(n_cell);
	cells.resize(n_part);

	grid_key_dx_iterator<dim> itg(gr);

	size_t start = 0;

	while (itg.isNext())
	{
		auto cell = itg.get();

		size_t clin = gr.LinId(cell);

		for (size_t j = 0 ; j < cellList.getNelements(clin) ; j++)
		{
			size_t p_id = cellList.get(clin,j);

			cellIndex_LocalIndex.template get<0>(p_id)[0] = clin;

			cellIndex_LocalIndex.template get<0>(p_id)[1] = j;

			cells.template get<0>(start+j) = p_id;
		}
		starts.template get<0>(clin) = start;
		start += cellList.getNelements(clin);

		++itg;
	}
}

template<typename sparse_vector_type>
__global__ void construct_cells(sparse_vector_type sv, grid_sm<3,void> gs)
{
	sv.init();

	grid_key_dx<3> key1({5,5,5});
	grid_key_dx<3> key2({5,5,6});
	grid_key_dx<3> key3({5,6,5});
	grid_key_dx<3> key4({5,6,6});
	grid_key_dx<3> key5({6,5,5});
	grid_key_dx<3> key6({6,5,6});
	grid_key_dx<3> key7({6,6,5});
	grid_key_dx<3> key8({6,6,6});

	grid_key_dx<3> key9({7,7,7});

	grid_key_dx<3> key10({9,9,9});

	sv.template insert<0>(gs.LinId(key1)) = gs.LinId(key1);
	sv.template insert<0>(gs.LinId(key2)) = gs.LinId(key2);
	sv.template insert<0>(gs.LinId(key3)) = gs.LinId(key3);
	sv.template insert<0>(gs.LinId(key4)) = gs.LinId(key4);
	sv.template insert<0>(gs.LinId(key5)) = gs.LinId(key5);
	sv.template insert<0>(gs.LinId(key6)) = gs.LinId(key6);
	sv.template insert<0>(gs.LinId(key7)) = gs.LinId(key7);
	sv.template insert<0>(gs.LinId(key8)) = gs.LinId(key8);
	sv.template insert<0>(gs.LinId(key9)) = gs.LinId(key9);
	sv.template insert<0>(gs.LinId(key10)) = gs.LinId(key10);

	sv.flush_block_insert();
}

void test_cell_count_n()
{
	openfpm::vector_sparse_gpu<aggregate<int>> vs;
	openfpm::vector_gpu<aggregate<unsigned int>> cells_nn;
	openfpm::vector_gpu<aggregate<int>> cells_nn_test;

	vs.template setBackground<0>(-1);

	vs.setGPUInsertBuffer(1,32);

	size_t sz[] = {17,17,17};
	grid_sm<3,void> gs(sz);

	CUDA_LAUNCH_DIM3(construct_cells,1,1,vs.toKernel(),gs);

	gpu::ofp_context_t gpuContext;

	vs.flush<sadd_<0>>(gpuContext,flush_type::FLUSH_ON_DEVICE);

	cells_nn.resize(11);
	cells_nn.fill<0>(0);

	grid_key_dx<3> start({0,0,0});
	grid_key_dx<3> stop({2,2,2});
	grid_key_dx<3> middle({1,1,1});

	int mid = gs.LinId(middle);

	grid_key_dx_iterator_sub<3> it(gs,start,stop);

	while (it.isNext())
	{
		auto p = it.get();

		cells_nn_test.add();
		cells_nn_test.get<0>(cells_nn_test.size()-1) = (int)gs.LinId(p) - mid;

		++it;
	}

	cells_nn_test.template hostToDevice<0>();

	auto itgg = vs.getGPUIterator();
	CUDA_LAUNCH((countNonEmptyNeighborCells),itgg,vs.toKernel(),cells_nn.toKernel(),cells_nn_test.toKernel());

	cells_nn.deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(0),8);
	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(1),8);
	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(2),8);
	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(3),8);
	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(4),8);
	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(5),8);
	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(6),8);
	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(7),9);
	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(8),2);
	BOOST_REQUIRE_EQUAL(cells_nn.template get<0>(9),1);

	// now we scan
	openfpm::scan((unsigned int *)cells_nn.template getDeviceBuffer<0>(), cells_nn.size(), (unsigned int *)cells_nn.template getDeviceBuffer<0>() , gpuContext);

	openfpm::vector_gpu<aggregate<unsigned int,unsigned int>> cell_nn_list;
	cell_nn_list.resize(7*8 + 9 + 2 + 1);

	CUDA_LAUNCH((fillNeighborCellList),itgg,vs.toKernel(),cells_nn.toKernel(),cells_nn_test.toKernel(),cell_nn_list.toKernel(),200);

	cell_nn_list.deviceToHost<0>();

	// 8 NN
	for (size_t i = 0 ; i < 7 ; i++)
	{
		BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*i+0),1535);
		BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*i+1),1536);
		BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*i+2),1552);
		BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*i+3),1553);
		BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*i+4),1824);
		BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*i+5),1825);
		BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*i+6),1841);
		BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*i+7),1842);
	}

	// 9 NN
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+0),1535);
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+1),1536);
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+2),1552);
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+3),1553);
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+4),1824);
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+5),1825);
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+6),1841);
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+7),1842);
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+8),2149);

	// 2 NN
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+9),1842);
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+9+1),2149);

	// 1 NN
	BOOST_REQUIRE_EQUAL(cell_nn_list.template get<0>(8*7+9+2),2763);
}

BOOST_AUTO_TEST_CASE( test_count_nn_cells )
{
	std::cout << "Test cell count nn" << std::endl;

	test_cell_count_n();
}

BOOST_AUTO_TEST_CASE( test_subindex_funcs )
{
	std::cout << "Test cell list GPU base func" << "\n";

	test_sub_index<3,float,int,unsigned char>();
	test_sub_index2<3,float,int,unsigned char>();

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

template<unsigned int dim, typename T, typename cnt_type, typename ids_type>
void test_reorder_parts(size_t n_part)
{
	// Create n_part
	openfpm::vector<aggregate<cnt_type>,CudaMemory,memory_traits_inte> cells;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,memory_traits_inte> cells_out;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,memory_traits_inte> starts;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,memory_traits_inte> sortToNonSort;
	openfpm::vector<aggregate<cnt_type>,CudaMemory,memory_traits_inte> NonSortToSort;
	openfpm::vector<aggregate<ids_type[2]>,CudaMemory,memory_traits_inte> cellIndex_LocalIndex;

	openfpm::vector<aggregate<float,float,float[3],float[3][3]>,CudaMemory,memory_traits_inte> vPrp;
	openfpm::vector<aggregate<float,float,float[3],float[3][3]>,CudaMemory,memory_traits_inte> vPrpReorder;

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

	CellList<dim,T, Mem_fast<>> cellList(domain,div_host,0);
	openfpm::vector<Point<dim,T>,CudaMemory,memory_traits_inte> vPos;

	create_n_part(n_part,vPos,cellList);
	vPrp.resize(n_part);
	vPrpReorder.resize(n_part);
	sortToNonSort.resize(n_part);
	NonSortToSort.resize(n_part);

	auto p_it = vPrp.getIterator();
	while (p_it.isNext())
	{
		auto p = p_it.get();

		vPrp.template get<0>(p) = 10000 + p;
		vPrp.template get<1>(p) = 20000 + p;

		vPrp.template get<2>(p)[0] = 30000 + p;
		vPrp.template get<2>(p)[1] = 40000 + p;
		vPrp.template get<2>(p)[2] = 50000 + p;

		vPrp.template get<3>(p)[0][0] = 60000 + p;
		vPrp.template get<3>(p)[0][1] = 70000 + p;
		vPrp.template get<3>(p)[0][2] = 80000 + p;
		vPrp.template get<3>(p)[1][0] = 90000 + p;
		vPrp.template get<3>(p)[1][1] = 100000 + p;
		vPrp.template get<3>(p)[1][2] = 110000 + p;
		vPrp.template get<3>(p)[2][0] = 120000 + p;
		vPrp.template get<3>(p)[2][1] = 130000 + p;
		vPrp.template get<3>(p)[0][2] = 140000 + p;

		++p_it;
	}

	grid_sm<dim,void> gr(div_host);

	create_starts_and_parts_ids(cellList,gr,vPos.size(),tot,starts,cellIndex_LocalIndex,cells_out);


	auto itgg = vPos.getGPUIterator();

	cells_out.template hostToDevice<0>();

	auto ite = vPos.getGPUIterator();

	vPrp.template hostToDevice<0,1,2,3>();

	CUDA_LAUNCH_DIM3((constructSortUnsortBidirectMap),
		ite.wthr,ite.thr,
		sortToNonSort.toKernel(),
		NonSortToSort.toKernel(),
		cells_out.toKernel()
	);

	CUDA_LAUNCH_DIM3(
		(reorderParticlesPrp<
			decltype(vPrp.toKernel()),
			decltype(NonSortToSort.toKernel()),
			0>),
		ite.wthr,ite.thr,
		vPrp.toKernel(),
		vPrpReorder.toKernel(),
		NonSortToSort.toKernel(),
		(size_t)0
	);

	bool check = true;
	vPrpReorder.template deviceToHost<0>();
	sortToNonSort.template deviceToHost<0>();
	NonSortToSort.template deviceToHost<0>();

	size_t st = 0;
	for (size_t i = 0 ; i < tot ; i++)
	{
		size_t n = cellList.getNelements(i);

		for (size_t j = 0 ; j < n ; j++)
		{
			size_t p = cellList.get(i,j);

			check &= vPrpReorder.template get<0>(st) == vPrp.template get<0>(p);
			check &= sortToNonSort.template get<0>(st) == p;
			check &= NonSortToSort.template get<0>(p) == st;

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

template<unsigned int dim, typename T, typename CellS> void Test_cell_gpu(Box<dim,T> & box)
{
	//Space where is living the Cell list
	//Box<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Subdivisions
	size_t div[dim] = {16,16,16};

	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cellList2(box,div);
	cellList2.setOpt(CL_NON_SYMMETRIC | CL_GPU_REORDER_PROPERTY | CL_GPU_RESTORE_PROPERTY);

	// vector of particles

	openfpm::vector<Point<dim,double>,CudaMemory,memory_traits_inte> vPos;
	openfpm::vector<aggregate<float,float[3],float[3][3]>,CudaMemory,memory_traits_inte> vPrp;

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

	vPos.add(p1);
	vPos.add(p2);
	vPos.add(p3);
	vPos.add(p4);
	vPos.add(p5);
	vPos.add(p6);
	vPos.add(p7);
	vPos.add(p8);
	vPos.add(p9);

	vPrp.resize(vPos.size());

	for (size_t i = 0 ; i < vPos.size() ; i++)
	{
		vPrp.template get<0>(i) = vPos.template get<0>(i)[0];

		vPrp.template get<1>(i)[0] = vPos.template get<0>(i)[0]+100.0;
		vPrp.template get<1>(i)[1] = vPos.template get<0>(i)[1]+100.0;
		vPrp.template get<1>(i)[2] = vPos.template get<0>(i)[2]+100.0;

		vPrp.template get<2>(i)[0][0] = vPos.template get<0>(i)[0]+1000.0;
		vPrp.template get<2>(i)[0][1] = vPos.template get<0>(i)[1]+1000.0;
		vPrp.template get<2>(i)[0][2] = vPos.template get<0>(i)[2]+1000.0;

		vPrp.template get<2>(i)[1][0] = vPos.template get<0>(i)[0]+2000.0;
		vPrp.template get<2>(i)[1][1] = vPos.template get<0>(i)[1]+3000.0;
		vPrp.template get<2>(i)[1][2] = vPos.template get<0>(i)[2]+4000.0;

		vPrp.template get<2>(i)[2][0] = vPos.template get<0>(i)[0]+5000.0;
		vPrp.template get<2>(i)[2][1] = vPos.template get<0>(i)[1]+6000.0;
		vPrp.template get<2>(i)[2][2] = vPos.template get<0>(i)[2]+7000.0;
	}

	vPrp.resize(vPos.size());

	vPos.template hostToDevice<0>();
	vPrp.template hostToDevice<0,1,2>();

	openfpm::vector<Point<dim,double>,CudaMemory,memory_traits_inte> vPosReorder(vPos.size());
	openfpm::vector<aggregate<float,float[3],float[3][3]>,CudaMemory,memory_traits_inte> vPrpReorder(vPrp.size());

	// create an gpu context
	gpu::ofp_context_t gpuContext(gpu::gpu_context_opt::no_print_props);

	cellList2.template construct<decltype(vPos), decltype(vPrp), 0, 1, 2>(
		vPos,
		vPrp,
		vPosReorder,
		vPrpReorder,
		gpuContext,
		vPos.size(),
		0,
		vPos.size()
	);

	vPrpReorder.template deviceToHost<0,1,2>();

	////// Correct order ////////////

	openfpm::vector<Point<dim,double>,CudaMemory,memory_traits_inte> pl_correct;

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
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<0>(i),(float)pl_correct.template get<0>(i)[0]);
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<1>(i)[0],(float)(pl_correct.template get<0>(i)[0]+100.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<1>(i)[1],(float)(pl_correct.template get<0>(i)[1]+100.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<1>(i)[2],(float)(pl_correct.template get<0>(i)[2]+100.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<2>(i)[0][0],(float)(pl_correct.template get<0>(i)[0] + 1000.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<2>(i)[0][1],(float)(pl_correct.template get<0>(i)[1] + 1000.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<2>(i)[0][2],(float)(pl_correct.template get<0>(i)[2] + 1000.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<2>(i)[1][0],(float)(pl_correct.template get<0>(i)[0] + 2000.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<2>(i)[1][1],(float)(pl_correct.template get<0>(i)[1] + 3000.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<2>(i)[1][2],(float)(pl_correct.template get<0>(i)[2] + 4000.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<2>(i)[2][0],(float)(pl_correct.template get<0>(i)[0] + 5000.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<2>(i)[2][1],(float)(pl_correct.template get<0>(i)[1] + 6000.0));
		BOOST_REQUIRE_EQUAL(vPrpReorder.template get<2>(i)[2][2],(float)(pl_correct.template get<0>(i)[2] + 7000.0));
	}

	// Check the sort to non sort buffer

	auto & vsrt = cellList2.getSortToNonSort();
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

	auto & vnsrt = cellList2.getNonSortToSort();

	BOOST_REQUIRE_EQUAL(vnsrt.size(),9);

	// Move to CPU

	vnsrt.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(vnsrt.template get<0>(8),0);
	BOOST_REQUIRE_EQUAL(vnsrt.template get<0>(0),1);
	BOOST_REQUIRE_EQUAL(vnsrt.template get<0>(1),2);
	BOOST_REQUIRE_EQUAL(vnsrt.template get<0>(2),3);
	BOOST_REQUIRE_EQUAL(vnsrt.template get<0>(4),4);
	BOOST_REQUIRE_EQUAL(vnsrt.template get<0>(3),5);
	BOOST_REQUIRE_EQUAL(vnsrt.template get<0>(5),6);
	BOOST_REQUIRE_EQUAL(vnsrt.template get<0>(6),7);
	BOOST_REQUIRE_EQUAL(vnsrt.template get<0>(7),8);
}


BOOST_AUTO_TEST_CASE( CellList_gpu_use)
{
	std::cout << "Test cell list GPU" << "\n";

	Box<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Box<3,double> box2({-1.0f,-1.0f,-1.0f},{1.0f,1.0f,1.0f});

	Test_cell_gpu<3,double,CellList_gpu<3,double,CudaMemory>>(box);

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( CellList_gpu_use_sparse )
{
	std::cout << "Test cell list GPU sparse" << "\n";

	Box<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Box<3,double> box2({-1.0f,-1.0f,-1.0f},{1.0f,1.0f,1.0f});

	Test_cell_gpu<3,double,CellList_gpu<3,double,CudaMemory,no_transform_only<3,double>,true>> (box);

	std::cout << "End cell list GPU sparse" << "\n";

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
__global__ void calc_force_number(vector_pos pos, vector_ns sortToNonSort, CellList_type cellList, vector_n_type n_out)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);

    auto it = cellList.getNNIteratorBox(cellList.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get_sort();
    	auto q_ns = it.get();

		int s1 = sortToNonSort.template get<0>(q);

		atomicAdd(&n_out.template get<0>(s1), 1);

    	++it;
    }
}

template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_number_noato(vector_pos pos, vector_ns sortToNonSort, CellList_type cellList, vector_n_type n_out)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);

    auto it = cellList.getNNIteratorBox(cellList.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get_sort();
    	auto q_ns = it.get();

		int s1 = sortToNonSort.template get<0>(q);

		++n_out.template get<0>(p);

    	++it;
    }
}

template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_number_box(vector_pos pos, vector_ns sortToNonSort, CellList_type cellList, vector_n_type n_out, unsigned int start)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x + start;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);

    auto it = cellList.getNNIteratorBox(cellList.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get_sort();

		int s1 = sortToNonSort.template get<0>(q);

		atomicAdd(&n_out.template get<0>(s1), 1);

    	++it;
    }
}

template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_number_box_noato(vector_pos pos, vector_ns sortToNonSort, CellList_type cellList, vector_n_type n_out, unsigned int start)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x + start;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);

    auto it = cellList.getNNIteratorBox(cellList.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get_sort();

		++n_out.template get<0>(p);

    	++it;
    }
}

template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_number_rad(vector_pos pos, vector_ns sortToNonSort, CellList_type cellList, vector_n_type n_out)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);

    auto it = cellList.getNNIteratorRadius(cellList.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get_sort();

		int s1 = sortToNonSort.template get<0>(q);

		atomicAdd(&n_out.template get<0>(s1), 1);

    	++it;
    }
}

template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_list_box(vector_pos pos, vector_ns sortToNonSort, CellList_type cellList, vector_n_type v_nscan ,vector_n_type v_list)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);
    int start_list = v_nscan.template get<0>(p);

    auto it = cellList.getNNIteratorBox(cellList.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get_sort();

		int s1 = sortToNonSort.template get<0>(q);

    	v_list.template get<0>(start_list) = s1;

    	++start_list;
    	++it;
    }
}

template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_list(vector_pos pos, vector_ns sortToNonSort, CellList_type cellList, vector_n_type v_nscan ,vector_n_type v_list)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);
    int start_list = v_nscan.template get<0>(p);

    auto it = cellList.getNNIteratorBox(cellList.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get_sort();

		int s1 = sortToNonSort.template get<0>(q);

    	v_list.template get<0>(start_list) = s1;

    	++start_list;
    	++it;
    }
}

template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_list_box_partial(vector_pos pos,
										vector_ns sortToNonSort,
										CellList_type cellList,
										vector_n_type v_nscan,
										vector_n_type v_nscan_part,
										vector_n_type v_list)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);
    int start_list = v_nscan.template get<0>(p) + v_nscan_part.template get<0>(p);

    auto it = cellList.getNNIteratorBox(cellList.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get_sort();

		int s1 = sortToNonSort.template get<0>(q);

    	v_list.template get<0>(start_list) = s1;

    	++start_list;
    	++it;
    }
}

template<typename vector_pos, typename vector_ns, typename CellList_type,typename vector_n_type>
__global__ void calc_force_list_rad(vector_pos pos, vector_ns sortToNonSort, CellList_type cellList, vector_n_type v_nscan ,vector_n_type v_list)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= pos.size()) return;

    Point<3,float> xp = pos.template get<0>(p);
    int start_list = v_nscan.template get<0>(p);

    auto it = cellList.getNNIteratorRadius(cellList.getCell(xp));

    while (it.isNext())
    {
    	auto q = it.get_sort();

		int s1 = sortToNonSort.template get<0>(q);

    	v_list.template get<0>(start_list) = s1;

    	++start_list;
    	++it;
    }
}

template<unsigned int impl>
struct execute_cl_test
{
	template<typename CellS, typename Cells_cpu_type, typename T>
	static void set_radius(CellS & cellList2, Cells_cpu_type & cl_cpu, T & radius)
	{
	}

	template<typename pl_type, typename sortToNonSort_type, typename cl2_type, typename n_out_type>
	static void calc_num(pl_type & vPos, sortToNonSort_type & sortToNonSort, cl2_type & cellList2, n_out_type & n_out, unsigned int start)
	{
		auto ite = vPos.getGPUIterator();

		CUDA_LAUNCH((calc_force_number),ite,vPos.toKernel(),
			sortToNonSort.toKernel(),
			cellList2.toKernel(),
			n_out.toKernel()
		);
	}

	template<typename pl_type, typename sortToNonSort_type, typename cl2_type, typename n_out_scan_type, typename nn_list_type>
	static void calc_list(pl_type & vPos, sortToNonSort_type & sortToNonSort, cl2_type & cellList2,n_out_scan_type & n_out_scan, nn_list_type & nn_list)
	{
		auto ite = vPos.getGPUIterator();

		CUDA_LAUNCH((calc_force_list),ite,vPos.toKernel(),
			sortToNonSort.toKernel(),
			cellList2.toKernel(),
			n_out_scan.toKernel(),
			nn_list.toKernel()
		);
	}

	template<typename NN_type>
	static auto getNN(NN_type & nn, size_t cell) -> decltype(nn.getNNIteratorBox(cell))
	{
		return nn.getNNIteratorBox(cell);
	}
};

template<>
struct execute_cl_test<1>
{
	template<typename CellS, typename Cells_cpu_type, typename T>
	static void set_radius(CellS & cellList2, Cells_cpu_type & cl_cpu, T & radius)
	{
		cellList2.setRadius(radius);
		cl_cpu.setRadius(radius);
	}

	template<typename pl_type, typename sortToNonSort_type, typename cl2_type, typename n_out_type>
	static void calc_num(pl_type & vPos, sortToNonSort_type & sortToNonSort, cl2_type & cellList2, n_out_type & n_out, unsigned int start)
	{
		auto ite = vPos.getGPUIterator();

		CUDA_LAUNCH((calc_force_number_rad<decltype(vPos.toKernel()),
				          decltype(sortToNonSort.toKernel()),
				          decltype(cellList2.toKernel()),
				          decltype(n_out.toKernel())>),
							   ite,vPos.toKernel(),
							   sortToNonSort.toKernel(),
							   cellList2.toKernel(),
							   n_out.toKernel());
	}

	template<typename pl_type, typename sortToNonSort_type, typename cl2_type, typename n_out_scan_type, typename nn_list_type>
	static void calc_list(pl_type & vPos, sortToNonSort_type & sortToNonSort, cl2_type & cellList2, n_out_scan_type & n_out_scan, nn_list_type & nn_list)
	{
		auto ite = vPos.getGPUIterator();

		CUDA_LAUNCH((calc_force_list_rad<decltype(vPos.toKernel()),
				          decltype(sortToNonSort.toKernel()),
				          decltype(cellList2.toKernel()),
				          decltype(nn_list.toKernel())>),
							   ite,vPos.toKernel(),
							   sortToNonSort.toKernel(),
							   cellList2.toKernel(),
							   n_out_scan.toKernel(),
							   nn_list.toKernel());
	}

	template<typename NN_type>
	static auto getNN(NN_type & nn, size_t cell) -> decltype(nn.getNNIteratorRadius(cell))
	{
		return nn.getNNIteratorRadius(cell);
	}
};

template<>
struct execute_cl_test<2>
{
	template<typename CellS, typename Cells_cpu_type, typename T>
	static void set_radius(CellS & cellList2, Cells_cpu_type & cl_cpu, T & radius)
	{
		cellList2.setRadius(radius);
		cl_cpu.setRadius(radius);
	}

	template<typename pl_type, typename sortToNonSort_type, typename cl2_type, typename n_out_type>
	static void calc_num_noato(pl_type & vPos, sortToNonSort_type & sortToNonSort, cl2_type & cellList2, n_out_type & n_out, unsigned int start)
	{
		auto ite = sortToNonSort.getGPUIterator();

		CUDA_LAUNCH((calc_force_number_box_noato<decltype(vPos.toKernel()),
				          decltype(sortToNonSort.toKernel()),
				          decltype(cellList2.toKernel()),
				          decltype(n_out.toKernel())>),
							   ite,vPos.toKernel(),
							   sortToNonSort.toKernel(),
							   cellList2.toKernel(),
							   n_out.toKernel(),
							   start);
	}

	template<typename pl_type, typename sortToNonSort_type, typename cl2_type, typename n_out_type>
	static void calc_num(pl_type & vPos, sortToNonSort_type & sortToNonSort, cl2_type & cellList2, n_out_type & n_out, unsigned int start)
	{
		auto ite = sortToNonSort.getGPUIterator();

		CUDA_LAUNCH((calc_force_number_box<decltype(vPos.toKernel()),
				          decltype(sortToNonSort.toKernel()),
				          decltype(cellList2.toKernel()),
				          decltype(n_out.toKernel())>),
							   ite,
							   vPos.toKernel(),
							   sortToNonSort.toKernel(),
							   cellList2.toKernel(),
							   n_out.toKernel(),
							   start);
	}

	template<typename pl_type, typename sortToNonSort_type, typename cl2_type, typename n_out_scan_type, typename nn_list_type>
	static void calc_list(pl_type & vPos, sortToNonSort_type & sortToNonSort, cl2_type & cellList2, n_out_scan_type & n_out_scan, nn_list_type & nn_list)
	{
		auto ite = sortToNonSort.getGPUIterator();

		CUDA_LAUNCH((calc_force_list_box<decltype(vPos.toKernel()),
				          decltype(sortToNonSort.toKernel()),
				          decltype(cellList2.toKernel()),
				          decltype(nn_list.toKernel())>),
							   ite,vPos.toKernel(),
							   sortToNonSort.toKernel(),
							   cellList2.toKernel(),
							   n_out_scan.toKernel(),
							   nn_list.toKernel());
	}

	template<typename pl_type, typename sortToNonSort_type, typename cl2_type, typename n_out_scan_type, typename nn_list_type>
	static void calc_list_partial(pl_type & vPos,
		sortToNonSort_type & sortToNonSort,
		cl2_type & cellList2,
		n_out_scan_type & n_out_scan,
		n_out_scan_type & n_out_scan_partial,
		nn_list_type & nn_list)
	{
		auto ite = sortToNonSort.getGPUIterator();

		CUDA_LAUNCH((calc_force_list_box_partial),ite,vPos.toKernel(),
							   	   	   sortToNonSort.toKernel(),
							   	   	   cellList2.toKernel(),
							   	   	   n_out_scan.toKernel(),
							   	   	   n_out_scan_partial.toKernel(),
							   	   	   nn_list.toKernel());
	}

	template<typename NN_type>
	static auto getNN(NN_type & nn, size_t cell) -> decltype(nn.getNNIteratorRadius(cell))
	{
		return nn.getNNIteratorRadius(cell);
	}
};

template<unsigned int dim, typename T, typename CellS, int impl>
void Test_cell_gpu_force(Box<dim,T> & box, size_t npart, const size_t (& div)[dim],int box_nn = 2)
{
	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cellList2(box,div,2);

	CellList<dim,T,Mem_fast<>,shift<dim,T>> cl_cpu(box,div,2);

	cellList2.setBoxNN(box_nn);

	T radius = (box.getHigh(0) - box.getLow(0))/div[0] * 2.0;
	execute_cl_test<impl>::set_radius(cellList2,cl_cpu,radius);

	// vector of particles

	openfpm::vector<Point<dim,T>,CudaMemory,memory_traits_inte> vPos;
	openfpm::vector<aggregate<T,T[3]>,CudaMemory,memory_traits_inte> vPrp;

	openfpm::vector<aggregate<unsigned int>,CudaMemory,memory_traits_inte> n_out;

	// create random particles

	fill_random_parts<3>(box,vPos,vPrp,npart);

	n_out.resize(vPos.size()+1);
	n_out.fill<0>(0);

	vPrp.resize(vPos.size());

	vPos.template hostToDevice<0>();
	vPrp.template hostToDevice<0,1>();

	// Construct an equivalent CPU cell-list

	// construct

	auto it2 = vPos.getIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,T> xp = vPos.get(p);

		cl_cpu.add(xp,p);

		++it2;
	}

	size_t ghostMarker = vPos.size() / 2;

	gpu::ofp_context_t gpuContext(gpu::gpu_context_opt::no_print_props);
	cellList2.construct(vPos, vPrp, gpuContext, ghostMarker, 0, vPos.size());

	auto & sortToNonSort = cellList2.getSortToNonSort();

	vPos.template hostToDevice<0>();

	execute_cl_test<impl>::calc_num(vPos,sortToNonSort,cellList2,n_out,0);

	// Domain particles

	auto & gdsi = cellList2.getDomainSortIds();
	gdsi.template deviceToHost<0>();
	sortToNonSort.template deviceToHost<0>();

	bool match = true;
	for (size_t i = 0 ; i < ghostMarker ; i++)
	{
		unsigned int p = gdsi.template get<0>(i);

		match &= (sortToNonSort.template get<0>(p) < ghostMarker);
	}

	BOOST_REQUIRE_EQUAL(match,true);

	// Check

	n_out.deviceToHost<0>();

	{
		bool check = true;
		auto it = vPos.getIterator();

		while(it.isNext())
		{
			auto p = it.get();

			Point<dim,T> xp = vPos.get(p);

			// Get NN iterator

			auto NN_it = execute_cl_test<impl>::getNN(cl_cpu,cl_cpu.getCell(xp)); /*cl_cpu.getNNIteratorBox(cl_cpu.getCell(xp))*/;

			size_t n_ele = 0;
			while (NN_it.isNext())
			{
				auto q = NN_it.get();

				n_ele++;

				++NN_it;
			}

			check &= n_ele == n_out.template get<0>(p);

			if (check == false)
			{
				std::cout << p << "  " << n_ele << "   " << n_out.template get<0>(p) << "   " << check << std::endl;
				break;
			}

			++it;
		}
		BOOST_REQUIRE_EQUAL(check,true);
	}

	// now we scan the buffer

	openfpm::vector<aggregate<unsigned int>,CudaMemory,memory_traits_inte> n_out_scan;
	openfpm::vector<aggregate<unsigned int>,CudaMemory,memory_traits_inte> nn_list;

	n_out_scan.resize(vPos.size()+1);

	openfpm::scan((unsigned int *)n_out.template getDeviceBuffer<0>(),n_out.size(),(unsigned int *)n_out_scan.template getDeviceBuffer<0>(),gpuContext);
	n_out_scan.template deviceToHost<0>();

	if (n_out_scan.template get<0>(vPos.size()) == 0)
	{return;}

	nn_list.resize(n_out_scan.template get<0>(vPos.size()));

	// Now for each particle we construct the list

	vPos.template hostToDevice<0>();

	execute_cl_test<impl>::calc_list(vPos,sortToNonSort,cellList2,n_out_scan,nn_list);

	nn_list.template deviceToHost<0>();

	// Check

	n_out.deviceToHost<0>();

	{
	bool check = true;
	auto it = vPos.getIterator();

	while(it.isNext())
	{
		auto p = it.get();

		Point<dim,T> xp = vPos.get(p);

		// Get NN iterator

		openfpm::vector<int> cpu_list;

		auto NN_it = execute_cl_test<impl>::getNN(cl_cpu,cl_cpu.getCell(xp)); /*cl_cpu.getNNIteratorBox(cl_cpu.getCell(xp));*/

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

template<unsigned int dim, typename T, typename CellS, int impl>
void Test_cell_gpu_force_split(Box<dim,T> & box, size_t npart, const size_t (& div)[dim],int box_nn = 2)
{
	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cl2_split1(box,div,2);
	CellS cl2_split2(box,div,2);

	CellList<dim,T,Mem_fast<>,shift<dim,T>> cl_cpu(box,div,2);

	cl2_split1.setBoxNN(box_nn);
	cl2_split2.setBoxNN(box_nn);

	T radius = (box.getHigh(0) - box.getLow(0))/div[0] * 2.0;
	execute_cl_test<impl>::set_radius(cl2_split1,cl_cpu,radius);
	execute_cl_test<impl>::set_radius(cl2_split2,cl_cpu,radius);

	// vector of particles

	openfpm::vector<Point<dim,T>,CudaMemory,memory_traits_inte> vPos;
	openfpm::vector<aggregate<T,T[3]>,CudaMemory,memory_traits_inte> vPrp;

	openfpm::vector<aggregate<unsigned int>,CudaMemory,memory_traits_inte> n_out;
	openfpm::vector<aggregate<unsigned int>,CudaMemory,memory_traits_inte> n_out_partial;

	// create random particles

	fill_random_parts<3>(box,vPos,vPrp,npart);

	n_out.resize(vPos.size()+1);
	n_out.fill<0>(0);

	vPrp.resize(vPos.size());

	vPos.template hostToDevice<0>();
	vPrp.template hostToDevice<0,1>();

	// Construct an equivalent CPU cell-list

	// construct

	auto it2 = vPos.getIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,T> xp = vPos.get(p);

		cl_cpu.add(xp,p);

		++it2;
	}

	size_t ghostMarker = vPos.size() / 2;

	gpu::ofp_context_t gpuContext(gpu::gpu_context_opt::no_print_props);
	cl2_split1.construct(vPos,vPrp,gpuContext,ghostMarker,0,vPos.size()/2);
	cl2_split2.construct(vPos,vPrp,gpuContext,ghostMarker,vPos.size()/2,vPos.size());
	auto & sortToNonSort_s1 = cl2_split1.getSortToNonSort();
	auto & sortToNonSort_s2 = cl2_split2.getSortToNonSort();

	execute_cl_test<impl>::calc_num_noato(vPos,sortToNonSort_s1,cl2_split1,n_out,0);
	n_out_partial = n_out;
	execute_cl_test<impl>::calc_num_noato(vPos,sortToNonSort_s2,cl2_split2,n_out,0);

	// Domain particles

	auto & gdsi_s1 = cl2_split1.getDomainSortIds();
	gdsi_s1.template deviceToHost<0>();
	sortToNonSort_s1.template deviceToHost<0>();

	bool match = true;
	for (size_t i = 0 ; i < ghostMarker ; i++)
	{
		unsigned int p = gdsi_s1.template get<0>(i);

		match &= (sortToNonSort_s1.template get<0>(p) < ghostMarker);
	}

	BOOST_REQUIRE_EQUAL(match,true);

	// Check

	n_out.deviceToHost<0>();

	{
		bool check = true;
		auto it = vPos.getIteratorTo(vPos.size()/2-1);

		while(it.isNext())
		{
			auto p = it.get();

			Point<dim,T> xp = vPos.get(p);

			// Get NN iterator

			auto NN_it = execute_cl_test<impl>::getNN(cl_cpu,cl_cpu.getCell(xp)); /*cl_cpu.getNNIteratorBox(cl_cpu.getCell(xp))*/;

			size_t n_ele = 0;
			while (NN_it.isNext())
			{
				auto q = NN_it.get();

				n_ele++;

				++NN_it;
			}

			check &= n_ele == n_out.template get<0>(p);

			if (check == false)
			{
				std::cout << p << "  " << n_ele << "   " << n_out.template get<0>(p) << "   " << check << std::endl;
				break;
			}

			++it;
		}
		BOOST_REQUIRE_EQUAL(check,true);
	}

	// now we scan the buffer

	openfpm::vector<aggregate<unsigned int>,CudaMemory,memory_traits_inte> n_out_scan;
	openfpm::vector<aggregate<unsigned int>,CudaMemory,memory_traits_inte> nn_list;

	n_out_scan.resize(n_out.size());

	openfpm::scan((unsigned int *)n_out.template getDeviceBuffer<0>(),n_out.size(),(unsigned int *)n_out_scan.template getDeviceBuffer<0>(),gpuContext);

	n_out_scan.template deviceToHost<0>();

	if (n_out_scan.template get<0>(vPos.size()) == 0)
	{return;}

	nn_list.resize(n_out_scan.template get<0>(vPos.size()));

	// Now for each particle we construct the list

	vPos.template hostToDevice<0>();

	execute_cl_test<impl>::calc_list(vPos,sortToNonSort_s1,cl2_split1,n_out_scan,nn_list);
	execute_cl_test<impl>::calc_list_partial(vPos,sortToNonSort_s2,cl2_split2,n_out_scan,n_out_partial,nn_list);

	nn_list.template deviceToHost<0>();

	// Check

	n_out.deviceToHost<0>();

	{
	bool check = true;
	auto it = vPos.getIteratorTo(vPos.size()/2-1);

	while(it.isNext())
	{
		auto p = it.get();

		Point<dim,T> xp = vPos.get(p);

		// Get NN iterator

		openfpm::vector<int> cpu_list;

		auto NN_it = execute_cl_test<impl>::getNN(cl_cpu,cl_cpu.getCell(xp)); /*cl_cpu.getNNIteratorBox(cl_cpu.getCell(xp));*/

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

		// sort both vector

		cpu_list.sort();
		gpu_list.sort();

		for (size_t j = 0 ; j < cpu_list.size() ; j++)
		{check &= cpu_list.get(j) == gpu_list.get(j);}

		if (check == false)
		{
			std::cout << "NPARTS: " << npart << std::endl;

			for (size_t j = 0 ; j < cpu_list.size() ; j++)
			{std::cout << cpu_list.get(j) << "  " << gpu_list.get(j) << std::endl;}

			break;
		}

		++it;
	}

	BOOST_REQUIRE_EQUAL(check,true);

	}
}

BOOST_AUTO_TEST_CASE( CellList_gpu_use_calc_force_box)
{
	std::cout << "Test cell list GPU" << "\n";

	Box<3,float> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Box<3,float> box2({-0.3f,-0.3f,-0.3f},{1.0f,1.0f,1.0f});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,2>(box,1000,{8,8,8});
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,2>(box,10000,{8,8,8});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,2>(box2,1000,{8,8,8});
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,2>(box2,10000,{8,8,8});

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( CellList_gpu_use_calc_force_box_split)
{
	std::cout << "Test cell list GPU split" << "\n";

	Box<3,float> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Box<3,float> box2({-0.3f,-0.3f,-0.3f},{1.0f,1.0f,1.0f});

	Test_cell_gpu_force_split<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,2>(box,1000,{32,32,32});
	Test_cell_gpu_force_split<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,2>(box,10000,{32,32,32});

	Test_cell_gpu_force_split<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,2>(box2,1000,{32,32,32});
	Test_cell_gpu_force_split<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,2>(box2,10000,{32,32,32});

	std::cout << "End cell list GPU split" << "\n";

	// Test the cell list
}

/*BOOST_AUTO_TEST_CASE( CellList_gpu_use_calc_force_box_split_test_performance)
{
	typedef float T;
	constexpr int dim = 3;
	typedef CellList_gpu<3,float,CudaMemory,shift_only<3,float>> CellS;

	std::cout << "Performance" << "\n";

	Box<3,float> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	size_t div[] = {64,64,64};

	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cl2_split1(box,div,2);
	CellS cl2_split2(box,div,2);

	cl2_split1.setBoxNN(2);
	cl2_split2.setBoxNN(2);

	T radius = (box.getHigh(0) - box.getLow(0))/div[0] * 2.0;

	// vector of particles

	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> vPos;
	openfpm::vector<Point<dim,T>,CudaMemory,typename memory_traits_inte<Point<dim,T>>::type,memory_traits_inte> vPosReorder;

	openfpm::vector<aggregate<T,T[3]>,CudaMemory,typename memory_traits_inte<aggregate<T,T[3]>>::type,memory_traits_inte> vPrp;
	openfpm::vector<aggregate<T,T[3]>,CudaMemory,typename memory_traits_inte<aggregate<T,T[3]>>::type,memory_traits_inte> vPrpReorder;

	openfpm::vector<aggregate<unsigned int>,CudaMemory,typename memory_traits_inte<aggregate<unsigned int>>::type,memory_traits_inte> n_out;
	openfpm::vector<aggregate<unsigned int>,CudaMemory,typename memory_traits_inte<aggregate<unsigned int>>::type,memory_traits_inte> n_out_partial;

	// create random particles

	fill_random_parts<3>(box,vPos,vPrp,1400000);

	vPrpReorder.resize(vPos.size());
	vPosReorder.resize(vPos.size());
	n_out.resize(vPos.size()+1);
	n_out.fill<0>(0);

	vPrp.resize(vPos.size());
	vPrpReorder.resize(vPos.size());

	vPos.hostToDevice<0>();
	vPrp.hostToDevice<0,1>();

	size_t ghostMarker = vPos.size() / 2;

	gpu::ofp_context_t gpuContext(gpu::gpu_context_opt::no_print_props);

	cl2_split1.construct(vPos,vPrp,gpuContext,ghostMarker,0,vPos.size()/2);
	cl2_split2.construct(vPos,vPrp,gpuContext,ghostMarker,vPos.size()/2,vPos.size());

	cudaDeviceSynchronize();

	timer t;
	t.start();

	cl2_split1.construct(vPos,vPrp,gpuContext,ghostMarker,0,vPos.size()/2);
	cl2_split2.construct(vPos,vPrp,gpuContext,ghostMarker,vPos.size()/2,vPos.size());

	t.stop();
	std::cout << "Time: " << t.getwct() << std::endl;

	cudaDeviceSynchronize();

	cl2_split1.construct(vPos,vPrp,gpuContext,ghostMarker,0,vPos.size());

	cudaDeviceSynchronize();

	timer t2;
	t2.start();

	cl2_split1.construct(vPos,vPrp,gpuContext,ghostMarker,0,vPos.size());

	t2.stop();
	std::cout << "Time: " << t2.getwct() << std::endl;

	std::cout << "Performance" << "\n";

	// Test the cell list
}*/


BOOST_AUTO_TEST_CASE( CellList_gpu_use_calc_force_box_sparse)
{
	std::cout << "Test cell list GPU" << "\n";

	Box<3,float> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Box<3,float> box2({-0.3f,-0.3f,-0.3f},{1.0f,1.0f,1.0f});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>,true>,2>(box,1000,{32,32,32},2);
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>,true>,2>(box,10000,{32,32,32},2);

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>,true>,2>(box2,1000,{32,32,32},2);
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>,true>,2>(box2,10000,{32,32,32},2);

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( CellList_gpu_use_calc_force_radius)
{
	std::cout << "Test cell list GPU" << "\n";

	Box<3,float> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Box<3,float> box2({-0.3f,-0.3f,-0.3f},{1.0f,1.0f,1.0f});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,1>(box,1000,{32,32,32});
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,1>(box,10000,{32,32,32});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,1>(box2,1000,{32,32,32});
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,1>(box2,10000,{32,32,32});

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

#if 0

BOOST_AUTO_TEST_CASE( CellList_gpu_use_calc_force)
{
	std::cout << "Test cell list GPU" << "\n";

	Box<3,float> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Box<3,float> box2({-0.3f,-0.3f,-0.3f},{1.0f,1.0f,1.0f});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,0>(box,1000,{16,16,16});
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,0>(box,10000,{16,16,16});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,0>(box2,1000,{16,16,16});
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>>,0>(box2,10000,{16,16,16});

	std::cout << "End cell list GPU" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( CellList_gpu_use_calc_force_sparse)
{
	std::cout << "Test cell list GPU force sparse" << "\n";

	Box<3,float> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Box<3,float> box2({-0.3f,-0.3f,-0.3f},{1.0f,1.0f,1.0f});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>,true>,0>(box,1000,{16,16,16});
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>,true>,0>(box,10000,{16,16,16});

	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>,true>,0>(box2,1000,{16,16,16});
	Test_cell_gpu_force<3,float,CellList_gpu<3,float,CudaMemory,shift_only<3,float>,true>,0>(box2,10000,{16,16,16});

	std::cout << "End cell list GPU force sparse" << "\n";

	// Test the cell list
}

#endif

template<typename CellList_type, typename Vector_type, typename Vector_out>
__global__ void cl_offload_gpu(CellList_type cellList, Vector_type parts, Vector_out output)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= parts.size()) return;

    Point<3,float> xp = parts.template get<0>(p);

    output.template get<0>(p) = cellList.getNelements(cellList.getCell(xp));
}

template<typename CellList_type, typename Vector_type, typename Vector_scan_type, typename Vector_list_type>
__global__ void cl_offload_gpu_list(CellList_type cellList, Vector_type parts, Vector_scan_type scan, Vector_list_type list)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= parts.size()) return;

    Point<3,float> xp = parts.template get<0>(p);

    int id = cellList.getCell(xp);
    int n_ele = cellList.getNelements(id);
    int start = scan.template get<0>(p);

    for (int j = 0 ; j < n_ele ; j++)
    {
		list.template get<0>(start+j) = cellList.get(id,j);
    }

}

#if 0

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

	CUDA_LAUNCH_DIM3((cl_offload_gpu<decltype(cl1.toKernel()),decltype(v.toKernel()),decltype(os.toKernel())>),ite.wthr,ite.thr,cl1.toKernel(),v.toKernel(),os.toKernel());

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

	gpu::ofp_context_t gpuContext;
	openfpm::scan((int *)os.template getDeviceBuffer<0>(),os.size(),(int *)os_scan.template getDeviceBuffer<0>(),gpuContext);

	os_scan.deviceToHost<0>();
	os.deviceToHost<0>(os.size()-1,os.size()-1);
	size_t size_list = os_scan.template get<0>(os_scan.size()-1) + os.template get<0>(os.size()-1);

	openfpm::vector_gpu<aggregate<int>> os_list;
	os_list.resize(size_list);

	CUDA_LAUNCH_DIM3((cl_offload_gpu_list<decltype(cl1.toKernel()),decltype(v.toKernel()),
			            decltype(os_scan.toKernel()),decltype(os_list.toKernel())>),ite.wthr,ite.thr,
			            cl1.toKernel(),v.toKernel(),os_scan.toKernel(),os_list.toKernel());

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

#endif

BOOST_AUTO_TEST_CASE( CellList_swap_test )
{
	size_t npart = 4096;

	Box<3,float> box({-1.0,-1.0,-1.0},{1.0,1.0,1.0});

	// Subdivisions
	size_t div[3] = {10,10,10};

	// Origin
	Point<3,float> org({0.0,0.0,0.0});

	// id Cell list
	CellList_gpu<3,float,CudaMemory,shift_only<3,float>> cellList2(box,div,2);
	CellList_gpu<3,float,CudaMemory,shift_only<3,float>> cellList3(box,div,2);
	CellList_gpu<3,float,CudaMemory,shift_only<3,float>> cellList4(box,div,2);

	// vector of particles

	openfpm::vector<Point<3,float>,CudaMemory,memory_traits_inte> vPos;
	openfpm::vector<aggregate<float,float[3]>,CudaMemory,memory_traits_inte> vPrp;

	// create random particles

	fill_random_parts<3>(box,vPos,vPrp,npart);

	vPrp.resize(vPos.size());

	vPos.template hostToDevice<0>();
	vPrp.template hostToDevice<0,1>();

	size_t ghostMarker = vPos.size() / 2;

	gpu::ofp_context_t gpuContext(gpu::gpu_context_opt::no_print_props);
	cellList2.construct(vPos,vPrp,gpuContext,ghostMarker);
	cellList4.construct(vPos,vPrp,gpuContext,ghostMarker);

	cellList3.swap(cellList2);

	// move device to host

	cellList3.debug_deviceToHost();
	cellList4.debug_deviceToHost();

	BOOST_REQUIRE_EQUAL(cellList3.getNCells(),cellList4.getNCells());

	openfpm::vector<size_t> s1;
	openfpm::vector<size_t> s2;

	bool check = true;
	for (size_t i = 0 ; i < cellList3.getNCells() ; i++)
	{
		check &= cellList3.getNelements(i) == cellList4.getNelements(i);

		for (size_t j = 0 ; j < cellList3.getNelements(i) ; j++)
		{
			s1.add(cellList3.get(i,j));
			s2.add(cellList4.get(i,j));
		}

		s1.sort();
		s2.sort();

		for (size_t j = 0 ; j < s1.size() ; j++)
		{
			check &= s1.get(j) == s2.get(j);
		}
	}

	BOOST_REQUIRE_EQUAL(check,true);

	//////////////// We check now that cellList3 and cellList4 match
}

BOOST_AUTO_TEST_SUITE_END()

