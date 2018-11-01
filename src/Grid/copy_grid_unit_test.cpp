/*
 * copy_grid_unit_test.hpp
 *
 *  Created on: Dec 4, 2017
 *      Author: i-bird
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Grid/map_grid.hpp"
#include "data_type/aggregate.hpp"
#include "Vector/map_vector.hpp"

BOOST_AUTO_TEST_SUITE( copy_grid_test )

template<unsigned int dim, typename grid>
void Test_copy_grid(grid & g_src, grid & g_dst,
					Box<dim,size_t> & bsrc_1, Box<dim,size_t> & bdst_1)
{
	auto gs1 = g_src.getGrid();
	auto gd1 = g_dst.getGrid();
	auto it = g_src.getIterator();

	grid_key_dx<dim> zero[1];
	zero[0].zero();

	while (it.isNext())
	{
		auto key = it.get();

		g_src.template get<0>(key) = gs1.LinId(key);

		++it;
	}

	copy_grid_fast<false,
				   dim,
				   grid_cpu<dim,aggregate<double>>,
				   grid_sm<dim,aggregate<double>>>::copy(gs1,gd1,
						                  bsrc_1,bdst_1,
										  g_src,g_dst,
										  zero);

	// Check

	bool match = true;

	grid_key_dx_iterator_sub<dim, no_stencil> its(gs1,bsrc_1.getKP1(), bsrc_1.getKP2());
	grid_key_dx_iterator_sub<dim, no_stencil> itd(gd1,bdst_1.getKP1(), bdst_1.getKP2());

	while (its.isNext())
	{
		auto key_s = its.get();
		auto key_d = itd.get();

		match &= g_src.template get<0>(key_s) == g_dst.template get<0>(key_d);

		++its;
		++itd;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

template<unsigned int dim, typename grid>
void Test_copy_grid_cmp(grid & g_src, grid & g_dst,
					Box<dim,size_t> & bsrc_1, Box<dim,size_t> & bdst_1)
{
	auto gs1 = g_src.getGrid();
	auto gd1 = g_dst.getGrid();
	auto it = g_src.getIterator();

	grid_key_dx<dim> zero[1];
	zero[0].zero();

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t k = 0 ; k < (size_t)gs1.LinId(key) % 4 ; k++)
		{
			g_src.template get<0>(key).add(gs1.LinId(key) + 1);
		}

		++it;
	}

	copy_grid_fast<true,
				   dim,
				   grid_cpu<dim,aggregate<openfpm::vector<double>>>,
				   grid_sm<dim,aggregate<openfpm::vector<double>>>>::copy(gs1,gd1,
						                  bsrc_1,bdst_1,
										  g_src,g_dst,
										  zero);

	// Check

	bool match = true;

	grid_key_dx_iterator_sub<dim, no_stencil> its(gs1,bsrc_1.getKP1(), bsrc_1.getKP2());
	grid_key_dx_iterator_sub<dim, no_stencil> itd(gd1,bdst_1.getKP1(), bdst_1.getKP2());

	while (its.isNext())
	{
		auto key_s = its.get();
		auto key_d = itd.get();

		match &= g_src.template get<0>(key_s).size() == g_dst.template get<0>(key_d).size();

		for (size_t i = 0 ; i < g_dst.template get<0>(key_d).size() ; i++)
		{
			match &= g_src.template get<0>(key_s).get(i) == g_dst.template get<0>(key_d).get(i);
		}

		++its;
		++itd;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( copy_grid_test_use)
{
	{
	size_t sz1[1] = {37};
	size_t sz2[2] = {37,37};
	size_t sz3[3] = {37,37,37};
	size_t sz4[4] = {37,37,37,37};

	grid_cpu<1,aggregate<double>> g1_src(sz1);
	grid_cpu<1,aggregate<double>> g1_dst(sz1);
	grid_cpu<2,aggregate<double>> g2_src(sz2);
	grid_cpu<2,aggregate<double>> g2_dst(sz2);
	grid_cpu<3,aggregate<double>> g3_src(sz3);
	grid_cpu<3,aggregate<double>> g3_dst(sz3);
	grid_cpu<4,aggregate<double>> g4_src(sz4);
	grid_cpu<4,aggregate<double>> g4_dst(sz4);
	g1_src.setMemory();
	g1_dst.setMemory();
	g2_src.setMemory();
	g2_dst.setMemory();
	g3_src.setMemory();
	g3_dst.setMemory();
	g4_src.setMemory();
	g4_dst.setMemory();

	// fill all grids

	Box<1,size_t> bsrc_1({4},{11});
	Box<1,size_t> bdst_1{{20},{27}};

	Test_copy_grid(g1_src,g1_dst,bsrc_1,bdst_1);

	Box<2,size_t> bsrc_2({4,7},{11,20});
	Box<2,size_t> bdst_2({20,5},{27,18});

	Test_copy_grid(g2_src,g2_dst,bsrc_2,bdst_2);

	Box<3,size_t> bsrc_3({4,7,1},{11,20,6});
	Box<3,size_t> bdst_3({20,5,10},{27,18,15});

	Test_copy_grid(g3_src,g3_dst,bsrc_3,bdst_3);

#ifdef SE_CLASS1
	Box<4,size_t> bsrc_4({4,7,1,3},{11,20,6,6});
	Box<4,size_t> bdst_4({20,5,10,13},{27,18,15,16});

	Test_copy_grid(g4_src,g4_dst,bsrc_4,bdst_4);
#else
	Box<4,size_t> bsrc_4({4,7,1,3},{11,20,6,6});
	Box<4,size_t> bdst_4({20,5,10,13},{27,18,15,16});

	Test_copy_grid(g4_src,g4_dst,bsrc_4,bdst_4);
#endif
	}
	///////////

	{
	size_t sz1[1] = {37};
	size_t sz2[2] = {37,37};
	size_t sz3[3] = {37,37,37};
	size_t sz4[4] = {37,37,37,37};

	grid_cpu<1,aggregate<openfpm::vector<double>>> g1_src(sz1);
	grid_cpu<1,aggregate<openfpm::vector<double>>> g1_dst(sz1);
	grid_cpu<2,aggregate<openfpm::vector<double>>> g2_src(sz2);
	grid_cpu<2,aggregate<openfpm::vector<double>>> g2_dst(sz2);
	grid_cpu<3,aggregate<openfpm::vector<double>>> g3_src(sz3);
	grid_cpu<3,aggregate<openfpm::vector<double>>> g3_dst(sz3);
	grid_cpu<4,aggregate<openfpm::vector<double>>> g4_src(sz4);
	grid_cpu<4,aggregate<openfpm::vector<double>>> g4_dst(sz4);
	g1_src.setMemory();
	g1_dst.setMemory();
	g2_src.setMemory();
	g2_dst.setMemory();
	g3_src.setMemory();
	g3_dst.setMemory();
	g4_src.setMemory();
	g4_dst.setMemory();

	// fill all grids

	Box<1,size_t> bsrc_1({4},{11});
	Box<1,size_t> bdst_1{{20},{27}};

	Test_copy_grid_cmp(g1_src,g1_dst,bsrc_1,bdst_1);

	Box<2,size_t> bsrc_2({4,7},{11,20});
	Box<2,size_t> bdst_2({20,5},{27,18});

	Test_copy_grid_cmp(g2_src,g2_dst,bsrc_2,bdst_2);

#ifndef SE_CLASS2

	Box<3,size_t> bsrc_3({4,7,1},{11,20,6});
	Box<3,size_t> bdst_3({20,5,10},{27,18,15});

	Test_copy_grid_cmp(g3_src,g3_dst,bsrc_3,bdst_3);

#ifdef SE_CLASS1
	Box<4,size_t> bsrc_4({4,7,1,3},{7,14,6,6});
	Box<4,size_t> bdst_4({20,5,10,13},{23,12,15,16});

	Test_copy_grid_cmp(g4_src,g4_dst,bsrc_4,bdst_4);
#else
	Box<4,size_t> bsrc_4({4,7,1,3},{11,20,6,6});
	Box<4,size_t> bdst_4({20,5,10,13},{27,18,15,16});

	Test_copy_grid_cmp(g4_src,g4_dst,bsrc_4,bdst_4);
#endif
#endif
	}
}

BOOST_AUTO_TEST_SUITE_END()


