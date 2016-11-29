/*
 * grid_iterators_unit_tests.hpp
 *
 *  Created on: Jun 27, 2016
 *      Author: i-bird
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Grid/iterators/grid_skin_iterator.hpp"
#include "Grid/map_grid.hpp"
#include "data_type/aggregate.hpp"

BOOST_AUTO_TEST_SUITE( grid_iterators_tests )

template <unsigned int dim> void test_skin_iterator(Box<3,size_t> & bx1, Box<3,size_t> & bx2, grid_sm<3,void> & g_sm,size_t (& bc)[dim] , size_t vol_test)
{
	grid_cpu<3,aggregate<size_t>> gtest(g_sm.getSize());
	gtest.setMemory();
	auto it = gtest.getSubIterator(0);

	while (it.isNext())
	{
		auto key = it.get();

		gtest.get<0>(key) = 0;

		++it;
	}

	size_t count = 0;
	grid_skin_iterator_bc<3> gsi(g_sm,bx1,bx2,bc);

	while (gsi.isNext() == true)
	{
		auto key = gsi.get();

		gtest.get<0>(key) += 1;

		count++;

		++gsi;
	}

	BOOST_REQUIRE_EQUAL(count,(size_t)vol_test);

	bool ret = true;
	auto it2 = gtest.getSubIterator(0);

	while (it2.isNext())
	{
		auto key = it2.get();

		ret &= gtest.get<0>(key) <= 1;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(ret,true);
}

BOOST_AUTO_TEST_CASE( grid_skin_iterator_test )
{
	// Test box inside box

	size_t sz[] = {112,112,112};
	Box<3,size_t> bx1({20,25,30},{90,93,98});
	Box<3,size_t> bx2({18,23,28},{92,95,100});
	Box<3,size_t> bx3({20,25,30},{90,93,30});
	Box<3,size_t> bx4({20,25,30},{90,93,50});
	Box<3,size_t> bx5({19,24,29},{90,93,50});

	grid_sm<3,void> g_sm(sz);
	size_t bc[] = {PERIODIC,PERIODIC,PERIODIC};

	// Volume calculation

	size_t tot = 1;

	tot *= bx2.getHigh(0) - bx2.getLow(0) + 1;
	tot *= bx2.getHigh(1) - bx2.getLow(1) + 1;
	tot *= bx2.getHigh(2) - bx2.getLow(2) + 1;

	size_t tot_ = 1;

	tot_ *= bx1.getHigh(0) - bx1.getLow(0) + 1 - 2;
	tot_ *= bx1.getHigh(1) - bx1.getLow(1) + 1 - 2;
	tot_ *= bx1.getHigh(2) - bx1.getLow(2) + 1 - 2;

	tot -= tot_;

	test_skin_iterator<3>(bx1,bx2,g_sm,bc,tot);
	test_skin_iterator<3>(bx2,bx1,g_sm,bc,0);
	test_skin_iterator<3>(bx3,bx3,g_sm,bc,4899);

	Point<3,size_t> p({60,60,60});

	bx2 += p;
	bx1 += p;

	test_skin_iterator<3>(bx1,bx2,g_sm,bc,tot);
	test_skin_iterator<3>(bx2,bx1,g_sm,bc,0);

	test_skin_iterator<3>(bx4,bx4,g_sm,bc,15042);
	test_skin_iterator<3>(bx5,bx4,g_sm,bc,7679);
}

BOOST_AUTO_TEST_SUITE_END()

