/*
 * SparseGrid_unit_tests.cpp
 *
 *  Created on: Oct 22, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_UNIT_TESTS_CPP_
#define OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_UNIT_TESTS_CPP_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "SparseGrid/SparseGrid.hpp"

BOOST_AUTO_TEST_SUITE( sparse_grid_test )

BOOST_AUTO_TEST_CASE( sparse_grid_use_test)
{
	size_t sz[3] = {10000,10000,10000};

	sgrid_cpu<3,aggregate<float>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	// We fill a sphere with a band

	grid_key_dx<3> key1({5000,5000,5000});
	grid_key_dx<3> key2({5001,5001,5001});
	grid_key_dx<3> key3({5002,5003,5003});

	grid.template insert<0>(key1) = 1.0;
	grid.template insert<0>(key2) = 2.0;
	grid.template insert<0>(key3) = 3.0;

	BOOST_REQUIRE_EQUAL(grid.template get<0>(key1),1.0);
	BOOST_REQUIRE_EQUAL(grid.template get<0>(key2),2.0);
	BOOST_REQUIRE_EQUAL(grid.template get<0>(key3),3.0);

	auto it = grid.getDomainIterator();

	size_t count = 0;

	while (it.isNext())
	{
		auto key = it.getKey();

		count++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(count,3);
}

BOOST_AUTO_TEST_CASE( sparse_grid_fill_all_test)
{
	size_t sz[3] = {171,171,171};

	sgrid_cpu<3,aggregate<float>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	grid_sm<3,void> g_sm(sz);

	grid_key_dx_iterator<3> kit(g_sm);

	while (kit.isNext())
	{
		auto key = kit.get();

		grid.template insert<0>(key) = g_sm.LinId(key);

		++kit;
	}

	auto it = grid.getDomainIterator();

	size_t count = 0;

	bool match = true;

	while (it.isNext())
	{
		auto key = it.getKey();

		// return a grid_key_dx
		auto key_pos = it.getKeyF();

		match &= (grid.template get<0>(key) == g_sm.LinId(key_pos));

		count++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(count,171*171*171);
	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_UNIT_TESTS_CPP_ */
