/*
 * compute_optimal_device_grid_unit_tests.hpp
 *
 *  Created on: Oct 1, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_TEST_COMPUTE_OPTIMAL_DEVICE_GRID_UNIT_TESTS_HPP_
#define OPENFPM_DATA_SRC_UTIL_TEST_COMPUTE_OPTIMAL_DEVICE_GRID_UNIT_TESTS_HPP_

#include "util/compute_optimal_device_grid.hpp"

BOOST_AUTO_TEST_SUITE( compute_optimal_device_grid_test )

BOOST_AUTO_TEST_CASE( compute_optimal_device_use_test )
{
	// Check 1D response

	size_t sz[1] = {1};
	device_grid<1> dg;
	calculate_optimal_device_grid<1>(dg,sz,1,1);

	BOOST_REQUIRE_EQUAL(dg.threads.x,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg.grids.x,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.z,1ul);

	// Fill with garbage
	dg.threads.x = 123;
	dg.threads.y = 123;
	dg.threads.z = 123;

	dg.grids.x = 123;
	dg.grids.y = 123;
	dg.grids.z = 123;

	calculate_optimal_device_grid<1>(dg,sz,8,4);

	BOOST_REQUIRE_EQUAL(dg.threads.x,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg.grids.x,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.z,1ul);

	sz[0] = 6;
	calculate_optimal_device_grid<1>(dg,sz,8,4);

	BOOST_REQUIRE_EQUAL(dg.threads.x,6ul);
	BOOST_REQUIRE_EQUAL(dg.threads.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg.grids.x,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.z,1ul);

	sz[0] = 6;

	// Fill with garbage
	dg.threads.x = 123;
	dg.threads.y = 123;
	dg.threads.z = 123;

	dg.grids.x = 123;
	dg.grids.y = 123;
	dg.grids.z = 123;
	calculate_optimal_device_grid<1>(dg,sz,28,10);

	BOOST_REQUIRE_EQUAL(dg.threads.x,6ul);
	BOOST_REQUIRE_EQUAL(dg.threads.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg.grids.x,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.z,1ul);

	sz[0] = 17;

	calculate_optimal_device_grid<1>(dg,sz,28,10);

	BOOST_REQUIRE_EQUAL(dg.threads.x,17ul);
	BOOST_REQUIRE_EQUAL(dg.threads.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg.grids.x,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.z,1ul);

	sz[0] = 17;

	calculate_optimal_device_grid<1>(dg,sz,1,1);

	BOOST_REQUIRE_EQUAL(dg.threads.x,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg.grids.x,17ul);
	BOOST_REQUIRE_EQUAL(dg.grids.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.z,1ul);

	sz[0] = 17;

	// Fill with garbage
	dg.threads.x = 123;
	dg.threads.y = 123;
	dg.threads.z = 123;

	dg.grids.x = 123;
	dg.grids.y = 123;
	dg.grids.z = 123;
	calculate_optimal_device_grid<1>(dg,sz,5,3);

	BOOST_REQUIRE_EQUAL(dg.threads.x,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg.grids.x,17ul);
	BOOST_REQUIRE_EQUAL(dg.grids.y,1ul);
	BOOST_REQUIRE_EQUAL(dg.grids.z,1ul);


	// 2D test

	size_t sz2[] = {64,64};
	device_grid<2> dg2;

	calculate_optimal_device_grid<2>(dg2,sz2,1536,512);

	BOOST_REQUIRE_EQUAL(dg2.threads.x,64ul);
	BOOST_REQUIRE_EQUAL(dg2.threads.y,16ul);
	BOOST_REQUIRE_EQUAL(dg2.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg2.grids.x,1ul);
	BOOST_REQUIRE_EQUAL(dg2.grids.y,4ul);
	BOOST_REQUIRE_EQUAL(dg2.grids.z,1ul);

	// 2D test

	sz2[0] = 126;
	sz2[1] = 126;

	calculate_optimal_device_grid<2>(dg2,sz2,1536,512);

	BOOST_REQUIRE_EQUAL(dg2.threads.x,18ul);
	BOOST_REQUIRE_EQUAL(dg2.threads.y,18ul);
	BOOST_REQUIRE_EQUAL(dg2.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg2.grids.x,7ul);
	BOOST_REQUIRE_EQUAL(dg2.grids.y,7ul);
	BOOST_REQUIRE_EQUAL(dg2.grids.z,1ul);

	sz2[0] = 162;
	sz2[1] = 126;

	calculate_optimal_device_grid<2>(dg2,sz2,1536,512);

	BOOST_REQUIRE_EQUAL(dg2.threads.x,162ul);
	BOOST_REQUIRE_EQUAL(dg2.threads.y,6ul);
	BOOST_REQUIRE_EQUAL(dg2.threads.z,1ul);

	BOOST_REQUIRE_EQUAL(dg2.grids.x,1ul);
	BOOST_REQUIRE_EQUAL(dg2.grids.y,21ul);
	BOOST_REQUIRE_EQUAL(dg2.grids.z,1ul);

	size_t sz3[3] = {126,126,126};

	device_grid<3> dg3;

	calculate_optimal_device_grid<3>(dg3,sz3,1536,512);

	BOOST_REQUIRE_EQUAL(dg3.threads.x,18ul);
	BOOST_REQUIRE_EQUAL(dg3.threads.y,18ul);
	BOOST_REQUIRE_EQUAL(dg3.threads.z,2ul);

	BOOST_REQUIRE_EQUAL(dg3.grids.x,7ul);
	BOOST_REQUIRE_EQUAL(dg3.grids.y,7ul);
	BOOST_REQUIRE_EQUAL(dg3.grids.z,63ul);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_UTIL_TEST_COMPUTE_OPTIMAL_DEVICE_GRID_UNIT_TESTS_HPP_ */
