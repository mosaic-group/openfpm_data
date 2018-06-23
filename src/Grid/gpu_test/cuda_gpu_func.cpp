/*
 * cuda_gpu_func.cpp
 *
 *  Created on: Jun 3, 2018
 *      Author: i-bird
 */

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Grid/map_grid.hpp"
#include "Point_test.hpp"

BOOST_AUTO_TEST_SUITE( grid_gpu_func_test )


BOOST_AUTO_TEST_CASE (gpu_computation_func)
{
#ifdef CUDA_GPU

	size_t sz[3] = {64,64,64};
	grid_gpu<3, Point_test<float> > c3(sz);

	grid_key_dx<3> k1({1,1,1});
	grid_key_dx<3> k2({62,62,62});

	c3.setMemory();

	auto gcf = c3.getGPUIterator(k1,k2);

	BOOST_REQUIRE_EQUAL(gcf.thr.x,16ul);
	BOOST_REQUIRE_EQUAL(gcf.thr.y,8ul);
	BOOST_REQUIRE_EQUAL(gcf.thr.z,8ul);

	BOOST_REQUIRE_EQUAL(gcf.wthr.x,4ul);
	BOOST_REQUIRE_EQUAL(gcf.wthr.y,8ul);
	BOOST_REQUIRE_EQUAL(gcf.wthr.z,8ul);

	grid_key_dx<3> k3({50,50,50});
	grid_key_dx<3> k4({62,62,62});
	grid_key_dx<3> k5({60,61,62});

	auto gcf2 = c3.getGPUIterator(k3,k4);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,13ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,8ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,8ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,2ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,2ul);

	gcf2 = c3.getGPUIterator(k3,k4,511);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,8ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,8ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,4ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,2ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,2ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,4ul);

	gcf2 = c3.getGPUIterator(k3,k4,1);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,1ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,13ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,13ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,13ul);

	gcf2 = c3.getGPUIterator(k3,k5,32);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,4ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,4ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,2ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,3ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,3ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,7ul);

	gcf2 = c3.getGPUIterator(k3,k5,1);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,1ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,11ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,12ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,13ul);

#endif
}


BOOST_AUTO_TEST_SUITE_END()
