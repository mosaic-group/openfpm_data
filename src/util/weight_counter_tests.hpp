/*
 * Weight_counter_tests.hpp
 *
 *  Created on: August 1, 2015
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef WEIGHT_COUNTER_TESTS_HPP_
#define WEIGHT_COUNTER_TESTS_HPP_

#include "Grid/map_grid.hpp"
#include "Space/SpaceBox.hpp"
#include "util/vectors_util.hpp"


BOOST_AUTO_TEST_SUITE( weight_counter_test )

BOOST_AUTO_TEST_CASE( weight_counter_1 )
{
	std::cout << "Weight counter test start" << "\n";

	grid_cpu<2,Point_test<float>> t({5,5});

	{
	SpaceBox<2,float> b({0.5,0.5},{3,3});

	weight_counter<2,grid_cpu<2,Point_test<float>>,float,IS_GRID> w;

	size_t result = w.weight(t,b);

	BOOST_REQUIRE_EQUAL(result,4);
	}

	{
	SpaceBox<2,float> b({0,0},{3,3});

	weight_counter<2,grid_cpu<2,Point_test<float>>,float,IS_GRID> w;

	size_t result = w.weight(t,b);

	BOOST_REQUIRE_EQUAL(result,9);
	}

	std::cout << "Weight counter test stop" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* WEIGHT_COUNTER_TESTS_HPP_ */
