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
#include "Vector/map_vector.hpp"



BOOST_AUTO_TEST_SUITE( weight_counter_test )

BOOST_AUTO_TEST_CASE( weight_counter_grid )
{
	std::cout << "Grid weight counter test start" << "\n";

	grid_cpu<3,Point<3,float>> t({5,5,5});
	const SpaceBox<3,float> domain({0.0,0.0,0.0},{5.0,5.0,5.0});

	//general case of box, create and check weight
	{
	SpaceBox<3,float> b({0.5,0.5,0.5},{3.5,3.5,3.5});

	weight_counter<3,grid_cpu<3,Point<3,float>>,float,IS_GRID> w;

	size_t result = w.weight(t,domain,b);

	BOOST_REQUIRE_EQUAL(result,27);
	}

	//case of box with integer Low and High points, create and check weight (positive exclude)
	{
	SpaceBox<3,float> b({0.0,0.0,0.0},{3.0,3.0,3.0});

	weight_counter<3,grid_cpu<3,Point<3,float>>,float,IS_GRID> w;

	size_t result = w.weight(t,domain,b);

	BOOST_REQUIRE_EQUAL(result,27);
	}

	std::cout << "Grid weight counter test stop" << "\n";
}

BOOST_AUTO_TEST_CASE( weight_counter_vector )
{
	std::cout << "Vector weight counter test start" << "\n";

	openfpm::vector<Point<3,float>> t;

	SpaceBox<3,float> b({1.0,1.0,1.0},{3.0,3.0,3.0});
	const SpaceBox<3,float> domain({0.0,0.0,0.0},{5.0,5.0,5.0});

	//points to test
	Point<3,float> p = b.rndPE();
	Point<3,float> p1({1.0,1.0,1.0});
	Point<3,float> p2({3.0,3.0,3.0});

	t.add(p);
	t.add(p1);
	t.add(p2);

	weight_counter<3,openfpm::vector<Point<3,float>>,float,IS_VECTOR> w;

	size_t result = w.weight(t,domain,b);

	BOOST_REQUIRE_EQUAL(result,2);

	std::cout << "Vector weight counter test stop" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* WEIGHT_COUNTER_TESTS_HPP_ */
