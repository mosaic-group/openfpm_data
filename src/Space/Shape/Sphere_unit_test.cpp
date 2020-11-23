/*
 * Sphere_unit_test.cpp
 *
 *  Created on: Feb 26, 2020
 *      Author: i-bird
 */


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Grid/map_grid.hpp"
#include "data_type/aggregate.hpp"
#include "Vector/map_vector.hpp"

BOOST_AUTO_TEST_SUITE( sphere_test )


BOOST_AUTO_TEST_CASE( Sphere_test_use)
{
	Point<3,double> p({0.1,0.1,0.1});
	Point<3,double> p1({0.12,0.12,0.12});
	Point<3,double> p3({0.25,0.25,0.25});

	Sphere<3,double> s(p,0.1);

	BOOST_REQUIRE_EQUAL(s.isInside(p1),true);
	BOOST_REQUIRE_EQUAL(s.isInside(p3),false);

	double dist = s.distance(p3);
	BOOST_REQUIRE_EQUAL(dist,0.15980762113533162);
}

BOOST_AUTO_TEST_SUITE_END()
