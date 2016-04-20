/*
 * Point_unit_test.hpp
 *
 *  Created on: Aug 11, 2015
 *      Author: i-bird
 */

#ifndef SRC_SPACE_SHAPE_POINT_UNIT_TEST_HPP_
#define SRC_SPACE_SHAPE_POINT_UNIT_TEST_HPP_

#include <iostream>
#include "Space/Shape/Point.hpp"
#include "Grid/comb.hpp"
#include "util/convert.hpp"

BOOST_AUTO_TEST_SUITE( Point_test_suite )

BOOST_AUTO_TEST_CASE( Point_use )
{
	std::cout << "Point unit test start" << "\n";

	//! [assign and operation]

	Point<3,float> p;

	p.get(0) = 1.0;
	p.get(1) = 2.0;
	p.get(2) = 3.0;

	p = p*p + p;

	BOOST_REQUIRE_EQUAL(p.get(0),2.0);
	BOOST_REQUIRE_EQUAL(p.get(1),6.0);
	BOOST_REQUIRE_EQUAL(p.get(2),12.0);

	// Combination 3D
	comb<3> c;
	c.c[0] = 1;
	c.c[1] = 1;
	c.c[2] = 1;
	Point<3,float> pc = toPoint<3,float>::convert(c);

	Point<3,float> one;
	// fill the point with one
	one.one();
	Point<3,float> middle = p / 2;

	// middle * (one + comb_p) = 0
	p = p + middle * (one + pc);

	BOOST_REQUIRE_EQUAL(p.get(0),4.0);
	BOOST_REQUIRE_EQUAL(p.get(1),12.0);
	BOOST_REQUIRE_EQUAL(p.get(2),24.0);

	//! [assign and operation]

	//! [norm]

	Point<3,float> pn({1.0,1.0,1.0});

	BOOST_REQUIRE_CLOSE(pn.norm(),1.732050,0.0001);

	//! [norm]

	// distance between two points

	Point<3,float> p1({1.0,1.0,1.0});
	Point<3,float> p2({2.0,3.0,4.0});

	BOOST_REQUIRE_EQUAL(p1.distance(p2),p2.distance(p1));
	BOOST_REQUIRE_CLOSE(p1.distance(p2),3.74165738677,0.0001);

	BOOST_REQUIRE_EQUAL(p1.distance2(p2),p2.distance2(p1));
	BOOST_REQUIRE_CLOSE(p1.distance2(p2),14.0,0.0001);

	std::cout << "Point unit test stop" << "\n";

	// check two point are equal

	BOOST_REQUIRE( p1 != p2 );
	BOOST_REQUIRE( p1 == p1 );

	// Check division and operator*

	p2 /= 2.0;
	p1 = p1 * 2.0;

	BOOST_REQUIRE_CLOSE(p2.get(0),1.0,0.001);
	BOOST_REQUIRE_CLOSE(p2.get(1),1.5,0.001);
	BOOST_REQUIRE_CLOSE(p2.get(2),2.0,0.001);

	BOOST_REQUIRE_CLOSE(p1.get(0),2.0,0.001);
	BOOST_REQUIRE_CLOSE(p1.get(1),2.0,0.001);
	BOOST_REQUIRE_CLOSE(p1.get(2),2.0,0.001);


}

BOOST_AUTO_TEST_SUITE_END()




#endif /* SRC_SPACE_SHAPE_POINT_UNIT_TEST_HPP_ */
