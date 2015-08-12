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

	//! [assign and operation]

	std::cout << "Point unit test stop" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()




#endif /* SRC_SPACE_SHAPE_POINT_UNIT_TEST_HPP_ */
