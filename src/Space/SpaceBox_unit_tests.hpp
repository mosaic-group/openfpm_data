/*
 * SpaceBox_unit_tests.hpp
 *
 *  Created on: May 8, 2015
 *      Author: i-bird
 */

#ifndef SPACEBOX_UNIT_TESTS_HPP_
#define SPACEBOX_UNIT_TESTS_HPP_

#include "SpaceBox.hpp"

BOOST_AUTO_TEST_SUITE( spacebox_test )

BOOST_AUTO_TEST_CASE( spacebox_use)
{
	std::cout << "SpaceBox unit test start" << "\n";

	float spacing[2] = {0.1,0.1};

	SpaceBox<2,float> sp({1.0,1.0},{2.0,2.0});
	sp.rescale(spacing);

	BOOST_REQUIRE_CLOSE(sp.getLow(0),1.0,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getLow(1),1.0,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getHigh(0),1.1,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getHigh(1),1.1,0.0001);

	std::cout << "SpaceBox unit test stop" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SPACEBOX_UNIT_TESTS_HPP_ */
