/*
 * Box_unit_tests.hpp
 *
 *  Created on: May 8, 2015
 *      Author: i-bird
 */

#ifndef BOX_UNIT_TESTS_HPP_
#define BOX_UNIT_TESTS_HPP_

BOOST_AUTO_TEST_SUITE( box_test )

BOOST_AUTO_TEST_CASE( box_use)
{
	std::cout << "Box unit test start" << "\n";

	float spacing[2] = {0.1,0.1};

	Box<2,float> sp({1.0,1.0},{2.0,3.0});
	sp.expand(spacing);

	BOOST_REQUIRE_CLOSE(sp.getLow(0),1.0,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getLow(1),1.0,0.0001);

	BOOST_REQUIRE_CLOSE(sp.getHigh(0),2.1,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getHigh(1),3.1,0.0001);

	std::cout << "Box unit test stop" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()



#endif /* BOX_UNIT_TESTS_HPP_ */
