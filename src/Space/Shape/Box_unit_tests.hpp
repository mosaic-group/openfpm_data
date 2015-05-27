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

	// expand the box with some spacing

	float spacing[2] = {0.1,0.1};

	Box<2,float> sp({1.0,1.0},{2.0,3.0});
	sp.expand(spacing);

	BOOST_REQUIRE_CLOSE(sp.getLow(0),1.0,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getLow(1),1.0,0.0001);

	BOOST_REQUIRE_CLOSE(sp.getHigh(0),2.1,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getHigh(1),3.1,0.0001);

	// create an enclosing box

	{
	Box<3,float> box1({0.1,0.2,0.3},{1.0,1.1,1.3});
	Box<3,float> box2({0.5,0.6,0.7},{2.0,2.1,2.2});

	box1.enclose(box2);

	BOOST_REQUIRE_EQUAL(box1.getLow(0),0.1f);
	BOOST_REQUIRE_EQUAL(box1.getLow(1),0.2f);
	BOOST_REQUIRE_EQUAL(box1.getLow(2),0.3f);

	BOOST_REQUIRE_EQUAL(box1.getHigh(0),2.0f);
	BOOST_REQUIRE_EQUAL(box1.getHigh(1),2.1f);
	BOOST_REQUIRE_EQUAL(box1.getHigh(2),2.2f);
	}

	// Create the smallest boxes between several boxes or
	// Create a box that is contained into more boxes centering them on p1

	{
	Box<3,float> box1({0.0,0.0,0.0},{1.0,1.1,1.3});
	Box<3,float> box2({0.5,2.0,0.5},{2.0,2.1,2.2});
	Box<3,float> box3({1.5,1.5,4.2},{5.0,5.1,5.2});

	box1.contained(box2);
	box1.contained(box3);

	BOOST_REQUIRE_CLOSE(box1.getHigh(0),1.0f,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1),0.1f,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2),1.0f,0.0001);
	}

	std::cout << "Box unit test stop" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()



#endif /* BOX_UNIT_TESTS_HPP_ */
