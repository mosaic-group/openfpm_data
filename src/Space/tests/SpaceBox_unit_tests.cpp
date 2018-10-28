/*
 * SpaceBox_unit_tests.hpp
 *
 *  Created on: May 8, 2015
 *      Author: i-bird
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


#include "Space/SpaceBox.hpp"

#define N_RANDOM_POINT 1024

BOOST_AUTO_TEST_SUITE( spacebox_test )

BOOST_AUTO_TEST_CASE( spacebox_use)
{
	std::cout << "SpaceBox unit test start" << "\n";

	//! [Definition of a spacebox and rescale]

	float spacing[2] = {0.1,0.1};

	{
	SpaceBox<2,float> sp({1.0,1.0},{2.0,2.0});
	sp.rescale(spacing);

	BOOST_REQUIRE_CLOSE(sp.getLow(0),1.0,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getLow(1),1.0,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getHigh(0),1.1,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getHigh(1),1.1,0.0001);
	}

	//! [Definition of a spacebox and rescale]

	{
	SpaceBox<2,float> sp({1.0,1.0},{2.0,2.0});
	sp.mul(spacing);
	sp.expand(spacing);

	BOOST_REQUIRE_CLOSE(sp.getLow(0),0.1,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getLow(1),0.1,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getHigh(0),0.3,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getHigh(1),0.3,0.0001);
	}

	//! [Definition of a spacebox and intersection between them]

	{
	SpaceBox<2,float> sp1({1.0,1.0},{2.0,2.0});
	SpaceBox<2,float> sp2({0.5,0.5},{1.5,1.5});
	SpaceBox<2,float> sp3;

	bool inte = sp1.Intersect(sp2,sp3);

	BOOST_REQUIRE_EQUAL(inte,true);
	BOOST_REQUIRE_EQUAL(sp3.getLow(0),1.0);
	BOOST_REQUIRE_EQUAL(sp3.getLow(1),1.0);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(0),1.5);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(1),1.5);

	//! [Definition of a spacebox and intersection between them]

	sp1.set({0.0,0.0},{1.0,1.0});
	sp2.set({0.2,-0.5},{0.4,1.5});

	inte = sp1.Intersect(sp2,sp3);

	BOOST_REQUIRE_EQUAL(inte,true);
	BOOST_REQUIRE_EQUAL(sp3.getLow(0),0.2f);
	BOOST_REQUIRE_EQUAL(sp3.getLow(1),0.0f);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(0),0.4f);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(1),1.0f);

	sp1.set({0.0,0.0},{1.0,1.0});
	sp2.set({0.2,0.2},{0.4,0.4});

	inte = sp1.Intersect(sp2,sp3);

	BOOST_REQUIRE_EQUAL(inte,true);
	BOOST_REQUIRE_EQUAL(sp3.getLow(0),0.2f);
	BOOST_REQUIRE_EQUAL(sp3.getLow(1),0.2f);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(0),0.4f);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(1),0.4f);

	sp1.set({0.0,0.0},{1.0,1.0});
	sp2.set({1.2,0.0},{2.4,1.0});

	inte = sp1.Intersect(sp2,sp3);

	BOOST_REQUIRE_EQUAL(inte,false);


	// Box line intersection

	sp1.set({0.0,0.0},{1.0,1.0});
	sp2.set({0.5,0.1},{0.5,1.1});

	inte = sp1.Intersect(sp2,sp3);

	BOOST_REQUIRE_EQUAL(inte,true);
	BOOST_REQUIRE_EQUAL(sp3.getLow(0),0.5f);
	BOOST_REQUIRE_EQUAL(sp3.getLow(1),0.1f);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(0),0.5f);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(1),1.0f);

	// Box Box with line result

	sp1.set({0.0,0.0},{1.0,1.0});
	sp2.set({1.0,0.0},{1.5,1.5});

	inte = sp1.Intersect(sp2,sp3);

	BOOST_REQUIRE_EQUAL(inte,true);
	BOOST_REQUIRE_EQUAL(sp3.getLow(0),1.0f);
	BOOST_REQUIRE_EQUAL(sp3.getLow(1),0.0f);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(0),1.0f);
	BOOST_REQUIRE_EQUAL(sp3.getHigh(1),1.0f);


	sp1.set({0.0,0.0},{1.0,1.0});
	Ghost<2,float> g(0.5);

	sp1.enlarge(g);

	BOOST_REQUIRE_EQUAL(sp1.getLow(0),-0.5f);
	BOOST_REQUIRE_EQUAL(sp1.getLow(1),-0.5f);
	BOOST_REQUIRE_EQUAL(sp1.getHigh(0),1.5f);
	BOOST_REQUIRE_EQUAL(sp1.getHigh(1),1.5f);
	}

	//! [Create random points inside the SpaceBox]

	SpaceBox<3,float> sp_box({0.0,0.0,0.0},{1.0,1.0,1.0});

	for (int i = 0 ; i < N_RANDOM_POINT ; i++)
	{
		Point<3,float> p = sp_box.rnd();

		BOOST_REQUIRE_EQUAL(sp_box.isInside(p),true);
	}

	//! [Create random points inside the SpaceBox]

	// Create random points outside the space box and check

	SpaceBox<3,float> sp_box_out({1.1,1.1,1.1},{2.1,2.1,2.1});

	for (int i = 0 ; i < N_RANDOM_POINT ; i++)
	{
		Point<3,float> p = sp_box_out.rnd();

		BOOST_REQUIRE_EQUAL(sp_box.isInside(p),false);
	}

	std::cout << "SpaceBox unit test stop" << "\n";
}


BOOST_AUTO_TEST_SUITE_END()


