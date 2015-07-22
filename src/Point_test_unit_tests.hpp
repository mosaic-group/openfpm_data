/*
 * Point_test_unit_tests.hpp
 *
 *  Created on: Jun 20, 2015
 *      Author: i-bird
 */

#ifndef POINT_TEST_UNIT_TESTS_HPP_
#define POINT_TEST_UNIT_TESTS_HPP_

#include "Point_test.hpp"

BOOST_AUTO_TEST_SUITE( Point_test_unit_tests )

BOOST_AUTO_TEST_CASE( Point_test_unit_tests )
{
	typedef Point_test<float> P;

	Point_test<float> p;

	// fill the point p with data

	p.setx(1.0);
	p.sety(2.0);
	p.setz(3.0);
	p.sets(4.0);

	p.template get<P::v>()[0] = 5.0;
	p.template get<P::v>()[1] = 6.0;
	p.template get<P::v>()[2] = 7.0;

	p.template get<P::t>()[0][0] = 8.0;
	p.template get<P::t>()[0][1] = 9.0;
	p.template get<P::t>()[0][2] = 10.0;
	p.template get<P::t>()[1][0] = 11.0;
	p.template get<P::t>()[1][1] = 12.0;
	p.template get<P::t>()[1][2] = 13.0;
	p.template get<P::t>()[2][0] = 14.0;
	p.template get<P::t>()[2][1] = 15.0;
	p.template get<P::t>()[2][2] = 16.0;

	BOOST_REQUIRE_EQUAL(p.template get<P::x>(),1.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::y>(),2.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::z>(),3.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::s>(),4.0);

	BOOST_REQUIRE_EQUAL(p.template get<P::v>()[0],5.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::v>()[1],6.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::v>()[2],7.0);

	BOOST_REQUIRE_EQUAL(p.template get<P::t>()[0][0],8.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::t>()[0][1],9.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::t>()[0][2],10.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::t>()[1][0],11.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::t>()[1][1],12.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::t>()[1][2],13.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::t>()[2][0],14.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::t>()[2][1],15.0);
	BOOST_REQUIRE_EQUAL(p.template get<P::t>()[2][2],16.0);

	// operator equal
	Point_test<float> p2 = p;

	BOOST_REQUIRE_EQUAL(p2.template get<P::x>(),1.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::y>(),2.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::z>(),3.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::s>(),4.0);

	BOOST_REQUIRE_EQUAL(p2.template get<P::v>()[0],5.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::v>()[1],6.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::v>()[2],7.0);

	BOOST_REQUIRE_EQUAL(p2.template get<P::t>()[0][0],8.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::t>()[0][1],9.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::t>()[0][2],10.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::t>()[1][0],11.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::t>()[1][1],12.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::t>()[1][2],13.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::t>()[2][0],14.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::t>()[2][1],15.0);
	BOOST_REQUIRE_EQUAL(p2.template get<P::t>()[2][2],16.0);

	// equal from vector

/*	openfpm::vector<Point_test<float>> v;

	v.add(p);

	Point_test<float> p3;

	p3 = v.get(0);

	BOOST_REQUIRE_EQUAL(p3.template get<P::x>(),1.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::y>(),2.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::z>(),3.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::s>(),4.0);

	BOOST_REQUIRE_EQUAL(p3.template get<P::v>()[0],5.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::v>()[1],6.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::v>()[2],7.0);

	BOOST_REQUIRE_EQUAL(p3.template get<P::t>()[0][0],8.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::t>()[0][1],9.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::t>()[0][2],10.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::t>()[1][0],11.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::t>()[1][1],12.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::t>()[1][2],13.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::t>()[2][0],14.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::t>()[2][1],15.0);
	BOOST_REQUIRE_EQUAL(p3.template get<P::t>()[2][2],16.0);*/
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* POINT_TEST_UNIT_TESTS_HPP_ */
