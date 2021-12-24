/*
 * grid_key_dx_expression_unit_tests.hpp
 *
 *  Created on: Jun 21, 2015
 *      Author: i-bird
 */

#ifndef GRID_KEY_DX_EXPRESSION_UNIT_TESTS_HPP_
#define GRID_KEY_DX_EXPRESSION_UNIT_TESTS_HPP_

#include "Grid/grid_key.hpp"

BOOST_AUTO_TEST_SUITE( grid_expression_test )


BOOST_AUTO_TEST_CASE( grid_expression_use)
{
	const grid_key_dx<3> key1(1,2,3);
	const comb<3> c({(char)1,0,(char)-1});
	const grid_key_dx<3> key2(4,5,6);
	const grid_key_dx<3> key3(12,11,10);

	grid_key_dx<3> res2 = key3 + key2;

	BOOST_REQUIRE_EQUAL(res2.get(0),16);
	BOOST_REQUIRE_EQUAL(res2.get(1),16);
	BOOST_REQUIRE_EQUAL(res2.get(2),16);

	grid_key_dx<3> res = key3 + c - key2;

	BOOST_REQUIRE_EQUAL(res.get(0),3);
	BOOST_REQUIRE_EQUAL(res.get(1),6);
	BOOST_REQUIRE_EQUAL(res.get(2),9);

	res = res - res;

	BOOST_REQUIRE_EQUAL(res.get(0),0);
	BOOST_REQUIRE_EQUAL(res.get(1),0);
	BOOST_REQUIRE_EQUAL(res.get(2),0);

	const Point<3,long int> p({1,2,3});

	res = res + p;

	BOOST_REQUIRE_EQUAL(res.get(0),1);
	BOOST_REQUIRE_EQUAL(res.get(1),2);
	BOOST_REQUIRE_EQUAL(res.get(2),3);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* GRID_KEY_DX_EXPRESSION_UNIT_TESTS_HPP_ */
