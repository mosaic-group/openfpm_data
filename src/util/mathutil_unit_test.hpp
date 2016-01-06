/*
 * mathutil_unit_test.hpp
 *
 *  Created on: Dec 23, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_MATHUTIL_UNIT_TEST_HPP_
#define OPENFPM_DATA_SRC_UTIL_MATHUTIL_UNIT_TEST_HPP_

#include "mathutil.hpp"

BOOST_AUTO_TEST_SUITE( mathutil_test )

BOOST_AUTO_TEST_CASE( math )
{
	//! [factorial usage]

	size_t f = openfpm::math::factorial(0);
	BOOST_REQUIRE_EQUAL(f,1ul);

	f = openfpm::math::factorial(1);
	BOOST_REQUIRE_EQUAL(f,1ul);

	f = openfpm::math::factorial(7);
	BOOST_REQUIRE_EQUAL(f,5040ul);

	//! [factorial usage]

	//! [Combination usage]

	f = openfpm::math::C(7,0);
	BOOST_REQUIRE_EQUAL(f,1ul);

	f = openfpm::math::C(7,7);
	BOOST_REQUIRE_EQUAL(f,1ul);

	f = openfpm::math::C(7,2);
	BOOST_REQUIRE_EQUAL(f,21ul);

	//! [Combination usage]

	//! [round to big pow]

	f = openfpm::math::round_big_2(3);
	BOOST_REQUIRE_EQUAL(f,4ul);

	f = openfpm::math::round_big_2(7);
	BOOST_REQUIRE_EQUAL(f,8ul);

	f = openfpm::math::round_big_2(21);
	BOOST_REQUIRE_EQUAL(f,32ul);

	//! [round to big pow]

	//! [pow]

	f = openfpm::math::pow(3,3);
	BOOST_REQUIRE_EQUAL(f,27ul);

	f = openfpm::math::pow(4,2);
	BOOST_REQUIRE_EQUAL(f,16ul);

	f = openfpm::math::pow(2,5);
	BOOST_REQUIRE_EQUAL(f,32ul);

	//! [pow]

	//! [positive modulo]

	f = openfpm::math::positive_modulo(-1,3);
	BOOST_REQUIRE_EQUAL(f,2ul);

	f = openfpm::math::positive_modulo(-5,2);
	BOOST_REQUIRE_EQUAL(f,1ul);

	f = openfpm::math::positive_modulo(-7,5);
	BOOST_REQUIRE_EQUAL(f,3ul);

	//! [positive modulo]

	//! [periodic]

	float ff = openfpm::math::periodic(10.9,1.1,0.1);
	BOOST_REQUIRE_CLOSE(ff,0.9,0.0001);

	ff = openfpm::math::periodic(0.0,1.1,0.1);
	BOOST_REQUIRE_CLOSE(ff,1.0,0.0001);

	ff = openfpm::math::periodic(-7.0,1.1,0.1);
	BOOST_REQUIRE_CLOSE(ff,1.0,0.0001);

	// in float representation 1.1 is a little bigger than 1.1
	ff = openfpm::math::periodic(-10.9,1.1,0.1);
	BOOST_REQUIRE_CLOSE(ff,1.1,0.0001);

	//! [periodic]

	//! [periodic_l]

	ff = openfpm::math::periodic(0.0,1.1,0.1);
	BOOST_REQUIRE_CLOSE(ff,1.0,0.0001);

	ff = openfpm::math::periodic(1.2,1.1,0.1);
	BOOST_REQUIRE_CLOSE(ff,0.2,0.0001);

	// in float representation 1.1 is a little bigger than 1.1
	ff = openfpm::math::periodic(0.9,1.1,0.1);
	BOOST_REQUIRE_CLOSE(ff,0.9,0.0001);

	//! [periodic_l]
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_UTIL_MATHUTIL_UNIT_TEST_HPP_ */
