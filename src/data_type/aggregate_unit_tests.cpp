/*
 * aggregate_unit_tests.cpp
 *
 *  Created on: Jul 20, 2020
 *      Author: i-bird
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config.h"
#include "data_type/aggregate.hpp"

BOOST_AUTO_TEST_SUITE( aggregate_test )

BOOST_AUTO_TEST_CASE( aggregate_at_c_get_test )
{
	aggregate<double,double[3]> data;

	openfpm::at_c<0>(data) = 1.0;

	openfpm::at_c<1>(data)[0] = 1.0;
	openfpm::at_c<1>(data)[1] = 2.0;
	openfpm::at_c<1>(data)[2] = 3.0;

	BOOST_REQUIRE_EQUAL(openfpm::get<0>(data),1.0);

	BOOST_REQUIRE_EQUAL(openfpm::get<1>(data)[0],1.0);
	BOOST_REQUIRE_EQUAL(openfpm::get<1>(data)[1],2.0);
	BOOST_REQUIRE_EQUAL(openfpm::get<1>(data)[2],3.0);
}

template<unsigned int integer>
struct value_function
{
	enum
	{
		value = integer
	};
};

template<typename arg_f1,typename arg_f2, unsigned int s>
struct function
{
	typedef value_function<arg_f1::value + arg_f2::value + s> value;
};

BOOST_AUTO_TEST_CASE( meta_function_check )
{
	BOOST_REQUIRE_EQUAL((function<value_function<5>,value_function<4>,3>::value::value),12);
}

BOOST_AUTO_TEST_SUITE_END()
