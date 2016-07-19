/*
 * map_vector_std_util_unit_test.hpp
 *
 *  Created on: May 14, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_VECTOR_MAP_VECTOR_STD_UTIL_UNIT_TEST_HPP_
#define OPENFPM_DATA_SRC_VECTOR_MAP_VECTOR_STD_UTIL_UNIT_TEST_HPP_

#include "map_vector_std_util.hpp"
#include "map_vector_std.hpp"

BOOST_AUTO_TEST_SUITE( util_test )

BOOST_AUTO_TEST_CASE( map_vector_std_util )
{
	//! [Check has_base_to_copy]

	bool ret = has_base_to_copy< openfpm::vector<float> >::value;

	BOOST_REQUIRE_EQUAL(ret, true);

	ret = has_base_to_copy< openfpm::vector<float,PtrMemory> >::value;

	BOOST_REQUIRE_EQUAL(ret, false);

	//! [Check has_base_to_copy]
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_VECTOR_MAP_VECTOR_STD_UTIL_UNIT_TEST_HPP_ */
