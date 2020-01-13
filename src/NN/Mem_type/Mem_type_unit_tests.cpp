/*
 * Mem_type_unit_tests.cpp
 *
 *  Created on: Dec 26, 2017
 *      Author: Pietro Incardona
 */
#include "config.h"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "NN/Mem_type/MemFast.hpp"
#include "NN/Mem_type/MemBalanced.hpp"
#include "NN/Mem_type/MemMemoryWise.hpp"

BOOST_AUTO_TEST_SUITE( Mem_type_test )

template<typename Mem_type>
void test_mem_type()
{
	Mem_type mem(128);

	//

	mem.init_to_zero(128,10);

	mem.add(0,5);

	BOOST_REQUIRE_EQUAL(mem.getNelements(0),1ul);

	BOOST_REQUIRE_EQUAL(mem.get(0,0),5ul);

	mem.init_to_zero(128,5);

	BOOST_REQUIRE_EQUAL(mem.getNelements(0),0ul);
}

BOOST_AUTO_TEST_CASE ( Mem_type_check )
{
	test_mem_type<Mem_fast<>>();
	test_mem_type<Mem_bal<>>();
	test_mem_type<Mem_mw<>>();
}

BOOST_AUTO_TEST_SUITE_END()
