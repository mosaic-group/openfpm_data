#ifndef VECTOR_UNIT_TESTS_HPP
#define VECTOR_UNIT_TESTS_HPP

#include "map_vector.hpp"
#include "Point.hpp"

BOOST_AUTO_TEST_SUITE( vector_test )

BOOST_AUTO_TEST_CASE( vector_use)
{
	 openfpm::vector<Point<float>> ofv;
}

BOOST_AUTO_TEST_SUITE_END()

#endif
