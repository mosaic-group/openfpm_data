#ifndef VECTOR_UNIT_TESTS_HPP
#define VECTOR_UNIT_TESTS_HPP

#include "map_vector.hpp"
#include "Point.hpp"

BOOST_AUTO_TEST_SUITE( vector_test )

// Test the openfpm vector

BOOST_AUTO_TEST_CASE( vector_use)
{
	 openfpm::vector<Point<float>> ofv;

	 // Point
	 Point<float> p;
	 p.setx(1.0);
	 p.sety(2.0);
	 p.setz(3.0);

	 ofv.push_back(p);

	 std::cout << "Vector at position: " << ofv.template get<0>(0) << "\n";
	 exit(-1);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
