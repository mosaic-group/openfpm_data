#ifndef VECTOR_UNIT_TESTS_HPP
#define VECTOR_UNIT_TESTS_HPP

#include "map_vector.hpp"
#include "Point_test.hpp"

BOOST_AUTO_TEST_SUITE( vector_test )

// Test the openfpm vector

BOOST_AUTO_TEST_CASE( vector_use)
{
	 openfpm::vector<Point_test<float>> ofv;

	 typedef Point_test<float> P;

	 // Point
	 Point_test<float> p;
	 p.setx(1.0);
	 p.sety(2.0);
	 p.setz(3.0);
	 p.sets(4.0);

	 p.get<P::v>()[0] = 5.0;
	 p.get<P::v>()[1] = 6.0;
	 p.get<P::v>()[2] = 7.0;

	 p.get<P::t>()[0][0] = 5.0;
	 p.get<P::t>()[0][1] = 6.0;
	 p.get<P::t>()[0][2] = 7.0;
	 p.get<P::t>()[1][0] = 5.0;
	 p.get<P::t>()[1][1] = 6.0;
	 p.get<P::t>()[1][2] = 7.0;
	 p.get<P::t>()[2][0] = 5.0;
	 p.get<P::t>()[2][1] = 6.0;
	 p.get<P::t>()[2][2] = 7.0;

	 // try to push 1000 object

	 for (size_t i = 0 ; i < 1000 ; i++)
		 ofv.push_back(p);

	 // try to push 2000 object

	 for (size_t i = 0 ; i < 2000 ; i++)
		 ofv.push_back(p);

	 // try to push 3000 object

	 for (size_t i = 0 ; i < 3000 ; i++)
		 ofv.push_back(p);

	 // try to push 10000 object

	 for (size_t i = 0 ; i < 10000 ; i++)
		 ofv.push_back(p);

	 // try to push 20000 object

	 for (size_t i = 0 ; i < 20000 ; i++)
		 ofv.push_back(p);

	 std::cout << "Vector at position: " << ofv.template get<0>(0) << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

#endif
