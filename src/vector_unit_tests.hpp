#ifndef VECTOR_UNIT_TESTS_HPP
#define VECTOR_UNIT_TESTS_HPP

#include "map_vector.hpp"
#include "Point_test.hpp"

BOOST_AUTO_TEST_SUITE( vector_test )

// Test the openfpm vector

BOOST_AUTO_TEST_CASE( vector_use)
{
	std::cout << "Vector unit test start" << "\n";

	 openfpm::vector<Point_test<float>> ofv;

	 typedef Point_test<float> P;

	 // Point
	 Point_test<float> p;
	 p.setx(1.0);
	 p.sety(2.0);
	 p.setz(3.0);
	 p.sets(4.0);

	 p.get<P::v>()[0] = 1.0;
	 p.get<P::v>()[1] = 2.0;
	 p.get<P::v>()[2] = 7.0;

	 p.get<P::t>()[0][0] = 10.0;
	 p.get<P::t>()[0][1] = 13.0;
	 p.get<P::t>()[0][2] = 8.0;
	 p.get<P::t>()[1][0] = 19.0;
	 p.get<P::t>()[1][1] = 23.0;
	 p.get<P::t>()[1][2] = 5.0;
	 p.get<P::t>()[2][0] = 4.0;
	 p.get<P::t>()[2][1] = 3.0;
	 p.get<P::t>()[2][2] = 11.0;

	 ofv.add(p);

	 // try to push 1000 object

	 for (size_t i = 0 ; i < 1000 ; i++)
	 {
		 // Modify p

		 p.get<P::v>()[0] = 1.0 + i;
		 p.get<P::v>()[1] = 2.0 + i;
		 p.get<P::v>()[2] = 7.0 + i;

		 p.get<P::t>()[0][0] = 10.0 + i;
		 p.get<P::t>()[0][1] = 13.0 + i;
		 p.get<P::t>()[0][2] = 8.0 + i;
		 p.get<P::t>()[1][0] = 19.0 + i;
		 p.get<P::t>()[1][1] = 23.0 + i;
		 p.get<P::t>()[1][2] = 5.0 + i;
		 p.get<P::t>()[2][0] = 4.0 + i;
		 p.get<P::t>()[2][1] = 3.0 + i;
		 p.get<P::t>()[2][2] = 11.0 + i;

		 ofv.add(p);
	 }

	 // try to push 2000 object

	 for (size_t i = 0 ; i < 2000 ; i++)
	 {
		 p.get<P::v>()[0] = 1.0 + i + 1359.0;
		 p.get<P::v>()[1] = 2.0 + i + 1700.0;
		 p.get<P::v>()[2] = 7.0 + i + 1511.0;

		 p.get<P::t>()[0][0] = 10.0 + i + 2333.0;
		 p.get<P::t>()[0][1] = 13.0 + i + 2345.0;
		 p.get<P::t>()[0][2] = 8.0 + i + 1234.0;
		 p.get<P::t>()[1][0] = 19.0 + i + 333.0;
		 p.get<P::t>()[1][1] = 23.0 + i + 4560.0;
		 p.get<P::t>()[1][2] = 5.0 + i + 5672.0;
		 p.get<P::t>()[2][0] = 4.0 + i + 1111.0;
		 p.get<P::t>()[2][1] = 3.0 + i + 3423.0;
		 p.get<P::t>()[2][2] = 11.0 + i + 7321.0;

		 ofv.add(p);
	 }

	 // try to push 3000 object

	 for (size_t i = 0 ; i < 3000 ; i++)
	 {
		 p.get<P::v>()[0] = 1.0 + i + 1569.0;
		 p.get<P::v>()[1] = 2.0 + i + 1720.0;
		 p.get<P::v>()[2] = 7.0 + i + 1211.0;

		 p.get<P::t>()[0][0] = 10.0 + i + 1333.0;
		 p.get<P::t>()[0][1] = 13.0 + i + 1345.0;
		 p.get<P::t>()[0][2] = 8.0 + i + 3234.0;
		 p.get<P::t>()[1][0] = 19.0 + i + 4333.0;
		 p.get<P::t>()[1][1] = 23.0 + i + 7560.0;
		 p.get<P::t>()[1][2] = 5.0 + i + 5872.0;
		 p.get<P::t>()[2][0] = 4.0 + i + 1311.0;
		 p.get<P::t>()[2][1] = 3.0 + i + 3123.0;
		 p.get<P::t>()[2][2] = 11.0 + i + 7821.0;

		 ofv.add(p);
	 }

	 // try to push 10000 object

	 for (size_t i = 0 ; i < 10000 ; i++)
	 {
		 p.get<P::v>()[0] = 1.0 + i + 8569.0;
		 p.get<P::v>()[1] = 2.0 + i + 2720.0;
		 p.get<P::v>()[2] = 7.0 + i + 5211.0;

		 p.get<P::t>()[0][0] = 10.0 + i + 2333.0;
		 p.get<P::t>()[0][1] = 13.0 + i + 4345.0;
		 p.get<P::t>()[0][2] = 8.0 + i + 1234.0;
		 p.get<P::t>()[1][0] = 19.0 + i + 2333.0;
		 p.get<P::t>()[1][1] = 23.0 + i + 5560.0;
		 p.get<P::t>()[1][2] = 5.0 + i + 4872.0;
		 p.get<P::t>()[2][0] = 4.0 + i + 2311.0;
		 p.get<P::t>()[2][1] = 3.0 + i + 1123.0;
		 p.get<P::t>()[2][2] = 11.0 + i + 5821.0;

		 ofv.add(p);
	 }

	 // try to push 20000 object

	 for (size_t i = 0 ; i < 20000 ; i++)
	 {
		 p.get<P::v>()[0] = 1.0 + i + 1569.0;
		 p.get<P::v>()[1] = 2.0 + i + 3720.0;
		 p.get<P::v>()[2] = 7.0 + i + 2211.0;

		 p.get<P::t>()[0][0] = 10.0 + i + 1333.0;
		 p.get<P::t>()[0][1] = 13.0 + i + 2345.0;
		 p.get<P::t>()[0][2] = 8.0 + i + 2234.0;
		 p.get<P::t>()[1][0] = 19.0 + i + 4333.0;
		 p.get<P::t>()[1][1] = 23.0 + i + 1560.0;
		 p.get<P::t>()[1][2] = 5.0 + i + 2872.0;
		 p.get<P::t>()[2][0] = 4.0 + i + 1311.0;
		 p.get<P::t>()[2][1] = 3.0 + i + 5123.0;
		 p.get<P::t>()[2][2] = 11.0 + i + 1821.0;

		 ofv.add(p);
	 }

	 //

	 // Now check the vector

	 BOOST_REQUIRE_EQUAL(ofv.template get<P::x>(0),1.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::y>(0),2.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::z>(0),3.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::s>(0),4.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(0)[0],1.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(0)[1],2.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(0)[2],7.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(0)[0][0],10.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(0)[0][1],13.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(0)[0][2],8.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(0)[1][0],19.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(0)[1][1],23.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(0)[1][2],5.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(0)[2][0],4.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(0)[2][1],3.0);
	 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(0)[2][2],11.0);

	 // Check next 1000 object

	 for (size_t i = 0 ; i < 1000 ; i++)
	 {
		 // Modify p

		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+1)[0],1.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+1)[1],2.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+1)[2],7.0 + i);

		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1)[0][0],10.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1)[0][1],13.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1)[0][2],8.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1)[1][0],19.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1)[1][1],23.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1)[1][2],5.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1)[2][0],4.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1)[2][1],3.0 + i);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1)[2][2],11.0 + i);
	 }

	 // try to push 2000 object

	 for (size_t i = 0 ; i < 2000 ; i++)
	 {
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+1001)[0], 1.0 + i + 1359.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+1001)[1], 2.0 + i + 1700.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+1001)[2], 7.0 + i + 1511.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1001)[0][0], 10.0 + i + 2333.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1001)[0][1], 13.0 + i + 2345.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1001)[0][2], 8.0 + i + 1234.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1001)[1][0], 19.0 + i + 333.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1001)[1][1], 23.0 + i + 4560.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1001)[1][2], 5.0 + i + 5672.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1001)[2][0], 4.0 + i + 1111.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1001)[2][1], 3.0 + i + 3423.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+1001)[2][2], 11.0 + i + 7321.0);
	 }

	 // try to push 3000 object

	 for (size_t i = 0 ; i < 3000 ; i++)
	 {
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+3001)[0], 1.0 + i + 1569.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+3001)[1], 2.0 + i + 1720.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+3001)[2], 7.0 + i + 1211.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+3001)[0][0], 10.0 + i + 1333.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+3001)[0][1], 13.0 + i + 1345.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+3001)[0][2], 8.0 + i + 3234.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+3001)[1][0], 19.0 + i + 4333.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+3001)[1][1], 23.0 + i + 7560.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+3001)[1][2], 5.0 + i + 5872.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+3001)[2][0], 4.0 + i + 1311.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+3001)[2][1], 3.0 + i + 3123.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+3001)[2][2], 11.0 + i + 7821.0);
	 }

	 // try to push 10000 object

	 for (size_t i = 0 ; i < 10000 ; i++)
	 {
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+6001)[0], 1.0 + i + 8569.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+6001)[1], 2.0 + i + 2720.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+6001)[2], 7.0 + i + 5211.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+6001)[0][0], 10.0 + i + 2333.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+6001)[0][1], 13.0 + i + 4345.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+6001)[0][2], 8.0 + i + 1234.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+6001)[1][0], 19.0 + i + 2333.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+6001)[1][1], 23.0 + i + 5560.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+6001)[1][2], 5.0 + i + 4872.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+6001)[2][0], 4.0 + i + 2311.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+6001)[2][1], 3.0 + i + 1123.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+6001)[2][2], 11.0 + i + 5821.0);
	 }

	 // try to push 20000 object

	 for (size_t i = 0 ; i < 20000 ; i++)
	 {
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+16001)[0], 1.0 + i + 1569.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+16001)[1], 2.0 + i + 3720.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::v>(i+16001)[2], 7.0 + i + 2211.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+16001)[0][0], 10.0 + i + 1333.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+16001)[0][1], 13.0 + i + 2345.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+16001)[0][2], 8.0 + i + 2234.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+16001)[1][0], 19.0 + i + 4333.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+16001)[1][1], 23.0 + i + 1560.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+16001)[1][2], 5.0 + i + 2872.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+16001)[2][0], 4.0 + i + 1311.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+16001)[2][1], 3.0 + i + 5123.0);
		 BOOST_REQUIRE_EQUAL(ofv.template get<P::t>(i+16001)[2][2], 11.0 + i + 1821.0);
	 }

	 std::cout << "Vector unit test end" << "\n";

//	 std::cout << "Vector at position: " << ofv.template get<0>(0) << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

#endif
