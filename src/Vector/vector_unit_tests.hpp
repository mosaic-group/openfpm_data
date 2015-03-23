#ifndef VECTOR_UNIT_TESTS_HPP
#define VECTOR_UNIT_TESTS_HPP

#include "map_vector.hpp"
#include "Point_test.hpp"

#define FIRST_PUSH 2000000
#define SECOND_PUSH 2000000

BOOST_AUTO_TEST_SUITE( vector_test )

// Test the openfpm vector

BOOST_AUTO_TEST_CASE( vector_use)
{
	std::cout << "Vector unit test start" << "\n";

	 std::vector<Point_orig<float>> v_stl_test;
	 openfpm::vector<Point_test<float>> v_ofp_test;

//	 v_stl_test.reserve(FIRST_PUSH);
//	 v_ofp_test.reserve(FIRST_PUSH);

	 typedef Point_test<float> P;

	 // Now fill the STL vector

	 //#ifdef VERBOSE_TEST
	 	  std::cout << "Vector test: " << "\n";

	 	   timespec ts_start;
	 	   // clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
	 	   clock_gettime(CLOCK_REALTIME, &ts_start); // Works on Linux
	 //#endif

	   clock_gettime(CLOCK_REALTIME, &ts_start);

		 Point_orig<float> po;
		 po.setx(1.0);
		 po.sety(2.0);
		 po.setz(3.0);
		 po.sets(4.0);

		 for (size_t i = 0 ; i < FIRST_PUSH ; i++)
		 {
			 // Modify p

			 po.v[0] = 1.0 + i;
			 po.v[1] = 2.0 + i;
			 po.v[2] = 7.0 + i;

			 po.t[0][0] = 10.0 + i;
			 po.t[0][1] = 13.0 + i;
			 po.t[0][2] = 8.0 + i;
			 po.t[1][0] = 19.0 + i;
			 po.t[1][1] = 23.0 + i;
			 po.t[1][2] = 5.0 + i;
			 po.t[2][0] = 4.0 + i;
			 po.t[2][1] = 3.0 + i;
			 po.t[2][2] = 11.0 + i;

			 // add p

			 v_stl_test.push_back(po);
		 }

		 //#ifdef VERBOSE_TEST
		 	  timespec end_time;
		 	  clock_gettime(CLOCK_REALTIME, &end_time); // Works on Linux
		 	   float time_dif =(float)( end_time.tv_sec - ts_start.tv_sec  + (double)(end_time.tv_nsec - ts_start.tv_nsec)/1000000000.0 );

		 	   std::cout << "OPENFPM STL : " << FIRST_PUSH << " add " << "  Time: " << time_dif << " s  " << "\n";
		 //#endif

		 clock_gettime(CLOCK_REALTIME, &ts_start);

	 // Point
	 Point_test<float> p;
	 p.setx(1.0);
	 p.sety(2.0);
	 p.setz(3.0);
	 p.sets(4.0);

	 // push objects

	 for (size_t i = 0 ; i < FIRST_PUSH ; i++)
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

		 // add p

		 v_ofp_test.add(p);
	 }

//#ifdef VERBOSE_TEST
	  end_time;
	  clock_gettime(CLOCK_REALTIME, &end_time); // Works on Linux
	   time_dif =(float)( end_time.tv_sec - ts_start.tv_sec  + (double)(end_time.tv_nsec - ts_start.tv_nsec)/1000000000.0 );

	   std::cout << "OPENFPM End : " << FIRST_PUSH << " add " << "  Time: " << time_dif << " s  " << "\n";
//#endif

	 // try to duplicate the vector

	 openfpm::vector<Point_test<float>> dv_ofp_test = v_ofp_test.duplicate();

	 // Check if the duplicated vector match

	 for (size_t i = 0 ; i < FIRST_PUSH ; i++)
	 {
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::v>(i)[0],v_ofp_test.template get<P::v>(i)[0]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::v>(i)[1],v_ofp_test.template get<P::v>(i)[1]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::v>(i)[2],v_ofp_test.template get<P::v>(i)[2]);

		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::t>(i)[0][0],v_ofp_test.template get<P::t>(i)[0][0]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::t>(i)[0][1],v_ofp_test.template get<P::t>(i)[0][1]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::t>(i)[0][2],v_ofp_test.template get<P::t>(i)[0][2]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::t>(i)[1][0],v_ofp_test.template get<P::t>(i)[1][0]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::t>(i)[1][1],v_ofp_test.template get<P::t>(i)[1][1]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::t>(i)[1][2],v_ofp_test.template get<P::t>(i)[1][2]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::t>(i)[2][0],v_ofp_test.template get<P::t>(i)[2][0]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::t>(i)[2][1],v_ofp_test.template get<P::t>(i)[2][1]);
		 BOOST_REQUIRE_EQUAL(dv_ofp_test.template get<P::t>(i)[2][2],v_ofp_test.template get<P::t>(i)[2][2]);
	 }

	 // try to push OTHER objects object

	 for (size_t i = 0 ; i < SECOND_PUSH ; i++)
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

		 v_ofp_test.add(p);

		 // Check the first element at every iteration

		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::x>(0),1.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::y>(0),2.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::z>(0),3.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::s>(0),4.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::v>(0)[0],1.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::v>(0)[1],2.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::v>(0)[2],7.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::t>(0)[0][0],10.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::t>(0)[0][1],13.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::t>(0)[0][2],8.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::t>(0)[1][0],19.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::t>(0)[1][1],23.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::t>(0)[1][2],5.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::t>(0)[2][0],4.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::t>(0)[2][1],3.0);
		 BOOST_REQUIRE_EQUAL(v_ofp_test.template get<P::t>(0)[2][2],11.0);
	 }

	 std::cout << "Vector unit test end" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

#endif
