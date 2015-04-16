#ifndef VECTOR_UNIT_TESTS_HPP
#define VECTOR_UNIT_TESTS_HPP

#include "map_vector.hpp"
#include "Point_test.hpp"
#include "memory/PreAllocHeapMemory.hpp"
#include <cstring>
#include "Space/Shape/Point.hpp"

#define FIRST_PUSH 1000000
#define SECOND_PUSH 1000000

BOOST_AUTO_TEST_SUITE( vector_test )

typedef Point_test<float> P;

std::vector<Point_orig<float>> allocate_stl()
{
	std::vector<Point_orig<float>> v_stl_test;

	// Now fill the STL vector

	#ifdef VERBOSE_TEST
	timespec ts_start;
	// clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
	clock_gettime(CLOCK_REALTIME, &ts_start); // Works on Linux
	#endif

	Point_orig<float> po;
	po.setx(1.0);
	po.sety(2.0);
	po.setz(3.0);
	po.sets(4.0);

	for (size_t i = 0 ; i < FIRST_PUSH ; i++)
	{
		// Modify po

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

	#ifdef VERBOSE_TEST
	timespec end_time;
	clock_gettime(CLOCK_REALTIME, &end_time); // Works on Linux
	float time_dif =(float)( end_time.tv_sec - ts_start.tv_sec  + (double)(end_time.tv_nsec - ts_start.tv_nsec)/1000000000.0 );

	std::cout << "STL : " << FIRST_PUSH << " add " << "  Time: " << time_dif << " s  " << "\n";
	#endif

	return v_stl_test;
}

openfpm::vector<Point_test<float>> allocate_openfpm()
{
	openfpm::vector<Point_test<float>> v_ofp_test;

	#ifdef VERBOSE_TEST
	timespec ts_start;
	// clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
	clock_gettime(CLOCK_REALTIME, &ts_start); // Works on Linux
	#endif

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

	#ifdef VERBOSE_TEST
	timespec end_time;
	clock_gettime(CLOCK_REALTIME, &end_time); // Works on Linux
	float time_dif =(float)( end_time.tv_sec - ts_start.tv_sec  + (double)(end_time.tv_nsec - ts_start.tv_nsec)/1000000000.0 );

	std::cout << "OPENFPM : " << FIRST_PUSH << " add " << "  Time: " << time_dif << " s  " << "\n";
	#endif

	return v_ofp_test;
}

// Test the openfpm vector

BOOST_AUTO_TEST_CASE( vector_use)
{
	std::cout << "Vector unit test start" << "\n";

	std::vector<Point_orig<float>> v_stl_test = allocate_stl();
	openfpm::vector<Point_test<float>> v_ofp_test = allocate_openfpm();

	// try to duplicate the vector

	openfpm::vector<Point_test<float>> dv_ofp_test = v_ofp_test.duplicate();

	// Check if the STL and openfpm match

	for (size_t i = 0; i < FIRST_PUSH; i++)
	{
		BOOST_REQUIRE_EQUAL(v_stl_test[i].v[0],v_ofp_test.template get<P::v>(i)[0]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].v[1],v_ofp_test.template get<P::v>(i)[1]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].v[2],v_ofp_test.template get<P::v>(i)[2]);

		BOOST_REQUIRE_EQUAL(v_stl_test[i].t[0][0],v_ofp_test.template get<P::t>(i)[0][0]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].t[0][1],v_ofp_test.template get<P::t>(i)[0][1]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].t[0][2],v_ofp_test.template get<P::t>(i)[0][2]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].t[1][0],v_ofp_test.template get<P::t>(i)[1][0]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].t[1][1],v_ofp_test.template get<P::t>(i)[1][1]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].t[1][2],v_ofp_test.template get<P::t>(i)[1][2]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].t[2][0],v_ofp_test.template get<P::t>(i)[2][0]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].t[2][1],v_ofp_test.template get<P::t>(i)[2][1]);
		BOOST_REQUIRE_EQUAL(v_stl_test[i].t[2][2],v_ofp_test.template get<P::t>(i)[2][2]);
	}

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

	std::cout << "Vector unit test end" << "\n";
}

// Pre alloc test

struct pre_test
{
	//! position vector
	openfpm::vector<Point<2,float>,openfpm::device_cpu<Point<2,float>>,PreAllocHeapMemory<2>,openfpm::grow_policy_identity> pos;
	//! properties vector
	openfpm::vector<Point_test<float>,openfpm::device_cpu<Point_test<float>>,PreAllocHeapMemory<2>,openfpm::grow_policy_identity> prp;
};


BOOST_AUTO_TEST_CASE( vector_prealloc )
{
	openfpm::vector<pre_test> pb(3);

	for (size_t i = 0 ;  i < 3 ; i++)
	{
		// Create the size required to store the particles position and properties to communicate
		size_t s1 = openfpm::vector<Point<2,float>>::calculateMem(1024,0);
		size_t s2 = openfpm::vector<Point_test<float>>::calculateMem(1024,0);

		// Preallocate the memory
		size_t sz[2] = {s1,s2};
		PreAllocHeapMemory<2> * mem = new PreAllocHeapMemory<2>(sz);

		// Set the memory allocator
		pb.get(i).pos.setMemory(*mem);
		pb.get(i).prp.setMemory(*mem);

		// set the size and allocate, using mem warant that pos and prp is contiguous
		pb.get(i).pos.resize(1024);
		pb.get(i).prp.resize(1024);
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif
