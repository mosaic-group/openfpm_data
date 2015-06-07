#ifndef VECTOR_UNIT_TESTS_HPP
#define VECTOR_UNIT_TESTS_HPP

#include "map_vector.hpp"
#include "Point_test.hpp"
#include "memory/PreAllocHeapMemory.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "memory/PtrMemory.hpp"
#include <cstring>
#include "Space/Shape/Point.hpp"
#include "util/object_creator.hpp"

#define FIRST_PUSH 1000000
#define SECOND_PUSH 1000000

BOOST_AUTO_TEST_SUITE( vector_test )

//! [typedef point]
typedef Point_test<float> P;
//! [typedef point]

std::vector<Point_orig<float>> allocate_stl()
{
	std::vector<Point_orig<float>> v_stl_test;

	// Now fill the vector

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
	#ifdef VERBOSE_TEST
	timespec ts_start;
	// clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
	clock_gettime(CLOCK_REALTIME, &ts_start); // Works on Linux
	#endif

	//! [Create add and access]
	openfpm::vector<Point_test<float>> v_ofp_test;

	//! [Point declaration]
	Point_test<float> p;
	//! [Point declaration]

	p.setx(1.0);
	p.sety(2.0);
	p.setz(3.0);
	p.sets(4.0);

	// push objects

	for (size_t i = 0 ; i < FIRST_PUSH / 2 ; i++)
	{
		// Modify the point

		//! [Point usage]
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
		//! [Point usage]

		// add p

		v_ofp_test.add(p);
	}

	for (size_t i = FIRST_PUSH / 2 ; i < FIRST_PUSH ; i++)
	{
		v_ofp_test.add();

		size_t last = v_ofp_test.size()-1;

		// Modify the point

		v_ofp_test.get<P::v>(last)[0] = 1.0 + i;
		v_ofp_test.get<P::v>(last)[1] = 2.0 + i;
		v_ofp_test.get<P::v>(last)[2] = 7.0 + i;

		v_ofp_test.get<P::t>(last)[0][0] = 10.0 + i;
		v_ofp_test.get<P::t>(last)[0][1] = 13.0 + i;
		v_ofp_test.get<P::t>(last)[0][2] = 8.0 + i;
		v_ofp_test.get<P::t>(last)[1][0] = 19.0 + i;
		v_ofp_test.get<P::t>(last)[1][1] = 23.0 + i;
		v_ofp_test.get<P::t>(last)[1][2] = 5.0 + i;
		v_ofp_test.get<P::t>(last)[2][0] = 4.0 + i;
		v_ofp_test.get<P::t>(last)[2][1] = 3.0 + i;
		v_ofp_test.get<P::t>(last)[2][2] = 11.0 + i;
	}

	//! [Create add and access]

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

BOOST_AUTO_TEST_CASE( vector_std_utility )
{
	//! [Create add and access stl]

	// Create a vector with 13 element
	openfpm::vector<size_t> pb(13);

	// add at the end some othe element
	pb.add(0);
	pb.add(1);
	pb.add(2);

	// access the vector
	for (size_t i = 0 ;  i < 16 ; i++)
	{
		pb.get(i) = i+1;
	}

	//! [Create add and access stl]

	pb.fill(0);

	// Check is zero
	for (size_t i = 0 ;  i < 16 ; i++)
	{
		BOOST_REQUIRE_EQUAL(pb.get(i),0);
	}

}

size_t alloc[] = {235,345,520};
size_t n_alloc = sizeof(alloc)/sizeof(size_t);

BOOST_AUTO_TEST_CASE ( vector_prealloc_ext )
{
	// Memory for the ghost sending buffer
	HeapMemory mem;

	// sequence of pre-allocation pattern
	std::vector<size_t> pap;

	size_t total = 0;

	// Calculate the total size required for the sending buffer
	for (size_t i = 0 ; i < n_alloc ; i++)
	{
		size_t alloc_ele = openfpm::vector<Point_test<float>>::calculateMem(alloc[i],0);
		pap.push_back(alloc_ele);
		total += alloc_ele;
	}

	// Create an object of preallocated memory
	ExtPreAlloc<HeapMemory> * prAlloc = new ExtPreAlloc<HeapMemory>(pap,mem);

	typedef openfpm::vector<Point_test<float>,openfpm::device_cpu<Point_test<float>>,ExtPreAlloc<HeapMemory>> send_vector;

	// create a vector of send_vector (ExtPreAlloc warrant that all the created vector are contiguous)
	openfpm::vector<send_vector> g_send;

	// resize
	g_send.resize(n_alloc);

	// Number of allocation
	for (size_t i = 0 ; i < n_alloc ; i++)
	{
		// set the preallocated memory to ensure contiguity
		g_send.get(i).setMemory(*prAlloc);

		// resize the sending vector (No allocation is produced)
		g_send.get(i).resize(alloc[i]);
	}

	// Fill the send buffer with one
	for (size_t i = 0 ; i < n_alloc ; i++)
	{
		auto it = g_send.get(i).getIterator();
		auto & v = g_send.get(i);

		while(it.isNext())
		{
			auto kk = it.get();

			v.template get<P::x>(kk) = 1.0f;
			v.template get<P::y>(kk) = 1.0f;
			v.template get<P::z>(kk) = 1.0f;
			v.template get<P::s>(kk) = 1.0f;

			v.template get<P::v>(kk)[0] = 1.0f;
			v.template get<P::v>(kk)[1] = 1.0f;
			v.template get<P::v>(kk)[2] = 1.0f;

			v.template get<P::t>(kk)[0][0] = 1.0f;
			v.template get<P::t>(kk)[0][1] = 1.0f;
			v.template get<P::t>(kk)[0][2] = 1.0f;
			v.template get<P::t>(kk)[1][0] = 1.0f;
			v.template get<P::t>(kk)[1][1] = 1.0f;
			v.template get<P::t>(kk)[1][2] = 1.0f;
			v.template get<P::t>(kk)[2][0] = 1.0f;
			v.template get<P::t>(kk)[2][1] = 1.0f;
			v.template get<P::t>(kk)[2][2] = 1.0f;

			++it;
		}
	}

	// check that HeapMemory contain ones in the right position
	float * ptr = (float *) mem.getPointer();
	size_t offset = 0;

	for (size_t i = 0 ; i < n_alloc ; i++)
	{
		for (size_t j = 0 ; j < alloc[i] ; j++)
			BOOST_REQUIRE_EQUAL(ptr[j + offset/sizeof(float)],1.0f);

		offset += pap[i];
	}
}

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

		// set the size and allocate, using mem warrant that pos and prp is contiguous
		pb.get(i).pos.resize(1024);
		pb.get(i).prp.resize(1024);
	}
}


BOOST_AUTO_TEST_CASE( object_test_creator )
{
	bool tst = std::is_same< typename object_creator<Point_test<float>::type,0,1,5>::type, typename boost::fusion::vector3<float,float,float[3][3]> >::value;

	BOOST_REQUIRE_EQUAL(tst , true);
}

#define V_REM_PUSH 1024

BOOST_AUTO_TEST_CASE(vector_remove )
{
	typedef Point_test<float> p;

	//! [Create push and multiple remove]

	openfpm::vector<Point_test<float>> v1;

	for (size_t i = 0 ; i < V_REM_PUSH ; i++)
	{
		// Point
		Point_test<float> p;
		p.setx(i);

		v1.add(p);
	}

	{
	openfpm::vector<size_t> rem;
	rem.add(0);
	rem.add(1);
	rem.add(2);
	rem.add(3);

	v1.remove(rem);
	}

	//! [Create push and multiple remove]

	BOOST_REQUIRE_EQUAL(v1.size(),1020);
	BOOST_REQUIRE_EQUAL(v1.template get<p::x>(0),4);

	{
	openfpm::vector<size_t> rem;
	rem.add(v1.size()-3);
	rem.add(v1.size()-2);
	rem.add(v1.size()-1);
	rem.add(v1.size());

	v1.remove(rem);
	}

	BOOST_REQUIRE_EQUAL(v1.size(),1016);
	BOOST_REQUIRE_EQUAL(v1.template get<p::x>(v1.size()-1),1019);

	{
	openfpm::vector<size_t> rem;
	for (size_t i = 0 ; i < (V_REM_PUSH - 8) / 2 ; i++)
		rem.add(i * 2);

	// remove all the even number
	v1.remove(rem);
	}

	BOOST_REQUIRE_EQUAL(v1.size(),508);

	// Check only odd
	for (size_t i = 0 ; i < v1.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL((size_t)v1.template get<p::x>(v1.size()-1) % 2, 1);
	}
}

BOOST_AUTO_TEST_CASE( vector_memory_repr )
{
	// create a vector
	openfpm::vector<Point_test<float>> v1;

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

		v1.add(p);
	}

	PtrMemory * ptr1 = new PtrMemory(v1.getPointer(),sizeof(Point_test<float>)*FIRST_PUSH);

	// create vector representation to a piece of memory already allocated

	openfpm::vector<Point_test<float>,openfpm::device_cpu<Point_test<float>>,PtrMemory,openfpm::grow_policy_identity> v2;

	v2.setMemory(*ptr1);

	v2.resize(FIRST_PUSH);

	// check

	// Check if the duplicated vector match

	for (size_t i = 0 ; i < FIRST_PUSH ; i++)
	{
		BOOST_REQUIRE_EQUAL(v1.template get<P::v>(i)[0],v2.template get<P::v>(i)[0]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::v>(i)[1],v2.template get<P::v>(i)[1]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::v>(i)[2],v2.template get<P::v>(i)[2]);

		BOOST_REQUIRE_EQUAL(v1.template get<P::t>(i)[0][0],v2.template get<P::t>(i)[0][0]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::t>(i)[0][1],v2.template get<P::t>(i)[0][1]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::t>(i)[0][2],v2.template get<P::t>(i)[0][2]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::t>(i)[1][0],v2.template get<P::t>(i)[1][0]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::t>(i)[1][1],v2.template get<P::t>(i)[1][1]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::t>(i)[1][2],v2.template get<P::t>(i)[1][2]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::t>(i)[2][0],v2.template get<P::t>(i)[2][0]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::t>(i)[2][1],v2.template get<P::t>(i)[2][1]);
		BOOST_REQUIRE_EQUAL(v1.template get<P::t>(i)[2][2],v2.template get<P::t>(i)[2][2]);
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif
