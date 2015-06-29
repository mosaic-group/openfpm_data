/*
 * vector_test_util.hpp
 *
 *  Created on: Jun 19, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_TEST_UTIL_HPP_
#define VECTOR_TEST_UTIL_HPP_

#include "Point_test.hpp"
#include "Vector/map_vector.hpp"

//! [typedef point]
typedef Point_test<float> P;
//! [typedef point]

#ifdef TEST_COVERAGE_MODE
#define FIRST_PUSH 1000
#define SECOND_PUSH 1000
#else
#define FIRST_PUSH 1000000
#define SECOND_PUSH 1000000
#endif



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

openfpm::vector<Point_test<float>> allocate_openfpm(size_t n_ele)
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

	for (size_t i = 0 ; i < n_ele / 2 ; i++)
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

	for (size_t i = n_ele / 2 ; i < n_ele ; i++)
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




#endif /* VECTOR_TEST_UTIL_HPP_ */
