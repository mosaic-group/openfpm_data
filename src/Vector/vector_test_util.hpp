/*
 * vector_test_util.hpp
 *
 *  Created on: Jun 19, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_TEST_UTIL_HPP_
#define VECTOR_TEST_UTIL_HPP_

#include "config.h"
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

#include "timer.hpp"

std::vector<Point_orig<float>> allocate_stl()
{
	std::vector<Point_orig<float>> v_stl_test;

	// Now fill the vector

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

	return v_stl_test;
}

template <typename T>
openfpm::vector<T> allocate_openfpm_primitive(size_t n, size_t fill)
{
	openfpm::vector<T> v;

	for (size_t i = 0 ; i < n ; i++)
		v.add(fill);

	return v;
}

openfpm::vector<Point_test<float>> allocate_openfpm_fill(size_t n, size_t fill)
{
	// Test point
	typedef Point_test<float> p;

	Point_test<float> pt;
	openfpm::vector<Point_test<float>> v_send;

	pt.setx(fill);
	pt.sety(fill);
	pt.setz(fill);
	pt.sets(fill);

	pt.setv(0,fill);
	pt.setv(1,fill);
	pt.setv(2,fill);

	pt.sett(0,0,fill);
	pt.sett(0,1,fill);
	pt.sett(0,2,fill);
	pt.sett(1,0,fill);
	pt.sett(1,1,fill);
	pt.sett(1,2,fill);
	pt.sett(2,0,fill);
	pt.sett(2,1,fill);
	pt.sett(2,2,fill);


	// ADD n elements
	for (size_t i = 0 ; i < n ; i++)
		v_send.add(pt);

	return v_send;
}


openfpm::vector<Point_test<float>> allocate_openfpm(size_t n_ele)
{
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

	return v_ofp_test;
}

#endif /* VECTOR_TEST_UTIL_HPP_ */
