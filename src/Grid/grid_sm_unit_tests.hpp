/*
 * grid_sm_test.hpp
 *
 *  Created on: Dec 14, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_SM_UNIT_TESTS_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_SM_UNIT_TESTS_HPP_

#include "iterators/grid_key_dx_iterator_sub_bc.hpp"
#include "grid_key_dx_iterator_hilbert.hpp"

BOOST_AUTO_TEST_SUITE( grid_sm_test )


BOOST_AUTO_TEST_CASE( grid_sm_linearization )
{
	const grid_key_dx<3> key1(1,2,3);
	const grid_key_dx<3> zero(0,0,0);
	const grid_key_dx<3> seven(7,7,7);
	const comb<3> c({1,0,-1});
	size_t sz[3] = {8,8,8};

	grid_sm<3,int> gs(sz);
	size_t bc[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	long int lin = gs.LinId<CheckExistence>(key1,c.getComb(),bc);
	BOOST_REQUIRE_EQUAL(lin,146);
	lin = gs.LinId<CheckExistence>(zero,c.getComb(),bc);
	BOOST_REQUIRE_EQUAL(lin,-1);
	lin = gs.LinId<CheckExistence>(seven,c.getComb(),bc);
	BOOST_REQUIRE_EQUAL(lin,-1);

	for (size_t i = 0 ; i < 3 ; i++)
		bc[i] = PERIODIC;

	lin = gs.LinId<CheckExistence>(key1,c.getComb(),bc);
	BOOST_REQUIRE_EQUAL(lin,146);
	lin = gs.LinId<CheckExistence>(zero,c.getComb(),bc);
	BOOST_REQUIRE_EQUAL(lin,71);
	lin = gs.LinId<CheckExistence>(seven,c.getComb(),bc);
	BOOST_REQUIRE_EQUAL(lin,62);
}


BOOST_AUTO_TEST_CASE( grid_iterator_sub_p )
{
	const grid_key_dx<3> key1(4,4,4);
	const grid_key_dx<3> key2(-1,-1,-1);
	const grid_key_dx<3> key3(9,9,9);
	size_t sz[3] = {8,8,8};

	grid_sm<3,int> gs(sz);

	grid_key_dx_iterator_sub_bc<3> it(gs,key2,key1,{PERIODIC,PERIODIC,PERIODIC});

	size_t cnt = 0;

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0 ; i < 3 ; i++)
		{
			BOOST_REQUIRE_EQUAL(key.get(i) >= (long int)0,true);
			BOOST_REQUIRE_EQUAL(key.get(i) < (long int)sz[i],true);
		}

		cnt++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(cnt,216ul);

	grid_key_dx_iterator_sub_bc<3> it2(gs,key2,key3,{PERIODIC,PERIODIC,PERIODIC});

	cnt = 0;

	while (it2.isNext())
	{
		auto key = it2.get();

		for (size_t i = 0 ; i < 3 ; i++)
		{
			BOOST_REQUIRE_EQUAL(key.get(i) >= (long int)0,true);
			BOOST_REQUIRE_EQUAL(key.get(i) < (long int)sz[i],true);
		}

		cnt++;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(cnt,1331ul);

	cnt = 0;

	const grid_key_dx<3> key4(0,-1,0);
	const grid_key_dx<3> key5(2,2,2);

	grid_key_dx_iterator_sub_bc<3> it3(gs,key4,key5,{NON_PERIODIC,PERIODIC,NON_PERIODIC});

	while (it3.isNext())
	{
		auto key = it3.get();

		for (size_t i = 0 ; i < 3 ; i++)
		{
			BOOST_REQUIRE_EQUAL(key.get(i) >= (long int)0,true);
			BOOST_REQUIRE_EQUAL(key.get(i) < (long int)sz[i],true);
		}

		cnt++;

		++it3;
	}

	BOOST_REQUIRE_EQUAL(cnt,36ul);

	// bc non periodic with out-of-bound

	grid_key_dx_iterator_sub_bc<3,do_not_print_warning_on_adjustment<3>> it4(gs,key4,key5,{NON_PERIODIC,NON_PERIODIC,NON_PERIODIC});

	cnt = 0;

	while (it4.isNext())
	{
		auto key = it4.get();

		for (size_t i = 0 ; i < 3 ; i++)
		{
			BOOST_REQUIRE_EQUAL(key.get(i) >= (long int)0,true);
			BOOST_REQUIRE_EQUAL(key.get(i) < (long int)sz[i],true);
		}

		cnt++;

		++it4;
	}

	BOOST_REQUIRE_EQUAL(cnt,27ul);

	// case with no key

	const grid_key_dx<3> key6(-1,-1,-1);
	const grid_key_dx<3> key7(-1,-1,8);

	grid_key_dx_iterator_sub_bc<3,do_not_print_warning_on_adjustment<3>> it5(gs,key6,key7,{NON_PERIODIC,NON_PERIODIC,NON_PERIODIC});

	cnt = 0;

	while (it5.isNext())
	{
		auto key = it5.get();

		for (size_t i = 0 ; i < 3 ; i++)
		{
			BOOST_REQUIRE_EQUAL(key.get(i) >= 0,true);
			BOOST_REQUIRE_EQUAL(key.get(i) < (long int)sz[i],true);
		}

		cnt++;

		++it5;
	}

	BOOST_REQUIRE_EQUAL(cnt,0ul);
}

BOOST_AUTO_TEST_CASE( grid_key_dx_iterator_hilbert_test )
{
	// 2D test
	{
		size_t count = 0;

		//An order of a hilberts curve
		int32_t m = 2;

		grid_key_dx<2> start (0,0);

		//Create an iterator
		grid_key_dx_iterator_hilbert<2> h_it(m);

		while (h_it.isNext())
		{
			count++;

			++h_it;
		}

		//(2^m)^dim
		BOOST_REQUIRE_EQUAL(count, (size_t)16);

		h_it.reset();

		bool val = h_it.get() == start;

		BOOST_REQUIRE_EQUAL(val,true);
	}

	// 3D test
	{
		size_t count = 0;

		//An order of a hilberts curve
		int32_t m = 2;

		grid_key_dx<3> start (0,0,0);

		//Create an iterator
		grid_key_dx_iterator_hilbert<3> h_it(m);

		while (h_it.isNext())
		{
			count++;

			++h_it;
		}

		//(2^m)^dim
		BOOST_REQUIRE_EQUAL(count, (size_t)64);

		h_it.reset();

		bool val = h_it.get() == start;

		BOOST_REQUIRE_EQUAL(val,true);
	}
}


BOOST_AUTO_TEST_CASE( grid_iterator_sp_test )
{
	size_t sz[3] = {16,16,16};

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	grid_key_dx<3> start(2,2,2);
	grid_key_dx<3> stop(10,10,10);

	auto info = c3.getGrid();

	grid_key_dx_iterator_sp<3> it(info,info.LinId(start),info.LinId(stop));

	size_t count = 0;

	while (it.isNext())
	{
		count++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(count,2185ul);
}

BOOST_AUTO_TEST_CASE( grid_iterator_test_use)
{
	{
	//! [Grid iterator test usage]
	size_t count = 0;

	// Subdivisions
	size_t div[3] = {16,16,16};

	// grid info
	grid_sm<3,void> g_info(div);

	// Create a grid iterator
	grid_key_dx_iterator<3> g_it(g_info);

	// Iterate on all the elements
	while (g_it.isNext())
	{
		grid_key_dx<3> key = g_it.get();

		// set the grid key to zero without any reason ( to avoid warning compilations )
		key.zero();

		count++;

		++g_it;
	}

	BOOST_REQUIRE_EQUAL(count, (size_t)16*16*16);
	//! [Grid iterator test usage]
	}

	{
	size_t count = 0;
	// Iterate only on the internal elements

	//! [Sub-grid iterator test usage]
	// Subdivisions
	size_t div[3] = {16,16,16};

	// grid info
	grid_sm<3,void> g_info(div);

	grid_key_dx<3> start(1,1,1);
	grid_key_dx<3> stop(14,14,14);

	// Create a grid iterator (start and stop included)
	grid_key_dx_iterator_sub<3> g_it(g_info,start,stop);

	// Iterate on all the elements
	while (g_it.isNext())
	{
		grid_key_dx<3> key = g_it.get();

		// set the grid key to zero without any reason ( to avoid warning compilations )
		key.zero();

		count++;

		++g_it;
	}

	BOOST_REQUIRE_EQUAL(count, (size_t)14*14*14);

	//! [Sub-grid iterator test usage]

	// reset the iterator and check that it start from gk_start
	g_it.reset();

	bool val = g_it.get() == start;

	BOOST_REQUIRE_EQUAL(val,true);
	}
}

BOOST_AUTO_TEST_CASE( grid_sub_iterator_test )
{
	//! [Sub-grid iterator test usage]
	// Subdivisions
	size_t count = 0;
	typedef Point_test<float> p;

	size_t div[3] = {16,16,16};

	// grid info
	grid_cpu<3,Point_test<float>> g(div);
	g.setMemory();

	grid_key_dx<3> start(1,1,1);
	grid_key_dx<3> stop(14,14,14);

	// Create a grid iterator (start and stop included)
	auto g_it =  g.getIterator(start,stop);

	// Iterate on all the elements
	while (g_it.isNext())
	{
		grid_key_dx<3> key = g_it.get();

		// set the x value
		g.template get<p::x>(key) = 1.0;

		count++;

		++g_it;
	}

	BOOST_REQUIRE_EQUAL(count, (size_t)14*14*14);

	//! [Sub-grid iterator test usage]
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_DATA_SRC_GRID_GRID_SM_UNIT_TESTS_HPP_ */
