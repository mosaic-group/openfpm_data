/*
 * CellList_test.hpp
 *
 *  Created on: Mar 23, 2015
 *      Author: Pietro Incardona
 */

#include "CellList.hpp"
#include "Grid/grid.hpp"

#ifndef CELLLIST_TEST_HPP_
#define CELLLIST_TEST_HPP_

BOOST_AUTO_TEST_SUITE( CellList_test )

BOOST_AUTO_TEST_CASE( CellList_use)
{
	//Space where is living the Cell list
	SpaceBox<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Subdivisions
	size_t div[3] = {16,16,16};

	// grid info
	grid<3,void> g_info(div);

	// Origin
	Point<3,double> org({0.0,0.0,0.0});

	// id Cell list
	CellList<3,double> cl1(box,div,org);

	// Create a grid iterator
	grid_key_dx_iterator<3> g_it(g_info);

	// Iterate through each element
	// Add 1 element for each cell

	Point<3,double> end = box.getP2();
	Point<3,double> middle = end / div / 2.0;

	Point<3,double> offset[3] = {middle,middle,middle};

	// Create offset shift vectors
	offset[0].get(0) += 1.0 / div[0] / 8.0;
	offset[1].get(1) += 1.0 / div[1] / 8.0;
	offset[2].get(2) += 1.0 / div[2] / 8.0;

	size_t id = 0;

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<3,double> key = g_it.get();
		key += offset[0];

		cl1.addElement(key,id);

		key = g_it.get();
		key += offset[1];

		cl1.addElement(key,id);

		++g_it;
		++id;
	}

	// check the cell are correctly filled

	// Create a grid iterator
	g_it.reset();

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<3,double> key = g_it.get();
		key += offset[2];

		size_t cell = cl1.getCell(key);
		size_t n_ele = cl1.getNelements(cell);

		BOOST_REQUIRE_EQUAL(n_ele,2);
		BOOST_REQUIRE_EQUAL(cl1.getElement(cell,1) - cl1.getElement(cell,0),1);

		++g_it;
		++id;
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* CELLLIST_TEST_HPP_ */
