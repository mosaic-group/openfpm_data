/*
 * CellList_test.hpp
 *
 *  Created on: Mar 23, 2015
 *      Author: Pietro Incardona
 */

#include "CellList.hpp"
#include "Grid/grid_sm.hpp"

#ifndef CELLLIST_TEST_HPP_
#define CELLLIST_TEST_HPP_


/*! \brief Test cell structure
 *
 * \tparam CellS
 *
 */
template<unsigned int dim, typename T, typename CellS> void Test_cell_s()
{
	//Space where is living the Cell list
	SpaceBox<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Subdivisions
	size_t div[dim] = {16,16,16};

	// grid info
	grid_sm<dim,void> g_info(div);

	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cl1(box,div,org);

	// Create a grid iterator
	grid_key_dx_iterator<dim> g_it(g_info);

	// Iterate through each element
	// Add 1 element for each cell

	// Usefull definition of points
	Point<dim,T> end = box.getP2();
	Point<dim,T> middle = end / div / 2.0;
	Point<dim,T> spacing = end / div;

	Point<dim,T> offset[dim] = {middle,middle,middle};

	// Create offset shift vectors
	for (int i = 0 ; i < dim ; i++)
	{
		offset[i].get(i) += (1.0 / div[i]) / 8.0;
	}

	size_t id = 0;

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<dim,T> key = g_it.get();
		key = key * spacing + offset[0];

		cl1.add(key,id);
		++id;

		key = g_it.get();
		key = key * spacing + offset[1];

		cl1.add(key,id);
		++id;

		++g_it;
	}

	// check the cell are correctly filled

	// reset iterator
	g_it.reset();

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<dim,T> key = g_it.get();
		key = key * spacing + offset[2];

		size_t cell = cl1.getCell(key);
		size_t n_ele = cl1.getNelements(cell);

		BOOST_REQUIRE_EQUAL(n_ele,2);
		BOOST_REQUIRE_EQUAL(cl1.get(cell,1) - cl1.get(cell,0),1);

		++g_it;
	}

	// reset itarator
	g_it.reset();

	// remove one particle from each cell

	while (g_it.isNext())
	{
		// remove 1 particle on each cell

		Point<dim,T> key = g_it.get();
		key = key * spacing + offset[0];

		auto cell = cl1.getCell(key);

		// Remove the first particle in the cell
		cl1.remove(cell,0);
		++g_it;
	}

	// Check we have 1 object per cell
	g_it.reset();

	while (g_it.isNext())
	{
		// remove 1 particle on each cell

		Point<dim,T> key = g_it.get();
		key = key * spacing + offset[0];

		auto cell = cl1.getCell(key);
		size_t n_ele = cl1.getNelements(cell);

		BOOST_REQUIRE_EQUAL(n_ele,1);
		++g_it;
	}


	// Check the neighborhood iterator on the internal grid (They do not wotk on external grid)

	// Check we have 1 object per cell

	// Create a grid iterator
	grid_key_dx<dim> p1(1,1,1);
	grid_key_dx<dim> p2(div[0]-2,div[1]-2,div[2]-2);
	grid_key_dx_iterator_sub<dim> g_it_s(g_info,p1,p2);

	while (g_it_s.isNext())
	{
		// remove 1 particle on each cell

		Point<dim,T> key = g_it_s.get();
		key = key * spacing + offset[0];

		auto NN = cl1.template getNNIterator<NO_CHECK>(cl1.getCell(key));
		size_t total = 0;

		while(NN.isNext())
		{
			size_t id = NN.get();

			// total

			total++;

			++NN;
		}

		BOOST_REQUIRE_EQUAL(total,openfpm::math::pow(3,dim));
		++g_it_s;
	}
}

BOOST_AUTO_TEST_SUITE( CellList_test )

BOOST_AUTO_TEST_CASE( CellList_use)
{
	std::cout << "Test cell list" << "\n";

	Test_cell_s<3,double,CellList<3,double>>();

	std::cout << "End cell list" << "\n";

	//Space where is living the Cell list
/*	SpaceBox<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Subdivisions
	size_t div[3] = {16,16,16};

	// grid info
	grid_sm<3,void> g_info(div);

	// Origin
	Point<3,double> org({0.0,0.0,0.0});

	// id Cell list
	CellList<3,double> cl1(box,div,org);

	// Create a grid iterator
	grid_key_dx_iterator<3> g_it(g_info);

	// Iterate through each element
	// Add 1 element for each cell

	// Usefull definition of points
	Point<3,double> end = box.getP2();
	Point<3,double> middle = end / div / 2.0;
	Point<3,double> spacing = end / div;

	Point<3,double> offset[3] = {middle,middle,middle};

	// Create offset shift vectors
	offset[0].get(0) += (1.0 / div[0]) / 8.0;
	offset[1].get(1) += (1.0 / div[1]) / 8.0;
	offset[2].get(2) += (1.0 / div[2]) / 8.0;

	size_t id = 0;

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<3,double> key = g_it.get();
		key = key * spacing + offset[0];

		cl1.add(key,id);
		++id;

		key = g_it.get();
		key = key * spacing + offset[1];

		cl1.add(key,id);
		++id;

		++g_it;
	}

	// check the cell are correctly filled

	// reset iterator
	g_it.reset();

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<3,double> key = g_it.get();
		key = key * spacing + offset[2];

		size_t cell = cl1.getCell(key);
		size_t n_ele = cl1.getNelements(cell);

		BOOST_REQUIRE_EQUAL(n_ele,2);
		BOOST_REQUIRE_EQUAL(cl1.get(cell,1) - cl1.get(cell,0),1);

		++g_it;
	}

	// reset itarator
	g_it.reset();

	// remove one particle from each cell

	while (g_it.isNext())
	{
		// remove 1 particle on each cell

		Point<3,double> key = g_it.get();
		key = key * spacing + offset[0];

		auto cell = cl1.getCell(key);

		// Remove the first particle in the cell
		cl1.remove(cell,0);
		++g_it;
	}

	// Check we have 1 object per cell
	g_it.reset();

	while (g_it.isNext())
	{
		// remove 1 particle on each cell

		Point<3,double> key = g_it.get();
		key = key * spacing + offset[0];

		auto cell = cl1.getCell(key);
		size_t n_ele = cl1.getNelements(cell);

		BOOST_REQUIRE_EQUAL(n_ele,1);
		++g_it;
	}


	// Check the neighborhood iterator on the internal grid (They do not wotk on external grid)

	// Check we have 1 object per cell

	// Create a grid iterator
	grid_key_dx<3> p1(1,1,1);
	grid_key_dx<3> p2(div[0]-1,div[1]-1,div[2]-1);
	grid_key_dx_iterator_sub<3> g_it_s(g_info,p1,p2);

	while (g_it_s.isNext())
	{
		// remove 1 particle on each cell

		Point<3,double> key = g_it_s.get();
		key = key * spacing + offset[0];

		auto NN = cl1.template getNNIterator<NO_CHECK>(cl1.getCell(key));
		size_t total = 0;

		while(NN.isNext())
		{
			size_t id = NN.get();

			// total

			total++;

			++NN;
		}

		BOOST_REQUIRE_EQUAL(total,openfpm::math::pow(3,3));
		++g_it;
	}*/

	// Test the cell list
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* CELLLIST_TEST_HPP_ */
