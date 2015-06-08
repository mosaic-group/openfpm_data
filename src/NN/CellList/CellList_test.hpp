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
	//! [Declare a cell list]
	//Space where is living the Cell list
	SpaceBox<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Subdivisions
	size_t div[dim] = {16,16,16};

	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cl2(box,div,org);
	//! [Declare a cell list]

	// grid info
	grid_sm<dim,void> g_info(div);

	// Test force reallocation in case of Cell list fast
	for (int i = 0 ; i < 32 ; i++)
		cl2.add(org,0);
	BOOST_REQUIRE_EQUAL(cl2.getNelements(cl2.getCell(org)),32);

	//! [Usage of cell list]

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
	for (size_t i = 0 ; i < dim ; i++)
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

	//! [Usage of cell list]

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

	//! [remove one particle from each cell]

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

	//! [remove one particle from each cell]

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


	//! [Usage of the neighborhood iterator]

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
			// total

			total++;

			++NN;
		}

		BOOST_REQUIRE_EQUAL(total,openfpm::math::pow(3,dim));
		++g_it_s;
	}

	//! [Usage of the neighborhood iterator]
}

BOOST_AUTO_TEST_SUITE( CellList_test )

BOOST_AUTO_TEST_CASE( CellDecomposer_use )
{
	//! [Cell decomposer use without shift]
	SpaceBox<3,double> box({0,0,0},{1,1,1});
	Point<3,double> p({0.5,0.5,0.5});

	// Number of cell on each dimension
	size_t div[] = {16,16,16};

	// Test Cell decomposer
	{
	CellDecomposer_sm<3,double> cd(box,div,0);
	size_t cell = cd.getCell(p);
	BOOST_REQUIRE_EQUAL(cell,8*16*16 + 8*16 + 8);
	auto key = cd.getCellGrid(p);
	BOOST_REQUIRE_EQUAL(key.get(0),8);
	BOOST_REQUIRE_EQUAL(key.get(1),8);
	BOOST_REQUIRE_EQUAL(key.get(2),8);
	}

	//! [Cell decomposer use without shift]

	//! [Test Cell decomposer with padding]
	{
	CellDecomposer_sm<3,double> cd(box,div,1);
	size_t cell = cd.getCell(p);
	BOOST_REQUIRE_EQUAL(cell,9*18*18 + 9*18 + 9);
	auto key = cd.getCellGrid(p);
	BOOST_REQUIRE_EQUAL(key.get(0),9);
	BOOST_REQUIRE_EQUAL(key.get(1),9);
	BOOST_REQUIRE_EQUAL(key.get(2),9);
	}

	//! [Test Cell decomposer with padding]

	//! [Test Cell decomposer with shift]
	{
	Point<3,double> sht({1.0,2.0,3.0});
	CellDecomposer_sm< 3,double,shift<3,double> > cd(box,div,sht,1);
	size_t cell = cd.getCell(p + sht);
	BOOST_REQUIRE_EQUAL(cell,9*18*18 + 9*18 + 9);
	auto key = cd.getCellGrid(p + sht);
	BOOST_REQUIRE_EQUAL(key.get(0),9);
	BOOST_REQUIRE_EQUAL(key.get(1),9);
	BOOST_REQUIRE_EQUAL(key.get(2),9);
	}

	//! [Test Cell decomposer with shift]
}

BOOST_AUTO_TEST_CASE( CellList_use)
{
	std::cout << "Test cell list" << "\n";

	Test_cell_s<3,double,CellList<3,double,FAST>>();
//	Test_cell_s<3,double,CellList<3,double,BALANCED>>();
//	Test_cell_s<3,double,CellList<3,double,MEMORY>>();

	std::cout << "End cell list" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* CELLLIST_TEST_HPP_ */
