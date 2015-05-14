/*
 * AdaptiveCellList_test.hpp
 *
 *  Created on: May 4, 2015
 *      Author: i-bird
 */

#ifndef ADAPTIVECELLLIST_TEST_HPP_
#define ADAPTIVECELLLIST_TEST_HPP_

#include "AdaptiveCellList.hpp"
#include "Grid/grid_sm.hpp"


template<unsigned int dim>
Point<dim+1,double> giveRadius(Point<dim,double> const & p)
{
  Point<dim+1, double> p_;
  for(unsigned int i=0; i<dim; ++i) p_.get(i) = p.get(i);
  p_[dim] = 0.5;
  return p_;
}

/*! \brief Test cell structure
 *
 * \tparam CellS
 *
 */
template<unsigned int dim, typename T, typename CellS> void Test_ar_cell_s()
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

		// particle one
		Point<dim,T> key = g_it.get();
		key = key * spacing + offset[0];

		cl1.add(giveRadius(key),id);
		++id;

		// particle two
		key = g_it.get();
		key = key * spacing + offset[1];

		cl1.add(giveRadius(key),id);
		++id;

		++g_it;
	}

	// check the cells are correctly filled

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
			// total

			total++;

			++NN;
		}

		BOOST_REQUIRE_EQUAL(total,openfpm::math::pow(3,dim));
		++g_it_s;
	}
}

BOOST_AUTO_TEST_SUITE( AdaptiveCellList_test )

BOOST_AUTO_TEST_CASE( AdaptiveCellList_use)
{
	std::cout << "Now testing the AR-List..." << std::endl;
	Test_ar_cell_s<3, double, AdaptiveCellList<3, double>>();
	std::cout << "Done testing the AR-List." << std::endl;
}

BOOST_AUTO_TEST_CASE( AdaptiveCellList_another_test)
{
	// Test case for your class

	//

}

BOOST_AUTO_TEST_SUITE_END()


#endif /* ADAPTIVECELLLIST_TEST_HPP_ */
