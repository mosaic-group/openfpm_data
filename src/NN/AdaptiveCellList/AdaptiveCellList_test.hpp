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
	/*
	//Space where is living the Cell list
	SpaceBox<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cl1(box, org);

	// Iterate through each element

	Point<dim,T> key({0.1,0.1,0.1});
	cl1.add(giveRadius(key),13);

	// particle two
	key = {0.9,0.9,0.9};
	cl1.add(giveRadius(key),42);

	// check the cells are correctly filled

	// Check we have 2 objects

	auto NN = cl1.template getNNIterator<NO_CHECK>(4711);
	size_t total = 0;

	while(NN.isNext())
	{
		// total

		total++;

		++NN;
	}

	BOOST_REQUIRE_EQUAL(total,2);
	
	//*/
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
