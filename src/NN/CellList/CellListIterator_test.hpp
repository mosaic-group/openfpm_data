/*
 * CellListIterator_test.hpp
 *
 *  Created on: May 7, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_TEST_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_TEST_HPP_

#include "NN/CellList/CellListIterator.hpp"
#include "NN/CellList/CellListFast_hilb.hpp"

BOOST_AUTO_TEST_SUITE( celllist_hilb_and_iterator_tests )

BOOST_AUTO_TEST_CASE( celllist_hilb_and_iterator_test )
{
	///////// INPUT DATA //////////

	const size_t dim = 3;

	size_t div[dim] = {4,5,6};

	//Number of particles
	size_t k = 300;

	///////////////////////////////

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Initialize a cell list
	CellList_hilb<dim,float> NN;

	NN.Initialize(box,div,k*0.9,1);

	float pos[dim];

	//Fill with particles
	for (size_t i = 0; i < k*0.9; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			pos[j] = rand()/double(RAND_MAX);
		}
		NN.add(pos,i);
	}

	//Test the iterator
	auto it_cl = NN.getIterator();

	size_t count = 0;

	while (it_cl.isNext())
	{
		auto p_key = it_cl.get();

		BOOST_REQUIRE(p_key < NN.get_gm());

		count++;
		++it_cl;
	}

	BOOST_REQUIRE_EQUAL(count,NN.get_gm());

	// Save cell keys

	NN.getKeys().save("NN_hilb_keys");

	// Load previous results and check equality

	openfpm::vector<size_t> keys_old;

	keys_old.load("NN_hilb_keys");

	for (size_t i = 0; i < keys_old.size(); i++)
	{
		size_t a1 = keys_old.get(i);
		size_t a2 = NN.getKeys().get(i);

		BOOST_REQUIRE_EQUAL(a1,a2);
	}

	size_t s1 = keys_old.size();
	size_t s2 = NN.getKeys().size();

	BOOST_REQUIRE_EQUAL(s1,s2);
}


BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_TEST_HPP_ */
