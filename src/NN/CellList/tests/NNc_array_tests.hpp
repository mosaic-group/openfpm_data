/*
 * NNc_array_tests.hpp
 *
 *  Created on: Feb 6, 2018
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_TESTS_NNC_ARRAY_TESTS_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_TESTS_NNC_ARRAY_TESTS_HPP_

#include "NN/CellList/NNc_array.hpp"

BOOST_AUTO_TEST_SUITE( NNc_array_tests )

BOOST_AUTO_TEST_CASE( NNc_array_tests_use )
{
	{
	NNc_array<2,9,false> nnc;
	NNc_array<2,9,true> nnc2;

	size_t sz[2] = {13,17};

	nnc.set_size(sz);
	nnc2.set_size(sz);
	nnc.init_full();
	nnc2.init_full();

	for (size_t i = 0 ; i < 9 ; i++)
	{
		BOOST_REQUIRE_EQUAL(nnc[i],nnc2[i]);
	}
	}

	{
	NNc_array<3,27,false> nnc;
	NNc_array<3,27,true> nnc2;

	size_t sz[3] = {13,17,11};

	nnc.set_size(sz);
	nnc2.set_size(sz);
	nnc.init_full();
	nnc2.init_full();

	for (size_t i = 0 ; i < 27 ; i++)
	{
		BOOST_REQUIRE_EQUAL(nnc[i],nnc2[i]);
	}
	}

	{
	NNc_array<5,243,false> nnc;
	NNc_array<5,243,true> nnc2;

	size_t sz[5] = {13,17,11,7,7};

	nnc.set_size(sz);
	nnc2.set_size(sz);
	nnc.init_full();
	nnc2.init_full();

	for (size_t i = 0 ; i < 243 ; i++)
	{
		BOOST_REQUIRE_EQUAL(nnc[i],nnc2[i]);
	}
	}


	{
	NNc_array<2,9,false> nnc;
	NNc_array<2,9,true> nnc2;

	size_t sz[2] = {13,17};

	nnc.set_size(sz);
	nnc2.set_size(sz);
	nnc.init_sym();
	nnc2.init_sym();

	for (size_t i = 0 ; i < 5 ; i++)
	{
		BOOST_REQUIRE_EQUAL(nnc[i],nnc2[i]);
	}
	}

	{
	NNc_array<3,27,false> nnc;
	NNc_array<3,27,true> nnc2;

	size_t sz[3] = {13,17,11};

	nnc.set_size(sz);
	nnc2.set_size(sz);
	nnc.init_full();
	nnc2.init_full();

	for (size_t i = 0 ; i < 27 ; i++)
	{
		BOOST_REQUIRE_EQUAL(nnc[i],nnc2[i]);
	}
	}

	{
	NNc_array<5,243,false> nnc;
	NNc_array<5,243,true> nnc2;

	size_t sz[5] = {13,17,11,7,7};

	nnc.set_size(sz);
	nnc2.set_size(sz);
	nnc.init_full();
	nnc2.init_full();

	for (size_t i = 0 ; i < 243 ; i++)
	{
		BOOST_REQUIRE_EQUAL(nnc[i],nnc2[i]);
	}
	}
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_TESTS_NNC_ARRAY_TESTS_HPP_ */
