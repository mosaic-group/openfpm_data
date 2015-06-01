/*
 * util_test.hpp
 *
 *  Created on: May 31, 2015
 *      Author: Pietro Incardona
 */

#ifndef UTIL_TEST_HPP_
#define UTIL_TEST_HPP_

#include "object_util.hpp"
#include "Point_test.hpp"

BOOST_AUTO_TEST_SUITE( util_test )

BOOST_AUTO_TEST_CASE( object_prop_copy )
{
	typedef Point_test<float> p;
	typedef Point_test<float>::type vboost;
	typedef object_creator<Point_test<float>::type,0,1,4>::type vboost_red;

	object<vboost> src;
	object<vboost_red> dst;

	// fill the source 0,1,4 with data

	boost::fusion::at_c<p::x>(src.data) = 1.0;
	boost::fusion::at_c<p::y>(src.data) = 2.0;

	for (size_t i = 0 ; i < 3 ;  i++)
		boost::fusion::at_c<p::v>(src.data)[i] = i + 5.0;

	object_copy<object<vboost>,object<vboost_red>,0,1,4>(src,dst);

	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<0>(dst.data),1.0);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<1>(dst.data),2.0);

	for (size_t i = 0 ; i < 3 ;  i++)
		BOOST_REQUIRE_EQUAL(boost::fusion::at_c<2>(dst.data)[i],i + 5.0);

}

BOOST_AUTO_TEST_SUITE_END()

#endif /* UTIL_TEST_HPP_ */
