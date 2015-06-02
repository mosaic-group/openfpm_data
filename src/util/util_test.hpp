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

	typedef encapc<1,Point_test<float>,openfpm::vector<Point_test<float>>::memory_t> encap_src;
	typedef encapc<1,object<vboost_red>,openfpm::vector<object<vboost_red>>::memory_t> encap_dst;

	object<vboost> src;
	object<vboost_red> dst;

	// fill the source 0,1,4 with data

	boost::fusion::at_c<p::x>(src.data) = 1.0;
	boost::fusion::at_c<p::y>(src.data) = 2.0;

	for (size_t i = 0 ; i < 3 ;  i++)
		boost::fusion::at_c<p::v>(src.data)[i] = i + 5.0;

	object_copy<object<vboost>,object<vboost_red>,NORMAL,0,1,4>(src,dst);

	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<0>(dst.data),1.0);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<1>(dst.data),2.0);

	for (size_t i = 0 ; i < 3 ;  i++)
		BOOST_REQUIRE_EQUAL(boost::fusion::at_c<2>(dst.data)[i],i + 5.0);

	openfpm::vector<p> v_point;
	openfpm::vector<object<vboost_red>> v_point_red;

	v_point.resize(2);
	v_point_red.resize(2);

	v_point.template get<p::x>(0) = 1.0;
	v_point.template get<p::y>(0) = 2.0;

	for (size_t i = 0 ; i < 3 ;  i++)
		v_point.template get<p::v>(0)[i] = i + 5.0;

	auto src_e = v_point.get(0);
	auto dst_e = v_point_red.get(0);

	object_copy<encap_src,encap_dst,ENCAP,0,1,4>(src_e,dst_e);

	BOOST_REQUIRE_EQUAL(v_point_red.get(0).template get<0>(),1.0);
	BOOST_REQUIRE_EQUAL(v_point_red.get(0).template get<1>(),2.0);

	for (size_t i = 0 ; i < 3 ;  i++)
		BOOST_REQUIRE_EQUAL(v_point_red.get(0).template get<2>()[i],i + 5.0);

}

BOOST_AUTO_TEST_SUITE_END()

#endif /* UTIL_TEST_HPP_ */
