/*
 * meta_cc_unit_tests.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_META_CC_UNIT_TESTS_HPP_
#define OPENFPM_DATA_SRC_UTIL_META_CC_UNIT_TESTS_HPP_

#include "meta_copy.hpp"
#include "meta_compare.hpp"
#include "data_type/aggregate.hpp"
#include "Point_test.hpp"

BOOST_AUTO_TEST_SUITE( util_test )

BOOST_AUTO_TEST_CASE( meta_copy_compare_test )
{
	{
	//! [Usage of meta copy and compare for primitives]

	float f_src = 1.0;
	float f_dst;

	meta_copy<float>::meta_copy_(f_src,f_dst);

	BOOST_REQUIRE_EQUAL(f_src,f_dst);

	bool ret = meta_compare<float>::meta_compare_f(f_src,f_dst);

	BOOST_REQUIRE_EQUAL(ret,true);

	//! [Usage of meta copy and compare for primitives]

	f_dst = 2.0;
	ret = meta_compare<float>::meta_compare_f(f_src,f_dst);
	BOOST_REQUIRE_EQUAL(ret,false);

	}

	{
	//! [Usage of meta copy and compare for array of primitives]

	float f_src[2][3] = {{1.0,2.9,4.0},{2.3,4.4,9.0}};
	float f_dst[2][3];

	meta_copy<float[2][3]>::meta_copy_(f_src,f_dst);

	bool ret = meta_compare<float[2][3]>::meta_compare_f(f_src,f_dst);

	BOOST_REQUIRE_EQUAL(ret,true);

	//! [Usage of meta copy and compare for array of primitives]

	f_dst[1][2] = 5.0;
	ret = meta_compare<float[2][3]>::meta_compare_f(f_src,f_dst);
	BOOST_REQUIRE_EQUAL(ret,false);

	}

	{
	//! [Usage of meta copy and compare for openfpm aggregates]

	aggregate<float,int,float[3]> agg1;
	aggregate<float,int,float[3]> agg2;

	boost::fusion::at_c<0>(agg1.data) = 1.0;
	boost::fusion::at_c<1>(agg1.data) = 2.0;
	boost::fusion::at_c<2>(agg1.data)[0] = 3.0;
	boost::fusion::at_c<2>(agg1.data)[1] = 4.0;
	boost::fusion::at_c<2>(agg1.data)[2] = 5.0;

	meta_copy<aggregate<float,int,float[3]>>::meta_copy_(agg1,agg2);

	bool ret = meta_compare<aggregate<float,int,float[3]>>::meta_compare_f(agg1,agg2);

	BOOST_REQUIRE_EQUAL(ret,true);

	//! [Usage of meta copy and compare for openfpm aggregates]

	boost::fusion::at_c<2>(agg2.data)[2] = 2.0;
	ret = meta_compare<aggregate<float,int,float[3]>>::meta_compare_f(agg1,agg2);

	BOOST_REQUIRE_EQUAL(ret,false);

	}

	{
	//! [Usage of meta copy and compare for complex object]

	std::string s_src("Test string");
	std::string s_dst;

	meta_copy<std::string>::meta_copy_(s_src,s_dst);

	BOOST_REQUIRE_EQUAL(s_src,s_dst);

	bool ret = meta_compare<std::string>::meta_compare_f(s_src,s_dst);

	BOOST_REQUIRE_EQUAL(ret,true);

	//! [Usage of meta copy and compare for complex object]

	s_dst = std::string("Test string2");
	ret = meta_compare<std::string>::meta_compare_f(s_src,s_dst);
	BOOST_REQUIRE_EQUAL(ret,false);

	}

	{
	//! [Usage of meta copy and compare for complex aggregates object]

	aggregate<std::string,std::vector<float>,std::map<size_t,std::string>,std::string[3]> a_src;

	// fill the complex aggregates
	a_src.template get<0>() = std::string("Test string");

	a_src.template get<1>().push_back(5.0);
	a_src.template get<1>().push_back(15.0);
	a_src.template get<1>().push_back(45.0);
	a_src.template get<1>().push_back(7.0);

	a_src.template get<2>()[0] = std::string("Test string 2");
	a_src.template get<2>()[10] = std::string("Test string 3");
	a_src.template get<2>()[9] = std::string("Test string 4");
	a_src.template get<2>()[1] = std::string("Test string 5");

	a_src.template get<3>()[0] = std::string("Last string 9");
	a_src.template get<3>()[1] = std::string("Last string 10");
	a_src.template get<3>()[2] = std::string("Last string 11");

	aggregate<std::string,std::vector<float>,std::map<size_t,std::string>,std::string[3]> a_dst;

	meta_copy<aggregate<std::string,std::vector<float>,std::map<size_t,std::string>,std::string[3]>>::meta_copy_(a_src,a_dst);

	bool ret = meta_compare<aggregate<std::string,std::vector<float>,std::map<size_t,std::string>,std::string[3]>>::meta_compare_f(a_src,a_dst);

	BOOST_REQUIRE_EQUAL(ret,true);

	//! [Usage of meta copy and compare for complex aggregates object]

	a_dst.template get<3>()[1] = std::string("Last string 20");
	ret = meta_compare<aggregate<std::string,std::vector<float>,std::map<size_t,std::string>,std::string[3]>>::meta_compare_f(a_src,a_dst);
	BOOST_REQUIRE_EQUAL(ret,false);

	}


	{

	//! [Usage of meta copy and compare for Point_test]

	typedef Point_test<float> p;

	Point_test<float> p_src;
	Point_test<float> p_dst;

	// fill p_src
	p_src.template get<p::x>() = 1;
	p_src.template get<p::y>() = 567;
	p_src.template get<p::z>() = 341;
	p_src.template get<p::s>() = 5670;
	p_src.template get<p::v>()[0] = 921;
	p_src.template get<p::v>()[1] = 5675;
	p_src.template get<p::v>()[2] = 117;
	p_src.template get<p::t>()[0][0] = 1921;
	p_src.template get<p::t>()[0][1] = 25675;
	p_src.template get<p::t>()[0][2] = 3117;
	p_src.template get<p::t>()[1][0] = 4921;
	p_src.template get<p::t>()[1][1] = 55675;
	p_src.template get<p::t>()[1][2] = 6117;
	p_src.template get<p::t>()[2][0] = 7921;
	p_src.template get<p::t>()[2][1] = 85675;
	p_src.template get<p::t>()[2][2] = 9117;

	meta_copy<Point_test<float>>::meta_copy_(p_src,p_dst);

	bool ret = meta_compare<Point_test<float>>::meta_compare_f(p_src,p_dst);
	BOOST_REQUIRE_EQUAL(ret,true);

	//! [Usage of meta copy and compare for Point_test]

	p_dst.template get<p::t>()[2][2] = 9317;
	ret = meta_compare<Point_test<float>>::meta_compare_f(p_src,p_dst);
	BOOST_REQUIRE_EQUAL(ret,false);

	}
}


BOOST_AUTO_TEST_CASE( meta_copy_test_op )
{
	{
	//! [Usage of meta copy with operation]

	float f_src = 1.0;
	float f_dst = 2.0;

	meta_copy_op<add_,float>::meta_copy_op_(f_src,f_dst);

	BOOST_REQUIRE_EQUAL(f_dst,3.0);

	//! [Usage of meta copy with operation]

	}

	{
	//! [Usage of meta copy with operation for array of primitives]

	float f_src[2][3] = {{1.0,2.5,4.0},{2.5,4.5,9.0}};
	float f_dst[2][3] = {{1.0,2.0,3.0},{1.0,2.0,3.0}};

	meta_copy_op<add_,float[2][3]>::meta_copy_op_(f_src,f_dst);

	BOOST_REQUIRE_EQUAL(f_dst[0][0],2.0);
	BOOST_REQUIRE_EQUAL(f_dst[0][1],4.5);
	BOOST_REQUIRE_EQUAL(f_dst[0][2],7.0);

	BOOST_REQUIRE_EQUAL(f_dst[1][0],3.5);
	BOOST_REQUIRE_EQUAL(f_dst[1][1],6.5);
	BOOST_REQUIRE_EQUAL(f_dst[1][2],12.0);

	//! [Usage of meta copy with operation for array of primitives]

	}

	{
	//! [Usage of meta copy with operation for openfpm aggregates]

	aggregate<float,int,float[3]> agg1;
	aggregate<float,int,float[3]> agg2;

	boost::fusion::at_c<0>(agg1.data) = 1.0;
	boost::fusion::at_c<1>(agg1.data) = 2.0;
	boost::fusion::at_c<2>(agg1.data)[0] = 3.0;
	boost::fusion::at_c<2>(agg1.data)[1] = 4.0;
	boost::fusion::at_c<2>(agg1.data)[2] = 5.0;

	boost::fusion::at_c<0>(agg2.data) = 1.0;
	boost::fusion::at_c<1>(agg2.data) = 2.0;
	boost::fusion::at_c<2>(agg2.data)[0] = 3.0;
	boost::fusion::at_c<2>(agg2.data)[1] = 4.0;
	boost::fusion::at_c<2>(agg2.data)[2] = 5.0;

	meta_copy_op<add_,aggregate<float,int,float[3]>>::meta_copy_op_(agg1,agg2);

	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<0>(agg2.data),2.0);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<1>(agg2.data),4.0);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<2>(agg2.data)[0],6.0);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<2>(agg2.data)[1],8.0);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<2>(agg2.data)[2],10.0);

	//! [Usage of meta copy with operation for openfpm aggregates]


	}

	{
	//! [Usage of meta copy with operation for complex object]

	std::string s_src("Test string");
	std::string s_dst("Test string2");

	meta_copy_op<add_,std::string>::meta_copy_op_(s_src,s_dst);

	BOOST_REQUIRE_EQUAL(s_dst,std::string("Test string2Test string"));

	//! [Usage of meta copy with operation for complex object]

	}


	{

	//! [Usage of meta copy with operation for Point_test]

	typedef Point_test<float> p;

	Point_test<float> p_src;
	Point_test<float> p_dst;

	// fill p_src
	p_src.template get<p::x>() = 1;
	p_src.template get<p::y>() = 567;
	p_src.template get<p::z>() = 341;
	p_src.template get<p::s>() = 5670;
	p_src.template get<p::v>()[0] = 921;
	p_src.template get<p::v>()[1] = 5675;
	p_src.template get<p::v>()[2] = 117;
	p_src.template get<p::t>()[0][0] = 1921;
	p_src.template get<p::t>()[0][1] = 25675;
	p_src.template get<p::t>()[0][2] = 3117;
	p_src.template get<p::t>()[1][0] = 4921;
	p_src.template get<p::t>()[1][1] = 55675;
	p_src.template get<p::t>()[1][2] = 6117;
	p_src.template get<p::t>()[2][0] = 7921;
	p_src.template get<p::t>()[2][1] = 85675;
	p_src.template get<p::t>()[2][2] = 9117;

	// fill p_dst
	p_dst.template get<p::x>() = 2;
	p_dst.template get<p::y>() = 568;
	p_dst.template get<p::z>() = 342;
	p_dst.template get<p::s>() = 5671;
	p_dst.template get<p::v>()[0] = 922;
	p_dst.template get<p::v>()[1] = 5676;
	p_dst.template get<p::v>()[2] = 118;
	p_dst.template get<p::t>()[0][0] = 1922;
	p_dst.template get<p::t>()[0][1] = 25676;
	p_dst.template get<p::t>()[0][2] = 3118;
	p_dst.template get<p::t>()[1][0] = 4922;
	p_dst.template get<p::t>()[1][1] = 55676;
	p_dst.template get<p::t>()[1][2] = 6118;
	p_dst.template get<p::t>()[2][0] = 7922;
	p_dst.template get<p::t>()[2][1] = 85676;
	p_dst.template get<p::t>()[2][2] = 9118;

	meta_copy_op<add_,Point_test<float>>::meta_copy_op_(p_src,p_dst);

	BOOST_REQUIRE_EQUAL(p_dst.template get<p::x>(),3);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::y>(),1135);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::z>(),683);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::v>()[0],1843);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::v>()[1],11351);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::v>()[2],235);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::t>()[0][0],3843);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::t>()[0][1],51351);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::t>()[0][2],6235);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::t>()[1][0],9843);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::t>()[1][1],111351);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::t>()[1][2],12235);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::t>()[2][0],15843);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::t>()[2][1],171351);
	BOOST_REQUIRE_EQUAL(p_dst.template get<p::t>()[2][2],18235);

	//! [Usage of meta copy with operation for complex aggregates object]

	}
}

BOOST_AUTO_TEST_CASE( meta_copy_d_compare_test )
{
	{
	//! [Usage of meta_copy_d for primitives]

	float f_src = 1.0;
	float f_dst;

	meta_copy_d<float,float>::meta_copy_d_(f_src,f_dst);

	BOOST_REQUIRE_EQUAL(f_src,f_dst);

	//! [Usage of meta_copy_d for primitives]

	}

	{
	//! [Usage of meta_copy_d for array of primitives]

	float f_src[2][3] = {{1.0,2.9,4.0},{2.3,4.4,9.0}};
	float f_dst[2][3];

	meta_copy_d<float[2][3],float[2][3]>::meta_copy_d_(f_src,f_dst);

	BOOST_REQUIRE_EQUAL(f_src[0][0],f_dst[0][0]);
	BOOST_REQUIRE_EQUAL(f_src[0][1],f_dst[0][1]);
	BOOST_REQUIRE_EQUAL(f_src[0][2],f_dst[0][2]);
	BOOST_REQUIRE_EQUAL(f_src[1][0],f_dst[1][0]);
	BOOST_REQUIRE_EQUAL(f_src[1][1],f_dst[1][1]);
	BOOST_REQUIRE_EQUAL(f_src[1][2],f_dst[1][2]);

	//! [Usage of meta_copy_d for array of primitives]

	}

	{
	//! [Usage of meta_copy_d for openfpm aggregates]

	aggregate<float,int,float[3]> agg1;
	aggregate<float,int,float[3]> agg2;

	boost::fusion::at_c<0>(agg1.data) = 1.0;
	boost::fusion::at_c<1>(agg1.data) = 2.0;
	boost::fusion::at_c<2>(agg1.data)[0] = 3.0;
	boost::fusion::at_c<2>(agg1.data)[1] = 4.0;
	boost::fusion::at_c<2>(agg1.data)[2] = 5.0;

	meta_copy_d<aggregate<float,int,float[3]>,aggregate<float,int,float[3]>>::meta_copy_d_(agg1,agg2);

	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<0>(agg1.data),boost::fusion::at_c<0>(agg2.data));
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<1>(agg1.data),boost::fusion::at_c<1>(agg2.data));
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<2>(agg1.data)[0],boost::fusion::at_c<2>(agg2.data)[0]);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<2>(agg1.data)[1],boost::fusion::at_c<2>(agg2.data)[1]);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<2>(agg1.data)[2],boost::fusion::at_c<2>(agg2.data)[2]);

	//! [Usage of meta_copy_d for openfpm aggregates]

	}

	{
	//! [Usage of meta_copy_d for complex object]

	std::string s_src("Test string");
	std::string s_dst;

	meta_copy_d<std::string,std::string>::meta_copy_d_(s_src,s_dst);

	BOOST_REQUIRE_EQUAL(s_src,s_dst);

	//! [Usage of meta_copy_d and compare for complex object]

	}

	{
	//! [Usage of meta_copy_d for complex aggregates object]

	aggregate<std::string,std::vector<float>,std::map<size_t,std::string>,std::string[3]> a_src;

	// fill the complex aggregates
	a_src.template get<0>() = std::string("Test string");

	a_src.template get<1>().push_back(5.0);
	a_src.template get<1>().push_back(15.0);
	a_src.template get<1>().push_back(45.0);
	a_src.template get<1>().push_back(7.0);

	a_src.template get<2>()[0] = std::string("Test string 2");
	a_src.template get<2>()[10] = std::string("Test string 3");
	a_src.template get<2>()[9] = std::string("Test string 4");
	a_src.template get<2>()[1] = std::string("Test string 5");

	a_src.template get<3>()[0] = std::string("Last string 9");
	a_src.template get<3>()[1] = std::string("Last string 10");
	a_src.template get<3>()[2] = std::string("Last string 11");

	aggregate<std::string,std::vector<float>,std::map<size_t,std::string>,std::string[3]> a_dst;

	meta_copy_d<aggregate<std::string,std::vector<float>,std::map<size_t,std::string>,std::string[3]>,
                aggregate<std::string,std::vector<float>,std::map<size_t,std::string>,std::string[3]>>::meta_copy_d_(a_src,a_dst);


	BOOST_REQUIRE_EQUAL(a_src.template get<0>(),a_dst.template get<0>());
	BOOST_REQUIRE_EQUAL(a_src.template get<1>()[0],a_dst.template get<1>()[0]);
	BOOST_REQUIRE_EQUAL(a_src.template get<1>()[1],a_dst.template get<1>()[1]);
	BOOST_REQUIRE_EQUAL(a_src.template get<1>()[2],a_dst.template get<1>()[2]);
	BOOST_REQUIRE_EQUAL(a_src.template get<1>()[3],a_dst.template get<1>()[3]);

	BOOST_REQUIRE_EQUAL(a_src.template get<2>()[0],a_dst.template get<2>()[0]);
	BOOST_REQUIRE_EQUAL(a_src.template get<2>()[10],a_dst.template get<2>()[10]);
	BOOST_REQUIRE_EQUAL(a_src.template get<2>()[9],a_dst.template get<2>()[9]);
	BOOST_REQUIRE_EQUAL(a_src.template get<2>()[1],a_dst.template get<2>()[1]);

	BOOST_REQUIRE_EQUAL(a_src.template get<3>()[0],a_dst.template get<3>()[0]);
	BOOST_REQUIRE_EQUAL(a_src.template get<3>()[1],a_dst.template get<3>()[1]);
	BOOST_REQUIRE_EQUAL(a_src.template get<3>()[2],a_dst.template get<3>()[2]);

	//! [Usage of meta_copy_d for complex aggregates object]

	}


	{

	//! [Usage of meta_copy_d and compare for Point_test]

	typedef Point_test<float> p;

	Point_test<float> p_src;
	Point_test<float> p_dst;

	// fill p_src
	p_src.template get<p::x>() = 1;
	p_src.template get<p::y>() = 567;
	p_src.template get<p::z>() = 341;
	p_src.template get<p::s>() = 5670;
	p_src.template get<p::v>()[0] = 921;
	p_src.template get<p::v>()[1] = 5675;
	p_src.template get<p::v>()[2] = 117;
	p_src.template get<p::t>()[0][0] = 1921;
	p_src.template get<p::t>()[0][1] = 25675;
	p_src.template get<p::t>()[0][2] = 3117;
	p_src.template get<p::t>()[1][0] = 4921;
	p_src.template get<p::t>()[1][1] = 55675;
	p_src.template get<p::t>()[1][2] = 6117;
	p_src.template get<p::t>()[2][0] = 7921;
	p_src.template get<p::t>()[2][1] = 85675;
	p_src.template get<p::t>()[2][2] = 9117;

	meta_copy_d<Point_test<float>,Point_test<float>>::meta_copy_d_(p_src,p_dst);

	BOOST_REQUIRE_EQUAL(p_src.template get<p::x>(),p_dst.template get<p::x>());
	BOOST_REQUIRE_EQUAL(p_src.template get<p::y>(),p_dst.template get<p::y>());
	BOOST_REQUIRE_EQUAL(p_src.template get<p::z>(),p_dst.template get<p::z>());
	BOOST_REQUIRE_EQUAL(p_src.template get<p::s>(),p_dst.template get<p::s>());

	BOOST_REQUIRE_EQUAL(p_src.template get<p::v>()[0],p_dst.template get<p::v>()[0]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::v>()[1],p_dst.template get<p::v>()[1]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::v>()[2],p_dst.template get<p::v>()[2]);

	BOOST_REQUIRE_EQUAL(p_src.template get<p::t>()[0][0],p_dst.template get<p::t>()[0][0]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::t>()[0][1],p_dst.template get<p::t>()[0][1]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::t>()[0][2],p_dst.template get<p::t>()[0][2]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::t>()[1][0],p_dst.template get<p::t>()[1][0]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::t>()[1][1],p_dst.template get<p::t>()[1][1]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::t>()[1][2],p_dst.template get<p::t>()[1][2]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::t>()[2][0],p_dst.template get<p::t>()[2][0]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::t>()[2][1],p_dst.template get<p::t>()[2][1]);
	BOOST_REQUIRE_EQUAL(p_src.template get<p::t>()[2][2],p_dst.template get<p::t>()[2][2]);

	//! [Usage of meta_copy_d and compare for Point_test]

	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_UTIL_META_CC_UNIT_TESTS_HPP_ */
