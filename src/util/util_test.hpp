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
#include "util/ct_array.hpp"
#include "Vector/map_vector.hpp"
#include "common.hpp"
#include "check_no_pointers.hpp"
#include "Grid/util.hpp"
#include "data_type/scalar.hpp"
#include "util/convert.hpp"

//! [Check has_posMask struct definition]

struct test_has_posMask
{
	float data;

	static constexpr bool stag_mask[] = {true,false,true};
};

struct test_no_has_posMask
{
	float data;
};

//! [Check has_posMask struct definition]

//! [Declaration of struct with attributes and without]

struct test_has_attributes
{
	struct attributes
	{
		static const std::string name[2];
	};
};

const std::string test_has_attributes::attributes::name[]={"attributes1","attributes2"};

struct test_no_attributes
{
};

//! [Declaration of struct with attributes and without]

BOOST_AUTO_TEST_SUITE( util_test )

BOOST_AUTO_TEST_CASE( object_prop_copy )
{
	{
	//! [object copy example]
	typedef Point_test<float> p;
	typedef Point_test<float>::type vboost;
	typedef object_creator<Point_test<float>::type,0,1,4>::type vboost_red;

	object<vboost> src;
	object<vboost_red> dst;

	// fill the source properties x,y,v = 0,1,4 with data

	boost::fusion::at_c<p::x>(src.data) = 1.0;
	boost::fusion::at_c<p::y>(src.data) = 2.0;

	for (size_t i = 0 ; i < 3 ;  i++)
		boost::fusion::at_c<p::v>(src.data)[i] = i + 5.0;

	// copy from src to dst
	object_si_d<object<vboost>,object<vboost_red>,OBJ_NORMAL,0,1,4>(src,dst);

	// Check the result
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<0>(dst.data),1.0);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<1>(dst.data),2.0);

	for (size_t i = 0 ; i < 3 ;  i++)
		BOOST_REQUIRE_EQUAL(boost::fusion::at_c<2>(dst.data)[i],i + 5.0);

	//! [object copy example]

	//! [object copy encap example]

	typedef encapc<1,Point_test<float>,openfpm::vector<Point_test<float>>::memory_conf> encap_src;
	typedef encapc<1,object<vboost_red>,openfpm::vector<object<vboost_red>>::memory_conf> encap_dst;

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

	object_si_d<encap_src,encap_dst,OBJ_ENCAP,0,1,4>(src_e,dst_e);

	BOOST_REQUIRE_EQUAL(v_point_red.get(0).template get<0>(),1.0);
	BOOST_REQUIRE_EQUAL(v_point_red.get(0).template get<1>(),2.0);

	for (size_t i = 0 ; i < 3 ;  i++)
		BOOST_REQUIRE_EQUAL(v_point_red.get(0).template get<2>()[i],i + 5.0);

	//! [object copy encap example]
	}

	{
	// Object write test

	//! [object write example]
	typedef Point_test<float> p;
	typedef Point_test<float>::type vboost;
	typedef object_creator<Point_test<float>::type,0,1,4>::type vboost_red;

	object<vboost> dst;
	object<vboost_red> src;

	// fill the source properties 0,1,2 with data

	boost::fusion::at_c<0>(src.data) = 1.0;
	boost::fusion::at_c<1>(src.data) = 2.0;

	for (size_t i = 0 ; i < 3 ;  i++)
		boost::fusion::at_c<2>(src.data)[i] = i + 5.0;

	// copy from src to dst
	object_s_di<object<vboost_red>,object<vboost>,OBJ_NORMAL,0,1,4>(src,dst);

	// Check the result
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<p::x>(dst.data),1.0);
	BOOST_REQUIRE_EQUAL(boost::fusion::at_c<p::y>(dst.data),2.0);

	for (size_t i = 0 ; i < 3 ;  i++)
		BOOST_REQUIRE_EQUAL(boost::fusion::at_c<p::v>(dst.data)[i],i + 5.0);

	//! [object write example]

	//! [object write encap example]

	typedef encapc<1,Point_test<float>,openfpm::vector<Point_test<float>>::memory_conf> encap_dst;
	typedef encapc<1,object<vboost_red>,openfpm::vector<object<vboost_red>>::memory_conf> encap_src;

	openfpm::vector<p> v_point;
	openfpm::vector<object<vboost_red>> v_point_red;

	v_point.resize(2);
	v_point_red.resize(2);

	v_point_red.template get<0>(0) = 11.0;
	v_point_red.template get<1>(0) = 12.0;

	for (size_t i = 0 ; i < 3 ;  i++)
		v_point_red.template get<2>(0)[i] = i + 15.0;

	auto dst_e = v_point.get(0);
	auto src_e = v_point_red.get(0);

	object_s_di<encap_src,encap_dst,OBJ_ENCAP,0,1,4>(src_e,dst_e);

	BOOST_REQUIRE_EQUAL(v_point.get(0).template get<p::x>(),11.0);
	BOOST_REQUIRE_EQUAL(v_point.get(0).template get<p::y>(),12.0);

	for (size_t i = 0 ; i < 3 ;  i++)
		BOOST_REQUIRE_EQUAL(v_point.get(0).template get<p::v>()[i],i + 15.0);

	//! [object write encap example]
	}

	{
	//! [object creator check for no pointers]
	struct no_method_pointer
	{
		float a;
		int b;
	};

	struct no_pointer
	{
		double c;
		int d;

		static bool noPointers() {return true;}
	};

	struct with_pointer
	{
		double * c;
		int d;

		static bool noPointers() {return false;}
	};

	typedef boost::fusion::vector<int,float,double,no_method_pointer,no_pointer,with_pointer,no_method_pointer> v;

	int val = boost::mpl::size< noPointers_sequence<v,0,2>::type >::value;
	BOOST_REQUIRE_EQUAL(val,0);

	val = boost::mpl::size< noPointers_sequence<v,0,2,4,5>::type >::value;
	BOOST_REQUIRE_EQUAL(val,0);

	val = boost::mpl::size< noPointers_sequence<v,0,1,2,3,4,5,6>::type >::value;
	BOOST_REQUIRE_EQUAL(val,2);

	typedef boost::fusion::vector<int *,float &,double *,no_method_pointer *,no_pointer &,with_pointer *,no_method_pointer &> vp;

	val = boost::mpl::size< noPointers_sequence<vp,0,1,2,3,4,5,6>::type >::value;
	BOOST_REQUIRE_EQUAL(val,2);

	//! [object creator check for no pointers]
	}

	// Check the the object respect the noPointers construction
	{

	struct no_method_pointer
	{
		float a;
		int b;
	};

	struct no_pointer
	{
		double c;
		int d;

		static bool noPointers() {return true;}
	};

	struct with_pointer
	{
		double * c;
		int d;

		static bool noPointers() {return false;}
	};

	typedef boost::fusion::vector<int,float,double> v;
	int val = object<v>::noPointers();
	BOOST_REQUIRE_EQUAL(val,PNP::NO_POINTERS);

	typedef boost::fusion::vector<int,float,double,no_pointer> va;
	val = object<va>::noPointers();
	BOOST_REQUIRE_EQUAL(val,PNP::NO_POINTERS);

	typedef boost::fusion::vector<int,float,double,no_pointer,with_pointer> vb;
	val = object<vb>::noPointers();
	BOOST_REQUIRE_EQUAL(val,PNP::POINTERS);

	typedef boost::fusion::vector<int,float,double,no_pointer,double *> vc;
	val = object<vc>::noPointers();
	BOOST_REQUIRE_EQUAL(val,PNP::POINTERS);

	typedef boost::fusion::vector<int,float,double,no_pointer,double &> vd;
	val = object<vd>::noPointers();
	BOOST_REQUIRE_EQUAL(val,PNP::POINTERS);

	typedef boost::fusion::vector<int,float,double,no_pointer,double[3]> ve;
	val = object<ve>::noPointers();
	BOOST_REQUIRE_EQUAL(val,PNP::NO_POINTERS);

	typedef boost::fusion::vector<int,float,double,no_pointer,no_pointer *> vf;
	val = object<vf>::noPointers();
	BOOST_REQUIRE_EQUAL(val,PNP::POINTERS);

	typedef boost::fusion::vector<int,float,double,no_pointer,no_pointer &> vg;
	val = object<vg>::noPointers();
	BOOST_REQUIRE_EQUAL(val,PNP::POINTERS);

	typedef boost::fusion::vector<int,float,double,no_pointer,no_pointer[3]> vh;
	val = object<vh>::noPointers();
	BOOST_REQUIRE_EQUAL(val,PNP::NO_POINTERS);
	}
}

//! [Metafunction definition]

template<size_t index, size_t N> struct MetaFunc {
   enum { value = index + N };
};

//! [Metafunction definition]

BOOST_AUTO_TEST_CASE( generate_array )
{
	{
	//! [compile time array]
	const size_t count = 5;
	typedef typename ::generate_array<size_t,count, MetaFunc>::result ct_test;

	// ct_test::data is equivalent to const size_t [5] = {5,6,7,8,9}

	for (size_t i = 0 ; i < count; ++i)
	{
		const size_t ct_val = ct_test::data[i];
		BOOST_REQUIRE_EQUAL(ct_val,count+i);
	}
	//! [compile time array]
	}

	// check constexpr compile time array as template parameters

#ifdef COVERTY_SCAN

	{
	//! [constexpr array]
	const size_t count = 5;
	typedef typename ::generate_array_constexpr<size_t,count, MetaFunc>::result ct_test_ce;

	// ct_test_ce::data is equivalent to constexpr size_t [5] = {5,6,7,8,9}

	const size_t ct_calc = MetaFunc<ct_test_ce::data[0],ct_test_ce::data[1]>::value;
	BOOST_REQUIRE_EQUAL(ct_calc,11);
	//! [constexpr array]
	}

#endif
}

BOOST_AUTO_TEST_CASE( check_templates_util_function )
{
	{
		{
		//! [Check no pointers]

		struct test_no_ptr
		{
			static bool noPointers()	{return true;}
		};

		struct test_ptr
		{
			static bool noPointers()	{return false;}
		};

		struct test_unknown
		{
		};

		BOOST_REQUIRE_EQUAL(has_noPointers<test_no_ptr>::type::value,true);
		BOOST_REQUIRE_EQUAL(has_noPointers<test_ptr>::type::value,true);
		BOOST_REQUIRE_EQUAL(has_noPointers<test_unknown>::type::value,false);

		//! [Check no pointers]

		}

		{
		//! [Declaration of an openfpm native structure]

		int val = is_openfpm_native<Point_test<float>>::value;
		BOOST_REQUIRE_EQUAL(val, true);
		val = is_openfpm_native<float>::value;
		BOOST_REQUIRE_EQUAL(val, false);

		//! [Declaration of an openfpm native structure]
		}

		{
		//! [Check is_typedef_and_data_same]

		struct test_typedef_same_data
		{
			typedef boost::fusion::vector<float,double,float[3]> type;

			type data;
		};

		struct test_typedef_not_same_data
		{
			typedef boost::fusion::vector<float,double,float[3]> type;

			boost::fusion::vector<float,double> data;
		};

		int val = is_typedef_and_data_same<true,test_typedef_same_data>::value;
		BOOST_REQUIRE_EQUAL(val, true);
		val = is_typedef_and_data_same<true,test_typedef_not_same_data>::value;
		BOOST_REQUIRE_EQUAL(val, false);

		//! [Check is_typedef_and_data_same]
		}

		{
		//! [Check has_data]

		struct test_has_data
		{
			float data;
		};

		struct test_no_has_data
		{
		};

		int val = has_data<test_has_data>::value;
		BOOST_REQUIRE_EQUAL(val, true);
		val = has_data<test_no_has_data>::value;
		BOOST_REQUIRE_EQUAL(val, false);

		//! [Check has_data]
		}

		{
		//! [Check has_posMask]

		int val = has_posMask<test_has_posMask>::value;
		BOOST_REQUIRE_EQUAL(val, true);
		val = has_posMask<test_no_has_posMask>::value;
		BOOST_REQUIRE_EQUAL(val, false);

		//! [Check has_posMask]
		}

		{
		//! [Check has_typedef_type]

		struct test_has_typedef
		{
			typedef float type;
		};

		struct test_no_has_data
		{
		};

		int val = has_typedef_type<test_has_typedef>::value;
		BOOST_REQUIRE_EQUAL(val, true);
		val = has_typedef_type<test_no_has_data>::value;
		BOOST_REQUIRE_EQUAL(val, false);

		//! [Check has_typedef_type]
		}

		{
		//! [Check has_attributes]

		int val = has_attributes<test_has_attributes>::value;
		BOOST_REQUIRE_EQUAL(val, true);
		val = has_attributes<test_no_attributes>::value;
		BOOST_REQUIRE_EQUAL(val, false);

		//! [Check has_typedef_type]
		}

		//! [Check no pointers in structure]
		struct test_no_ptr
		{
			static bool noPointers()	{return PNP::NO_POINTERS;}
		};

		struct test_ptr
		{
			static bool noPointers()	{return PNP::POINTERS;}
		};

		struct test_unknown
		{
		};

		int val = check_no_pointers<test_no_ptr>::value();
		BOOST_REQUIRE_EQUAL(val,PNP::NO_POINTERS);
		val = check_no_pointers<test_ptr>::value();
		BOOST_REQUIRE_EQUAL(val,PNP::POINTERS);
		val = check_no_pointers_impl<test_unknown,false>::value();
		BOOST_REQUIRE_EQUAL(val,PNP::UNKNOWN);

		//! [Check no pointers in structure]

		{
		//! [Check is_grid]

		struct stub_object
		{
			float a;
			double b;
		};

		bool val = is_grid< grid_cpu<2,scalar<float> > >::value;
		BOOST_REQUIRE_EQUAL( val ,true);
		val = is_grid< grid_cpu<3,object< boost::fusion::vector<float,double> > > >::value;
		BOOST_REQUIRE_EQUAL( val , true);
		val = is_grid< grid_cpu<4,Point_test<float>> >::value;
		BOOST_REQUIRE_EQUAL( val ,true);

		val = is_grid< grid_gpu<2,scalar<float> > >::value;
		BOOST_REQUIRE_EQUAL(val ,true);
		val = is_grid< grid_gpu<3,object< boost::fusion::vector<float,double>>> >::value;
		BOOST_REQUIRE_EQUAL(val, true);
		val = is_grid< grid_gpu<4,Point_test<float>> >::value;
		BOOST_REQUIRE_EQUAL(val,true);

		val = is_grid< float >::value;
		BOOST_REQUIRE_EQUAL( val, false);
		val = is_grid< stub_object > ::value;
		BOOST_REQUIRE_EQUAL( val, false);

		//! [Check is_grid]
		}

		{
		//! [Check is_vector]

		struct stub_object
		{
			float a;
			double b;
		};

		bool val = is_vector< openfpm::vector<scalar<float> > >::value;
		BOOST_REQUIRE_EQUAL( val ,true);
		val = is_vector< openfpm::vector<object< boost::fusion::vector<float,double> > > >::value;
		BOOST_REQUIRE_EQUAL( val , true);
		val = is_vector< openfpm::vector<Point_test<float>> >::value;
		BOOST_REQUIRE_EQUAL( val ,true);

		val = is_vector< openfpm::vector<scalar<float> > >::value;
		BOOST_REQUIRE_EQUAL(val ,true);
		val = is_vector< openfpm::vector<object< boost::fusion::vector<float,double>>> >::value;
		BOOST_REQUIRE_EQUAL(val, true);
		val = is_vector< openfpm::vector<Point_test<float>> >::value;
		BOOST_REQUIRE_EQUAL(val,true);

		val = is_vector< float >::value;
		BOOST_REQUIRE_EQUAL( val, false);
		val = is_vector< stub_object > ::value;
		BOOST_REQUIRE_EQUAL( val, false);

		//! [Check is_vector]
		}

		{
		//! [Check is_encap]

		struct stub_object
		{
			float a;
			double b;
		};

		bool val = is_encap< encapc<2,scalar<float>,memory_traits_lin<scalar<float>>> >::value;
		BOOST_REQUIRE_EQUAL( val ,true);
		val = is_encap< encapc<3,object<boost::fusion::vector<float,double> >,memory_traits_lin<object< boost::fusion::vector<float,double>>>>  >::value;
		BOOST_REQUIRE_EQUAL( val , true);
		val = is_encap< encapc<4,Point_test<float>,memory_traits_lin<Point_test<float>>> >::value;
		BOOST_REQUIRE_EQUAL( val ,true);

		val = is_encap< encapg<2,scalar<float>,memory_traits_lin<scalar<float>>> >::value;
		BOOST_REQUIRE_EQUAL( val ,true);
		val = is_encap< encapg<3,object<boost::fusion::vector<float,double> >,memory_traits_lin<object< boost::fusion::vector<float,double>>>>  >::value;
		BOOST_REQUIRE_EQUAL( val , true);
		val = is_encap< encapg<4,Point_test<float>,memory_traits_lin<Point_test<float>>> >::value;
		BOOST_REQUIRE_EQUAL( val ,true);

		val = is_encap< float >::value;
		BOOST_REQUIRE_EQUAL( val, false);
		val = is_encap< stub_object > ::value;
		BOOST_REQUIRE_EQUAL( val, false);

		//! [Check is_vector]
		}
	}
}

BOOST_AUTO_TEST_CASE( check_convert_function )
{
	{
	//! [Convert combination to Point]

	comb<2> c;
	c.c[0] = 1;
	c.c[1] = -1;

	Point<2,size_t> p = toPoint<2,size_t>::convert(c);

	BOOST_REQUIRE_EQUAL(p.get(0),1);
	BOOST_REQUIRE_EQUAL(p.get(1),-1);

	//! [Convert combination to Point]
	}

	{
	//! [Convert combination to Point3]

	comb<3> c;
	c.c[0] = 1;
	c.c[1] = -1;
	c.c[2] = 0;

	Point<3,float> p = toPoint<3,float>::convert(c);

	BOOST_REQUIRE_EQUAL(p.get(0),1);
	BOOST_REQUIRE_EQUAL(p.get(1),-1);
	BOOST_REQUIRE_EQUAL(p.get(2),0);

	//! [Convert combination to Point3]
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* UTIL_TEST_HPP_ */
