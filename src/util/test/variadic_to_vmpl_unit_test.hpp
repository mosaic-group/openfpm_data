/*
 * variadic_to_vmpl_unit_test.hpp
 *
 *  Created on: Aug 19, 2015
 *      Author: i-bird
 */

#ifndef SRC_UTIL_VARIADIC_TO_VMPL_UNIT_TEST_HPP_
#define SRC_UTIL_VARIADIC_TO_VMPL_UNIT_TEST_HPP_

#include "util/variadic_to_vmpl.hpp"
#include "util/util_debug.hpp"
#include <typeinfo>

//! [v_transform metafunction]
template <typename T>
struct F
{
	//! meta-function implementation
	typedef aggregate<T> type;
};
//! [v_transform metafunction]

//! [v_transform_two metafunction]
template <typename arg0, typename T>
struct Ftwo
{
	//! meta-function implementation
	typedef aggregate<T> type;
};
//! [v_transform_two metafunction]

BOOST_AUTO_TEST_CASE( variadic_to_vmpl_test)
{
	{
	//! [v_transform usage]

	typedef boost::mpl::vector<float,float,float[3]> bfv;

	// tbvf is boost::fusion::vector<scalar<float>,scalar<float>,scalar<float[3]>>
	typedef v_transform<F,bfv>::type tbfv;

	bool val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<0>>::type,aggregate<float>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<1>>::type,aggregate<float>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<2>>::type,aggregate<float[3]>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	//! [v_transform usage]
	}

	{
	//! [v_transform_two usage]

	typedef boost::mpl::vector<float,float,float[3]> bfv;

	// tbvf is boost::fusion::vector<scalar<float>,scalar<float>,scalar<float[3]>>
	typedef v_transform_two<Ftwo,float,bfv>::type tbfv;

	bool val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<0>>::type,aggregate<float>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<1>>::type,aggregate<float>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<2>>::type,aggregate<float[3]>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	//! [v_transform_two usage]
	}

	{
	//! [to_boost_vmpl usage]

	typedef to_boost_vmpl<1,4,5,9>::type bfv;

	bool val = std::is_same<boost::mpl::at<bfv,boost::mpl::int_<0>>::type,boost::mpl::int_<1>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<bfv,boost::mpl::int_<1>>::type,boost::mpl::int_<4>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<bfv,boost::mpl::int_<2>>::type,boost::mpl::int_<5>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<bfv,boost::mpl::int_<3>>::type,boost::mpl::int_<9>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	//! [to_boost_vmpl usage]
	}

	{
	//! [vmpl_sum_constant usage]

	typedef to_boost_vmpl<1,4,5,9>::type bfv;

	typedef vmpl_sum_constant<5,bfv>::type vsc;

	BOOST_REQUIRE_EQUAL(boost::mpl::size<vsc>::type::value,4);

	bool val = std::is_same<boost::mpl::at<vsc,boost::mpl::int_<0>>::type,boost::mpl::int_<6>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<vsc,boost::mpl::int_<1>>::type,boost::mpl::int_<9>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<vsc,boost::mpl::int_<2>>::type,boost::mpl::int_<10>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<vsc,boost::mpl::int_<3>>::type,boost::mpl::int_<14>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	//! [vmpl_sum_constant usage]
	}
}

BOOST_AUTO_TEST_CASE( lin_vmpl_test )
{
	typedef boost::mpl::vector<boost::mpl::int_<16>,boost::mpl::int_<17>,boost::mpl::int_<18>> vector;

	typedef boost::mpl::vector<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3>> offset;

	int lino = Lin_vmpl_off<vector,offset>(0,0,0);
	int lin = Lin_vmpl<vector>(0,0,0);

	BOOST_REQUIRE_EQUAL(lino,1+2*16+3*16*17);
	BOOST_REQUIRE_EQUAL(lin,0);

	lino = Lin_vmpl_off<vector,offset>(0,1,0);
	lin = Lin_vmpl<vector>(0,1,0);

	BOOST_REQUIRE_EQUAL(lino,1+3*16+3*16*17);
	BOOST_REQUIRE_EQUAL(lin,16);

	lino = Lin_vmpl_off<vector,offset>(0,0,1);
	lin = Lin_vmpl<vector>(0,0,1);

	BOOST_REQUIRE_EQUAL(lino,1+2*16+4*16*17);
	BOOST_REQUIRE_EQUAL(lin,16*17);
}

#endif /* SRC_UTIL_VARIADIC_TO_VMPL_UNIT_TEST_HPP_ */
