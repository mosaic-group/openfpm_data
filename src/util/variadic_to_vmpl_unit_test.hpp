/*
 * variadic_to_vmpl_unit_test.hpp
 *
 *  Created on: Aug 19, 2015
 *      Author: i-bird
 */

#ifndef SRC_UTIL_VARIADIC_TO_VMPL_UNIT_TEST_HPP_
#define SRC_UTIL_VARIADIC_TO_VMPL_UNIT_TEST_HPP_

#include "util/variadic_to_vmpl.hpp"
#include "data_type/scalar.hpp"
#include "util/util_debug.hpp"
#include <typeinfo>

//! [v_transform metafunction]
template <typename T>
struct F
{
	typedef scalar<T> type;
};
//! [v_transform_two metafunction]

//! [v_transform_two metafunction]
template <typename arg0, typename T>
struct Ftwo
{
	typedef scalar<T> type;
};
//! [v_transform_two metafunction]

BOOST_AUTO_TEST_CASE( variadic_to_vmpl_test)
{
	{
	//! [v_transform usage]

	typedef boost::mpl::vector<float,float,float[3]> bfv;

	// tbvf is boost::fusion::vector<scalar<float>,scalar<float>,scalar<float[3]>>
	typedef v_transform<F,bfv>::type tbfv;

	bool val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<0>>::type,scalar<float>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<1>>::type,scalar<float>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<2>>::type,scalar<float[3]>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	//! [v_transform usage]
	}

	{
	//! [v_transform_two usage]

	typedef boost::mpl::vector<float,float,float[3]> bfv;

	// tbvf is boost::fusion::vector<scalar<float>,scalar<float>,scalar<float[3]>>
	typedef v_transform_two<Ftwo,float,bfv>::type tbfv;

	bool val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<0>>::type,scalar<float>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<1>>::type,scalar<float>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<tbfv,boost::mpl::int_<2>>::type,scalar<float[3]>>::value;
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

	val = std::is_same<boost::mpl::at<bfv,boost::mpl::int_<3>>::type,boost::mpl::int_<5>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	val = std::is_same<boost::mpl::at<bfv,boost::mpl::int_<4>>::type,boost::mpl::int_<9>>::value;
	BOOST_REQUIRE_EQUAL(val,true);

	//! [to_boost_vmpl usage]
	}
}


#endif /* SRC_UTIL_VARIADIC_TO_VMPL_UNIT_TEST_HPP_ */
