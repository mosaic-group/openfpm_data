/*
 * has_max_prop.hpp
 *
 *  Created on: Jul 27, 2016
 *      Author: yaroslav
 */

#ifndef OPENFPM_DATA_SRC_PACKER_UNPACKER_HAS_MAX_PROP_HPP_
#define OPENFPM_DATA_SRC_PACKER_UNPACKER_HAS_MAX_PROP_HPP_

#include "util/common.hpp"

template<typename T, typename Sfinae = void>
struct has_max_prop_nn: std::false_type {};

/*! \brief has_max_prop check if a type has defined a member max_prop
 *
 * ### Example
 *
 * \snippet util_test.hpp Check has_max_prop
 *
 * return true if T::value_type::max_prop is a valid type
 *
 */
template<typename T>
struct has_max_prop_nn<T, typename Void<decltype( T::max_prop )>::type> : std::true_type
{};

template<typename T, bool has_max_prop>
struct max_prop_nn
{
	enum
	{
		number = T::max_prop
	};
};

template<typename T>
struct max_prop_nn<T,false>
{
	enum
	{
		number = 0
	};
};


template<typename T, bool hvt>
struct has_max_prop
{
	typedef has_max_prop<typename T::value_type,has_value_type< typename T::value_type >::value> hmp;

	enum
	{
		value = hmp::value
	};

	enum
	{
		number = hmp::number
	};
};

template<typename T>
struct has_max_prop<T,false>
{
	enum
	{
		value = has_max_prop_nn<T>::value
	};

	enum
	{
		number = max_prop_nn<T,has_max_prop_nn<T>::value>::number
	};
};




#endif /* OPENFPM_DATA_SRC_PACKER_UNPACKER_HAS_MAX_PROP_HPP_ */
