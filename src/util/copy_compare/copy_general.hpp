/*
 * aggregate_copy.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_COPY_GENERAL_HPP_
#define OPENFPM_DATA_SRC_UTIL_COPY_GENERAL_HPP_

#include "util/common.hpp"
#include "util/util_debug.hpp"
#include "copy_compare_aggregates.hpp"
#include "util/for_each_ref.hpp"
#include <boost/mpl/range_c.hpp>
#include <iostream>

/*! \brief structure to copy aggregates
 *
 * \tparam T type to copy
 *
 */
template<typename T, unsigned int agg=2 * is_openfpm_native<T>::value + std::is_copy_assignable<T>::value>
struct copy_general
{
	/*! \brief Spacialization when there is unknown copy method
	 *
	 * \tparam src source object to copy
	 * \tparam dst destination object
	 *
	 */
	inline copy_general(const T & src, T & dst)
	{
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << demangle(typeid(T).name()) << " does not have an operator= and is not an aggregate or an openfpm native structure, copy is not possible" << "\n";
	}
};

/*! \brief
 *
 *
 */
template<typename T>
struct copy_general<T,1>
{
	/*! \brief copy objects that has an operator= (implicit or explicit)
	 *
	 * \tparam src source object to copy
	 * \tparam dst destination object
	 *
	 */
	inline copy_general(const T & src, T & dst)
	{
		dst = src;
	}
};

template<typename T>
struct copy_general<T,2>
{
	/*! \brief copy objects that are aggregates
	 *
	 * \tparam src source object to copy
	 * \tparam dst destination object
	 *
	 */
	inline copy_general(const T & src, T & dst)
	{
		copy_aggregate<T> cp(src,dst);

		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::max_prop>>(cp);
	}
};

template<typename T>
struct copy_general<T,3>
{
	/*! \brief copy objects that are aggregates but define an operator=
	 *
	 * \tparam src source object to copy
	 * \tparam dst destination object
	 *
	 */
	inline copy_general(const T & src, T & dst)
	{
		dst = src;
	}
};

#endif /* OPENFPM_DATA_SRC_UTIL_COPY_GENERAL_HPP_ */
