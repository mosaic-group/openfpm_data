/*
 * compare_general.hpp
 *
 *  Created on: Nov 1, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_COMPARE_GENERAL_HPP_
#define OPENFPM_DATA_SRC_UTIL_COMPARE_GENERAL_HPP_


#include "util/common.hpp"
#include "util/util_debug.hpp"
#include "util/for_each_ref.hpp"
#include <boost/mpl/range_c.hpp>

/*! \brief structure to copy aggregates
 *
 * \tparam T type to copy
 *
 */
template<typename T, unsigned int agg=2 * is_openfpm_native<T>::value + std::is_copy_assignable<T>::value>
struct compare_general
{
	/*! \brief Spacialization when there is unknown compare method
	 *
	 * \tparam src source object to copy
	 * \tparam dst destination object
	 *
	 */
	static inline bool compare_general_f(const T & src, const T & dst)
	{
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << demangle(typeid(T).name()) << " does not have an operator== and is not an aggregate or an openfpm native structure, comparation is not possible" << "\n";
		return false;
	}
};

/*! \brief
 *
 *
 */
template<typename T>
struct compare_general<T,1>
{
	/*! \brief compare objects that has an operator== (implicit or explicit)
	 *
	 * \tparam src source object to copy
	 * \tparam dst destination object
	 *
	 */
	static inline bool compare_general_f(const T & src, const T & dst)
	{
		return dst == src;
	}
};

template<typename T>
struct compare_general<T,2>
{
	/*! \brief compare objects that are aggregates
	 *
	 * \tparam src source object to copy
	 * \tparam dst destination object
	 *
	 */
	static inline bool compare_general_f(const T & src, const T & dst)
	{
		compare_aggregate<T> cp(src,dst);

		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::max_prop>>(cp);

		return cp.result();
	}
};

template<typename T>
struct compare_general<T,3>
{
	/*! \brief compare objects that are aggregates but define an operator=
	 *
	 * \tparam src source object to copy
	 * \tparam dst destination object
	 *
	 */
	static inline bool compare_general_f(const T & src, const T & dst)
	{
		return dst == src;
	}
};


#endif /* OPENFPM_DATA_SRC_UTIL_COMPARE_GENERAL_HPP_ */
