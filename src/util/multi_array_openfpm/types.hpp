/*
 * types.hpp
 *
 *  Created on: Jun 30, 2018
 *      Author: i-bird
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <boost/mpl/at.hpp>
#include <boost/mpl/size.hpp>

//
// types.hpp - supply types that are needed by several headers
//
#include <cstddef>

#ifdef MULTI_ARRAY_INT_INDEX

namespace openfpm {
namespace detail {
namespace multi_array{

// needed typedefs
typedef int size_type;
typedef int index;

} // namespace multi_array
} // namespace detail
} // namespace openfpm

#else


namespace openfpm {
namespace detail {
namespace multi_array{

// needed typedefs
typedef std::size_t size_type;
typedef std::ptrdiff_t index;

} // namespace multi_array
} // namespace detail
} // namespace openfpm

#endif

template<typename vector, unsigned int p, bool = (p < boost::mpl::size<vector>::type::value) >
struct at_impl
{
	 typedef boost::mpl::int_<0> type;
};

template<typename vector, unsigned int p>
struct at_impl<vector,p,true>
{
	 typedef boost::mpl::int_<boost::mpl::at<vector,boost::mpl::int_<p>>::type::value> type;
};

#endif /* TYPES_HPP_ */
