/*
 * util.hpp
 *
 *  Created on: Jul 16, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_UTIL_HPP_
#define SRC_GRID_UTIL_HPP_

#include "util/common.hpp"

template<typename T, typename Sfinae = void>
struct is_grid: std::false_type {};


/*! \brief is_grid check if the type is a grid
 *
 * ### Example
 *
 * \snippet util.hpp Check is_grid
 *
 * return true if T is a grid
 *
 */
template<typename T>
struct is_grid<T, typename Void< typename T::yes_i_am_grid>::type> : std::true_type
{};


#endif /* SRC_GRID_UTIL_HPP_ */
