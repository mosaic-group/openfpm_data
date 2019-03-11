/*
 * util.hpp
 *
 *  Created on: Jul 16, 2015
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_UTIL_HPP_
#define SRC_VECTOR_UTIL_HPP_

#include "util/common.hpp"

template<typename T, typename Sfinae = void>
struct is_vector: std::false_type {};


/*! \brief is_grid check if the type is a vector
 *
 * ### Example
 *
 * \snippet util.hpp Check is_vector
 *
 * return true if T is a vector
 *
 */
template<typename T>
struct is_vector<T, typename Void< typename T::yes_i_am_vector>::type > : std::true_type
{};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, typename Sfinae = void>
struct is_vector_native: std::false_type {};


/*! \brief is_grid check if the type is a vector
 *
 * ### Example
 *
 * \snippet util.hpp Check is_vector
 *
 * return true if T is a vector
 *
 */
template<typename T>
struct is_vector_native<T, typename Void< typename T::yes_i_am_vector_native>::type > : std::true_type
{};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, typename Sfinae = void>
struct is_vector_dist: std::false_type {};


/*! \brief is_grid check if the type is a vector
 *
 * ### Example
 *
 * \snippet util.hpp Check is_vector
 *
 * return true if T is a vector
 *
 */
template<typename T>
struct is_vector_dist<T, typename Void< typename T::yes_i_am_vector_dist>::type > : std::true_type
{};

#endif /* SRC_VECTOR_UTIL_HPP_ */
