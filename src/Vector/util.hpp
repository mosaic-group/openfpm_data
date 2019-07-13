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

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*! \brief Check this is a gpu or cpu type cell-list
 *
 */
template<typename T, typename Sfinae = void>
struct is_gpu_celllist: std::false_type {};


template<typename T>
struct is_gpu_celllist<T, typename Void<typename T::yes_is_gpu_celllist>::type> : std::true_type
{};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*! \brief Check this is a gpu or cpu type cell-list
 *
 */
template<typename T, typename Sfinae = void>
struct is_gpu_ker_celllist: std::false_type {};


template<typename T>
struct is_gpu_ker_celllist<T, typename Void<typename T::yes_is_gpu_ker_celllist>::type> : std::true_type
{};

// structure to check the device pointer

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. It check if the
 * pointer ptr match one of the pointer properties
 *
 */
template<typename data_type>
struct check_device_ptr
{
	//! pointer to check
	void * ptr;

	//! Data to check
	data_type & data;

	mutable int prp;

	mutable bool result;

	/*! \brief constructor
	 *
	 * \param ptr pointer to check
	 * \param data data structure
	 *
	 */
	inline check_device_ptr(void * ptr, data_type & data)
	:ptr(ptr),data(data),result(false)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		if (data.template getPointer<T::value>() == ptr)
		{
			prp = T::value;
			result = true;
		}
	}
};

#endif /* SRC_VECTOR_UTIL_HPP_ */
