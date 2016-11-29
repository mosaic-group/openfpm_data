/*
 * check_no_pointers.hpp
 *
 *  Created on: Jul 7, 2015
 *      Author: i-bird
 */

#ifndef CHECK_NO_POINTERS_HPP_
#define CHECK_NO_POINTERS_HPP_

#include "common.hpp"

/*! \brief return type of check_no_pointers
 *
 *
 */
enum PNP
{
	POINTERS,
	NO_POINTERS,
	UNKNOWN
};

//! Check if the type T has pointers inside
template<typename T, bool has_pointer=has_noPointers<T>::type::value >
struct check_no_pointers_impl
{
	//! Return true if the structure T has a pointer
	static size_t value()	{return T::noPointers();}
};

//! Check if the type T has pointers inside
template<typename T>
struct check_no_pointers_impl<T,false>
{
	//! Return PNP::UNKNOWN if the structure T does not specify if it has a pointer
	static size_t value()	{return PNP::UNKNOWN;};
};

/*! \brief This class check if the type T has pointers inside
 *
 * It basically check if exist a function bool T::noPointers(), if not return UNKOWN
 * if exist it return the noPointers return status
 *
 * \tparam T type to check
 *
 * ### Example
 *
 * \snippet util_test.hpp Check no pointers in structure
 *
 */
template<typename T>
struct check_no_pointers
{
	static size_t value() {return check_no_pointers_impl<T>::value();}
};




#endif /* CHECK_NO_POINTERS_HPP_ */
