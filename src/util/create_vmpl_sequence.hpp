/*
 * mpl_sequence.hpp
 *
 * Set of classes to create an mpl::vector with a filled sequence
 *
 *  Created on: Aug 27, 2014
 *      Author: Pietro Incardona
 */

#ifndef MPL_SEQUENCE_HPP
#define MPL_SEQUENCE_HPP

#include <boost/mpl/int.hpp>
#include <boost/mpl/reverse.hpp>

template<int ...> struct index_tuple_sq{};

/*! \brief Exit condition
 *
 * Exit condition, when c equal to end return true
 *
 * \param c actual position in the counter
 * \param end end of the counter
 *
 */

template<int c,int end>
struct exit_impl_sq : boost::mpl::equal_to<boost::mpl::int_<c>,boost::mpl::int_<end>>
{};

/*! \brief Recursive specialization of to_variadic
 *
 * \param H suppose to be the vector where we are pushing numbers
 * \param c suppose to be the counter
 * \param L suppose to be the end of the cycle
 * \param exit, when true it say to terminate the sequence
 *
 */
template<int c , int end, bool exit,int ... vars>
struct to_variadic_impl
{
   typedef typename exit_impl_sq<c,end>::type exit_;
   typedef typename to_variadic_impl<c+1,end,exit_::value,vars ..., c>::type type;
};


//! Terminator of to_variadic
template<int c, int end, int ... vars>
struct to_variadic_impl<c,end,true, vars ...>
{
   typedef index_tuple_sq<vars ...> type;
};


/*!
*
* It create a mpl::vector sequence of integer from N to M
*
* usage:
*
* to_int_sequence<3,5>::type
*
* is mapped to
*
* boost::mpl::vector<int_<3>,int_<4>,int<5>>
*
* \param N metafunction (is a structure with typedef type)
* \param M boost::fusion::vector
*
*/

template<unsigned int N, unsigned int M>
struct to_int_sequence
{
	//! end condition
	typedef typename exit_impl_sq<N,M>::type exit_;

	//! generate the boost::fusion::vector apply H on each term
	typedef typename to_variadic_impl<N+1,M,exit_::value,N>::type type;
};


#endif
