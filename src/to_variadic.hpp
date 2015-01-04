/*
 * to_variadic.hpp
 *
 * Set of classes to convert from a boost::mpl::vector into an S variadic template
 * class appling a metafunction F on each element
 *
 * boost::mpl::vector<T,A,B> is converted into
 *
 * S<F<T>,F<A>,F<B>>
 *
 * \see to_variadic
 *
 *  Created on: Aug 27, 2014
 *      Author: Pietro Incardona
 */

#ifndef TO_VARIADIC_HPP
#define TO_VARIADIC_HPP

#include <boost/mpl/int.hpp>
#include <boost/mpl/reverse.hpp>
#include <boost/mpl/vector.hpp>

/*! \brief Exit condition
 *
 * Exit condition, when the distance between F and L is 0 return true
 *
 * \param F suppose to be the end of the mpl::vector
 * \param L suppose to be the actual position of the mpl::vector
 *
 */

template<typename F,typename L>
struct exit_impl : boost::mpl::equal_to<typename boost::mpl::distance<F,L>::type,boost::mpl::int_<0>>
{};

/*! \brief Recursive specialization of to_variadic
 *
 * \param F suppose to be the original end of boost::mpl::vector
 * \param L suppose to be the actual position of the boost::mpl::vector
 * \param exit, when true it say to terminate the sequence
 *
 */
template<template<typename> class H, typename F,typename L, bool exit,typename ...Args>
struct to_variadic_impl
{
   typedef typename boost::mpl::deref<F>::type front_;
   typedef typename boost::mpl::next<F>::type next_;
   typedef typename exit_impl<next_,L>::type exit_;
   typedef typename to_variadic_impl<H,next_,L,exit_::value,typename H<front_>::type,Args...>::type type;
};


//! Terminator of to_variadic
template<template<typename> class H,typename F,typename L,typename ...Args>
struct to_variadic_impl<H,F,L,true,Args...>
{
   typedef boost::fusion::vector<Args...> type;
};

template<typename Seq>
struct seq_traits_impl
{
   typedef typename boost::mpl::begin<Seq>::type first_;
   typedef typename boost::mpl::end<Seq>::type last_;
   typedef typename exit_impl<first_,last_>::type exit_;
};

/*!
*
* It convert from a boost::fusion::vector into another boost::fusion::vector
* applying a meta-function F on each element
*
* usage:
*
* boost::fusion::vector<float,float,float[3]>::type is converted into
*
* S<F<float>::type,F<float>::type,F<float[3]>::type
*
* \param H metafunction (is a structure with typedef type)
* \param L boost::fusion::vector
*
*/

template<template<typename> class H,typename L>
struct to_variadic
{
	//! reverse the sequence
	typedef typename boost::mpl::reverse<L>::type reversed_;

	//! first element
	typedef typename seq_traits_impl<reversed_>::first_ first;

	//! last element
	typedef typename seq_traits_impl<reversed_>::last_ last;

	//! calculate the exit condition
	typedef typename exit_impl<first,last>::type exit_;

	//! generate the boost::fusion::vector apply H on each term
	typedef typename to_variadic_impl<H,first,last,exit_::value >::type type;
};

/*!
*
* It convert a variadic template into a boost::mpl::vector
*
* usage:
*
* to_boost_mpl<3,4,7,10>::type is converted into
*
* boost::mpl::vector<int<3>,int<4>,int<7>,int<10>>
*
*
*/

template <int a, int... id>
struct to_boost_mpl_impl
{
	//! push in front the next number
	typedef typename boost::mpl::push_front<typename to_boost_mpl_impl<id...>::type,boost::mpl::int_<a>>::type type;
};

//! terminator for to_boost_mpl with last parameter
template <int a>
struct to_boost_mpl_impl<a>
{
	typedef boost::mpl::vector<boost::mpl::int_<a>> type;
};

//! terminator for to_boost_mpl with last parameter
template <int... id>
struct to_boost_mpl
{
	typedef typename to_boost_mpl_impl<id...>::type type;
};

//! terminator for to_boost_mpl with last parameter
template <>
struct to_boost_mpl<>
{
	typedef typename boost::mpl::vector<>::type type;
};

#endif
