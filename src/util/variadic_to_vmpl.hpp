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

#ifndef V_TRANSFORM_HPP
#define V_TRANSFORM_HPP

#include <boost/fusion/container/vector.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/reverse.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>

////////////////////////////////////////////////////////////////////////////////////////////////////

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

/*! \brief Recursive specialization of v_transform
 *
 * \param F suppose to be the original end of boost::mpl::vector
 * \param L suppose to be the actual position of the boost::mpl::vector
 * \param exit, when true it say to terminate the sequence
 *
 */
template<template<typename> class H, typename F,typename L, bool exit,typename ...Args>
struct v_transform_impl
{
   typedef typename boost::mpl::deref<F>::type front_;
   typedef typename boost::mpl::next<F>::type next_;
   typedef typename exit_impl<next_,L>::type exit_;
   typedef typename v_transform_impl<H,next_,L,exit_::value,typename H<front_>::type,Args...>::type type;
};


//! Terminator of to_variadic
template<template<typename> class H,typename F,typename L,typename ...Args>
struct v_transform_impl<H,F,L,true,Args...>
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
* It transform a boost::fusion::vector to another boost::fusion::vector
* applying a meta-function H on each element
*
* \param H metafunction (is a structure with typedef type)
* \param L boost::fusion::vector
*
* ### Meta-function definition
* \snippet variadic_to_vmpl_unit_test.hpp v_transform metafunction
*
* ### Usage
* \snippet variadic_to_vmpl_unit_test.hpp v_transform usage
*
*/

template<template<typename> class H,typename L>
struct v_transform
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
	typedef typename v_transform_impl<H,first,last,exit_::value >::type type;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*! \brief Recursive specialization of v_transform in case of metafunction with 2 argument
 *
 * \param H the metafunction
 * \param F suppose to be the original end of boost::mpl::vector
 * \param L suppose to be the actual position of the boost::mpl::vector
 * \param exit, when true it say to terminate the sequence
 *
 */
template<template<typename,typename> class H, typename arg0, typename F,typename L, bool exit,typename ...Args>
struct v_transform_two_impl
{
   typedef typename boost::mpl::deref<F>::type front_;
   typedef typename boost::mpl::next<F>::type next_;
   typedef typename exit_impl<next_,L>::type exit_;
   typedef typename v_transform_two_impl<H,arg0,next_,L,exit_::value,typename H<arg0,front_>::type,Args...>::type type;
};


//! Terminator of to_variadic
template<template<typename,typename> class H,typename arg0, typename F,typename L,typename ...Args>
struct v_transform_two_impl<H,arg0,F,L,true,Args...>
{
   typedef boost::fusion::vector<Args...> type;
};

/*!
*
* It transform a boost::fusion::vector to another boost::fusion::vector
* applying a meta-function H on each element
*
* \param H 2-argument metafunction (is a structure with typedef type)
* \param L boost::fusion::vector
*
* ### Meta-function definition
* \snippet variadic_to_vmpl_unit_test.hpp v_transform_two metafunction
*
* ### Usage
* \snippet variadic_to_vmpl_unit_test.hpp v_transform_two usage
*
*/

template<template<typename,typename> class H,typename arg0, typename L>
struct v_transform_two
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
	typedef typename v_transform_two_impl<H,arg0,first,last,exit_::value >::type type;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <int a, int... id>
struct to_boost_vmpl_impl
{
	//! push in front the next number
	typedef typename boost::mpl::push_front<typename to_boost_vmpl_impl<id...>::type,boost::mpl::int_<a>>::type type;
};

//! terminator for to_boost_mpl with last parameter
template <int a>
struct to_boost_vmpl_impl<a>
{
	typedef boost::mpl::vector<boost::mpl::int_<a>> type;
};

/*!
*
* It convert a variadic template into a boost::mpl::vector
*
* usage:
*
* to_boost_mpl<3,4,7,10>::type is converted into
*
* boost::mpl::vector<int_<3>,int_<4>,int_<7>,int_<10>>
*
*
*/
template <int... id>
struct to_boost_vmpl
{
	typedef typename to_boost_vmpl_impl<id...>::type type;
};

//! terminator for to_boost_mpl with last parameter
template <>
struct to_boost_vmpl<>
{
	typedef typename boost::mpl::vector<>::type type;
};

#endif
