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

/*! \brief Exit condition
 *
 * Exit condition, when the distance between F and L is 0 return true
 *
 * \param F suppose to be the end of the mpl::vector
 * \param L suppose to be the actual position of the mpl::vector
 *
 */

template<typename F,typename L>
struct exit : boost::mpl::equal_to<typename boost::mpl::distance<F,L>::type,boost::mpl::int_<0>>
{};

/*! \brief Recursive specialization of to_variadic
 *
 * \param F suppose to be the original end of boost::mpl::vector
 * \param L suppose to be the actual position of the boost::mpl::vector
 * \param exit, when true it say to terminate the sequence
 *
 */
template<typename S,typename H, typename F,typename L, bool exit,typename ...Args>
struct to_variadic_impl
{
   typedef typename boost::mpl::deref<F>::type front_;
   typedef typename boost::mpl::next<F>::type next_;
   typedef typename exit<next_,L>::type exit_;
   typedef typename to_variadic<S,next_,L,exit_::value,H<front_>::type,Args...>::type type;
};


//! Terminator of to_variadic
template<typename S,typename F,typename L,typename ...Args>
struct to_variadic_impl<F,L,true,Args...>
{
   typedef S<Args...> type;
};

/*!
*
* It convert from a boost::mpl::vector into an S variadic template
* class appling a metafunction F on each element
*
* usage:
*
* boost::mpl::vector<float,float,float[3]>::type is converted into
*
* S<F<float>::type,F<float>::type,F<float[3]>::type
*
* \param S class
* \param H metafunction (is a structure with typedef type)
* \param L boost::mpl::vector
*
*/

template<typename S,typename H, typename F,typename L>
struct to_variadic
{
	//! reverse the sequence
	typedef typename boost::mpl::reverse<L>::type reversed_;

	typedef to_variadic_impl<S,H,typename seq_traits<reversed_>::first_,
		      typename seq_traits<reversed_>::last_,
		      typename seq_traits<L>::exit_::value> type;
};

template<typename Seq>
struct seq_traits{
   typedef typename boost::mpl::begin<Seq>::type first_;
   typedef typename boost::mpl::end<Seq>::type last_;
   typedef typename exit<first_,last_>::type exit_;
};


#endif
