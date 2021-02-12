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
* index_tuple_sq<3,5,6>
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


//////////////////////////////////////////////////////////////////////////

/*! \brief Recursive specialization of to_variadic
 *
 * \param H suppose to be the vector where we are pushing numbers
 * \param c suppose to be the counter
 * \param L suppose to be the end of the cycle
 * \param exit, when true it say to terminate the sequence
 *
 */
template<int c , int end, int ele, bool exit,int ... vars>
struct to_variadic_const_impl
{
   typedef typename exit_impl_sq<c,end-1>::type exit_;
   typedef typename to_variadic_const_impl<c+1,end,ele,exit_::value,vars ..., ele>::type type;
};


//! Terminator of to_variadic
template<int c, int end, int ele, int ... vars>
struct to_variadic_const_impl<c,end,ele,true, vars ...>
{
   typedef boost::mpl::vector<boost::mpl::int_<vars> ...> type;
};

/*!
*
* It create a mpl::vector sequence of integer with a constant A
*
* usage:
*
* vmpl_create_constant<3,5>::type
*
* is mapped to
*
* boost::mpl::vector<int_<5>,int_<5>,int<5>>
*
*/
template<unsigned int N, unsigned int M>
struct vmpl_create_constant
{
	//! end condition
	typedef typename exit_impl_sq<N,M>::type exit_;

	//! generate the boost::fusion::vector apply H on each term
	typedef typename to_variadic_const_impl<1,N,M,exit_::value,M>::type type;
};

//////////////////

template<typename bint, unsigned int ele>
struct sum_ele
{
};

template<int op1, unsigned int ele>
struct sum_ele<boost::mpl::int_<op1>,ele>
{
	typedef boost::mpl::int_<op1+ele> type;
};

template<unsigned int ele>
struct sum_ele<boost::mpl::na,ele>
{
	typedef boost::mpl::na type;
};


template <unsigned int ele, typename ... vmpl>
struct vmpl_sum_constant_impl
{
	//! construct an mpl vector from the variadic
	typedef boost::mpl::vector<typename sum_ele<vmpl,ele>::type ...> type;
};

template <unsigned int ele, typename vmpl>
struct vmpl_sum_constant
{
	typedef int type;
};
/*!
*
* Take a boost::mpl::vector of boost::mpl::int_ and sum the element ele
*
* vmpl_sum_constant<boost::mpl::vector<boost::mpl::int_<2>,boost::mpl::int_<4>,boost::mpl::int<8>>, 5>::type is converted into
*
* boost::mpl::vector<int_<7>,int_<9>,int_<13>>
*
* \snippet variadic_to_vmpl_unit_test.hpp vmpl_sum_constant usage
*
*/
template <unsigned int ele, typename ... vars>
struct vmpl_sum_constant<ele, boost::mpl::vector<vars ... > >
{
	//! construct an mpl vector from the variadic
	typedef boost::mpl::vector<typename sum_ele<vars,ele>::type ...> type;
};



////////////////////////////////////////////////////////////////////////

template<int c, int accu, int stop, typename vmpl, bool exit>
struct vmpl_reduce_prod_impl
{
   typedef typename exit_impl_sq<c,stop>::type exit_;
   typedef typename boost::mpl::at<vmpl,boost::mpl::int_<c>>::type ele;
   typedef typename vmpl_reduce_prod_impl<c+1,accu*ele::value,stop,vmpl,exit_::value>::type type;
};


//! Terminator of to_variadic
template<int c, int accu, int stop, typename vmpl>
struct vmpl_reduce_prod_impl<c,accu,stop,vmpl,true>
{
   typedef boost::mpl::int_<accu> type;
};


template<typename vmpl>
struct vmpl_reduce_prod
{
	typedef typename vmpl_reduce_prod_impl<0,1,(int)vmpl::size::value-1,vmpl,false>::type type;
};

template<typename vmpl, int stop>
struct vmpl_reduce_prod_stop
{
	typedef typename vmpl_reduce_prod_impl<0,1,stop,vmpl,false>::type type;
};

template<typename vmpl>
struct vmpl_reduce_prod_stop<vmpl,-1>
{
	typedef typename boost::mpl::int_<1> type;
};

//! Linearize a set of index
template<typename vmpl, typename a> __device__ __host__ inline unsigned int Lin_vmpl(a v)
{
	return v*vmpl_reduce_prod_stop<vmpl,((int)vmpl::size::value) - 2>::type::value;
}

/*! \brief linearize an arbitrary set of index
 *
 * linearize an arbitrary set of index
 *
 */
template<typename vmpl, typename a, typename ...lT>
__device__ __host__ inline unsigned int Lin_vmpl(a v,lT...t)
{
		return v*vmpl_reduce_prod_stop<vmpl,(int)vmpl::size::value - (int)sizeof...(t) - 2>::type::value + Lin_vmpl<vmpl>(t...);
}

//! Linearize a set of index
template<typename vmpl, typename vmpl_off, typename a> __device__ __host__ inline unsigned int Lin_vmpl_off(a v)
{
	return (v + boost::mpl::at<vmpl_off,boost::mpl::int_< ((int)vmpl::size::value) - 1 > >::type::value )*vmpl_reduce_prod_stop<vmpl,(int)vmpl::size::value - 2>::type::value;
}

/*! \brief linearize an arbitrary set of index
 *
 * linearize an arbitrary set of index
 *
 */
template<typename vmpl, typename vmpl_off, typename a, typename ...lT>
__device__ __host__ inline unsigned int Lin_vmpl_off(a v,lT...t)
{
	return (v + boost::mpl::at<vmpl_off,boost::mpl::int_<((int)vmpl::size::value) - (int)sizeof...(t) - 1> >::type::value)*vmpl_reduce_prod_stop<vmpl,(int)((int)vmpl::size::value - sizeof...(t)  - 2)>::type::value + Lin_vmpl_off<vmpl,vmpl_off>(t...);
}

#endif
