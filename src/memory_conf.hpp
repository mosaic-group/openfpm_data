/*
 * memory_conf.hpp
 *
 *  Created on: Aug 28, 2014
 *      Author: Pietro Incardona
 */

#ifndef MEMORY_CONF_HPP_
#define MEMORY_CONF_HPP_

#include "to_variadic.hpp"
#include "t_to_memory_c.hpp"

/*! \brief This class convert a boost::mpl::fusion/vector to a boost::mpl::fusion/vector with memory_c interleaved
 *
 * This class convert a boost::mpl::fusion/vector to a boost::mpl::fusion/vector with memory_c interleaved
 *
 * Example:
 *
 * typedef boost::mlp::vector<float,float,float[3][3], ... > A
 *
 * inter_memc<A>
 *
 * is equivalent to
 *
 * boost::fusion::vector<memory_c<float>,memory_c<float>, memory_c<multi_array<boost::mpl::vector<float,3,3>, ...... >
 *
 * \param Seq Is suppose to be an boost::mpl::vector/fusion
 *
 */
template<typename Seq>
struct inter_memc
{
	typedef typename to_variadic<t_to_memory_c,Seq,typename boost::mpl::begin<Seq>::type ,typename boost::mpl::end<Seq>::type >::type type;
};

/*! \brief Transform the boost::fusion::vector into memory specification (memory_traits)
 *
 * Transform the boost::fusion::vector into memory specification (memory_traits).
 * In this implementation we interleave each property of the base type with memory_c
 *
 * We basically create a buffer for each property
 *
 * \param T base type (T::type must define a boost::mpl::vector )
 *
 *
 */

template<typename T>
struct memory_traits_inte
{
	//! for each element in the vector interleave memory_c
	typedef typename inter_memc<T>::type type;
};

/*! \brief Transform the boost::fusion::vector into memory specification (memory_traits)
 *
 * Transform the boost::fusion::vector into memory specification (memory_traits).
 * In this implementation we create a buffer of base type with memory_c
 *
 * We basically create a buffer for each property
 *
 * \param T base type (T::type must define a boost::mpl::vector )
 *
 *
 */

template<typename T>
struct memory_traits_lin
{
	//! for each element in the vector interleave memory_c
	typedef memory_c<boost::fusion::vector<float, float, float, float, float [3], float [3][3]>> type;
};

#endif
