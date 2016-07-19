/*
 * memory_conf.hpp
 *
 *  Created on: Aug 28, 2014
 *      Author: Pietro Incardona
 */

#ifndef MEMORY_CONF_HPP_
#define MEMORY_CONF_HPP_

#include "util/variadic_to_vmpl.hpp"
#include "t_to_memory_c.hpp"
#include "Vector/vect_isel.hpp"

/*! \brief This class convert a boost::mpl::fusion/vector to a boost::mpl::fusion/vector with memory_c interleaved
 *
 * This class convert a boost::mpl::fusion/vector to a boost::mpl::fusion/vector with memory_c interleaved
 *
 * Example:
 *
 * typedef boost::mpl::vector<float,float,float[3][3], ... > A
 *
 * inter_memc<A>
 *
 * produce
 *
 * boost::fusion::vector<memory_c<float>,memory_c<float>, memory_c<multi_array<boost::mpl::vector<float,3,3>, ...... >
 *
 * \param Seq Is suppose to be an boost::mpl::vector/fusion
 *
 */
template<typename Seq>
struct inter_memc
{
	typedef typename v_transform<t_to_memory_c,Seq>::type type;
};

/*! \brief Transform the boost::fusion::vector into memory specification (memory_traits)
 *
 * Transform the boost::fusion::vector into memory_traits.
 * In this implementation we interleave each property of the base type with memory_c
 *
 * We basically create a buffer for each property
 *
 * \see see inter_mem_c for detail
 *
 * \param T base type (T::type must define a boost::fusion::vector )
 *
 *
 */
template<typename T>
struct memory_traits_inte
{
	//! for each element in the vector interleave memory_c
	typedef typename inter_memc<typename T::type>::type type;

	typedef int yes_is_inte;
};

/*! \brief small meta-function to get the type of the memory
 *
 *
 */
template<typename T, bool is_agg>
struct memory_traits_lin_type
{
	typedef memory_c<typename T::type> type;
};

/*! \brief small meta-function to get the type of the memory
 *
 *
 */
template<typename T>
struct memory_traits_lin_type<T,false>
{
	typedef void type;
};

/*! \brief Transform the boost::fusion::vector into memory specification (memory_traits)
 *
 * Transform the boost::fusion::vector into memory specification (memory_traits).
 * In this implementation we create a buffer of base type with memory_c
 *
 * We basically create a buffer for each property
 *
 * \param T base type (T::type must define a boost::fusion::vector )
 *
 *
 */

template<typename T>
struct memory_traits_lin
{
	//! for each element in the vector interleave memory_c
	typedef typename memory_traits_lin_type<T,openfpm::vect_isel<T>::value == OPENFPM_NATIVE>::type type;

	typedef int yes_is_tlin;
};


//////////////////////////////////////////////////////////////

template<typename T, typename Sfinae = void>
struct is_layout_mlin: std::false_type {};


/*! \brief is_layout_mlin
 *
 * ### Example
 *
 * \snippet util.hpp Check if the memory layout is memory_traits_lin
 *
 * return true if T is a memory_traits_lin
 *
 */
template<typename T>
struct is_layout_mlin<T, typename Void< typename T::yes_is_tlin>::type> : std::true_type
{};


template<typename T, typename Sfinae = void>
struct is_layout_inte: std::false_type {};


/*! \brief is_layout_inte
 *
 * ### Example
 *
 * \snippet util.hpp Check if the memory layout is memory_traits_inte
 *
 * return true if T is a memory_traits_inte
 *
 */
template<typename T>
struct is_layout_inte<T, typename Void< typename T::yes_is_inte>::type> : std::true_type
{};

#endif
