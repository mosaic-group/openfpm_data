/*
 * t_to_memory_c.hpp
 *
 *  Created on: Aug 27, 2014
 *      Author: Pietro Incardona
 */

#ifndef T_TO_MEMORY_C_HPP_
#define T_TO_MEMORY_C_HPP_

#include <boost/mpl/int.hpp>
#include <memory_c.hpp>

/*! \brief t_to_memory_c is a metafunction that given T it convert it into
 *
 * t_to_memory_c<T>::type to memory_c<T> for scalar
 *
 * or
 *
 * t_to_memory_c<T[N1][N2]>::type to memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int<N1>,boost::mpl::int<N2>>>>
 *
 * for vector
 *
 *
 * Here we specialize t_to_memory_c for several dimensionalities
 *
 * \param T type
 * \param N1 .... N10 dimensions
 *
 */

#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>

// First thing the metafunction has to do, is to remove eventually reference and const attribute

template<typename T>
struct remove_attributes_const_ref
{
	typedef typename boost::remove_const <typename boost::remove_reference< T >::type >::type type;
};

//! Partial specialization for scalar N=0
template<typename T>
struct t_to_memory_c_impl
{
	typedef memory_c<T> type;
};

//! Partial specialization for N=1
template<typename T,size_t N1>
struct t_to_memory_c_impl<T[N1]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>>>> type;
};

//! Partial specialization for N=2
template<typename T,size_t N1,size_t N2>
struct t_to_memory_c_impl<T[N1][N2]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>,
			                                          boost::mpl::int_<N2>>>> type;
};


//! Partial specialization for N=3
template<typename T,size_t N1,size_t N2,size_t N3>
struct t_to_memory_c_impl<T[N1][N2][N3]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>,
			                                          boost::mpl::int_<N2>,
			                                          boost::mpl::int_<N3>>>> type;
};

//! Partial specialization for N=4
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4>
struct t_to_memory_c_impl<T[N1][N2][N3][N4]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>,
			                                          boost::mpl::int_<N2>,
			                                          boost::mpl::int_<N3>,
			                                          boost::mpl::int_<N4>>>> type;
};

//! Partial specialization for N=5
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5>
struct t_to_memory_c_impl<T[N1][N2][N3][N4][N5]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>,
			                                          boost::mpl::int_<N2>,
			                                          boost::mpl::int_<N3>,
			                                          boost::mpl::int_<N4>,
			                                          boost::mpl::int_<N5>>>> type;
};

//! Partial specialization for N=6
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6>
struct t_to_memory_c_impl<T[N1][N2][N3][N4][N5][N6]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>,
			                                          boost::mpl::int_<N2>,
			                                          boost::mpl::int_<N3>,
			                                          boost::mpl::int_<N4>,
			                                          boost::mpl::int_<N5>,
			                                          boost::mpl::int_<N6>>>> type;
};

//! Partial specialization for N=7
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7>
struct t_to_memory_c_impl<T[N1][N2][N3][N4][N5][N6][N7]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>,
			                                          boost::mpl::int_<N2>,
			                                          boost::mpl::int_<N3>,
			                                          boost::mpl::int_<N4>,
			                                          boost::mpl::int_<N5>,
			                                          boost::mpl::int_<N6>,
			                                          boost::mpl::int_<N7>>>> type;
};

//! Partial specialization for N=8
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8>
struct t_to_memory_c_impl<T[N1][N2][N3][N4][N5][N6][N7][N8]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>,
			                                          boost::mpl::int_<N2>,
			                                          boost::mpl::int_<N3>,
			                                          boost::mpl::int_<N4>,
			                                          boost::mpl::int_<N5>,
			                                          boost::mpl::int_<N6>,
			                                          boost::mpl::int_<N7>,
			                                          boost::mpl::int_<N8>>>> type;
};

//! Partial specialization for N=9
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9>
struct t_to_memory_c_impl<T[N1][N2][N3][N4][N5][N6][N7][N8][N9]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>,
			                                          boost::mpl::int_<N2>,
			                                          boost::mpl::int_<N3>,
			                                          boost::mpl::int_<N4>,
			                                          boost::mpl::int_<N5>,
			                                          boost::mpl::int_<N6>,
			                                          boost::mpl::int_<N7>,
			                                          boost::mpl::int_<N8>,
			                                          boost::mpl::int_<N9>>>> type;
};

//! Partial specialization for N=10
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9, size_t N10>
struct t_to_memory_c_impl<T[N1][N2][N3][N4][N5][N6][N7][N8][N9][N10]>
{
	typedef memory_c<multi_array<boost::mpl::vector<T,boost::mpl::int_<N1>,
			                                          boost::mpl::int_<N2>,
			                                          boost::mpl::int_<N3>,
			                                          boost::mpl::int_<N4>,
			                                          boost::mpl::int_<N5>,
			                                          boost::mpl::int_<N6>,
			                                          boost::mpl::int_<N7>,
			                                          boost::mpl::int_<N8>,
			                                          boost::mpl::int_<N9>,
			                                          boost::mpl::int_<N10>>>> type;
};

/*! \brief Meta-function t_to_memory_c
 *
 * Meta-function t_to_memory_c, convert the type T into an appropriate memory_c type
 *
 * basically it convert t_to_memory_c<float> into memory_c<float> and array specification like
 * t_to_memory_c<float[3][3]> into memory_c<multy_array<float,3,3>>, and so on for higher dimensionality
 *
 */
template<typename T>
struct t_to_memory_c
{
	typedef typename t_to_memory_c_impl<typename remove_attributes_const_ref<T>::type>::type type;
};


#endif /* T_TO_MEMORY_C_HPP_ */
