/*
 * tokernel_transformation.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: i-bird
 */

#ifndef TOKERNEL_TRANSFORMATION_HPP_
#define TOKERNEL_TRANSFORMATION_HPP_

#include "data_type/aggregate.hpp"

/*! \brief this set of meta-functions traverse at compile time the tree-structure of types in Depth-first search.
 *         and transform any root node of type vector into vector_gpu_ker
 *
 * Consider
 *
 * vector_gpu<aggregate<int,vector_gpu<aggregate<int,float>>>>
 *
 * is a tree in this form
 *
 * \verbatim
 *
 *                           * vector_gpu<aggregate<...>>
 *                          / \
 *                         /   \
 *                        /     \
 *                      int      * vector_gpu<aggregate<...>>
 *                              / \
 *                             /   \
 *                            /     \
 *                          int    float
 *
 * \endverbatim
 *
 * The vector is transformed at compile-time into
 *
 * is a tree in this form
 *
 * \verbatim
 *
 *                           * vector_gpu_ker<aggregate<...>>
 *                          / \
 *                         /   \
 *                        /     \
 *                      int      * vector_gpu_ker<aggregate<...>>
 *                              / \
 *                             /   \
 *                            /     \
 *                          int    float
 *
 * \endverbatim
 *
 *
 */

namespace openfpm
{

	/*! \brief grid interface available when on gpu
	 *
	 * \tparam n_buf number of template buffers
	 *
	 */

	template<typename T, template <typename> class layout_base>
	struct vector_gpu_ker;
}

// Definition of the box
template<unsigned int dim , typename T> class Box;

template<template <typename> class layout_base, typename T, int = is_vector_native<T>::value + 2*is_vector_dist<T>::value >
struct toKernel_transform;

template<template <typename> class layout_base, typename T, typename ... args>
struct apply_trasform_impl
{
	typedef void type;
};

template<template <typename> class layout_base, typename T, int impl, typename ... args>
struct aggregate_or_known_type
{
	typedef aggregate<typename toKernel_transform<layout_base,args>::type ... > type;
};

template<template <typename> class layout_base, typename T, typename ... args>
struct apply_trasform_impl<layout_base,T,boost::fusion::vector<args...>>
{
	static const int impl = is_aggregate<T>::value + is_Box<T>::value * 2 + is_Point<T>::value * 4;

	typedef typename aggregate_or_known_type<layout_base,T,impl,args...>::type type;
};


template<template <typename> class layout_base,typename T>
struct apply_transform
{
	typedef typename apply_trasform_impl<layout_base,T,typename T::type>::type type;
};

/////////////////////////////////////////////// TRANSFORMER NODE /////////////////////////////////////////////////

template<template <typename> class layout_base, typename T >
struct toKernel_transform<layout_base,T,0>
{
	typedef T type;
};


template<template <typename> class layout_base, typename T>
struct toKernel_transform<layout_base,T,1>
{
	typedef typename apply_transform<layout_base,typename T::value_type>::type aggr;

	typedef openfpm::vector_gpu_ker<aggr,layout_base> type;
};

/////////////////////////////////////////////////// KNOWN TYPE SPECIALIZATION TERMINATORS //////////////////////

template<template <typename> class layout_base,typename T, typename ... args>
struct aggregate_or_known_type<layout_base,T,2,args ...>
{
	typedef Box<T::dims,typename T::btype > type;
};

template<template <typename> class layout_base,typename T, typename ... args>
struct aggregate_or_known_type<layout_base,T,4,args ...>
{
	typedef Point<T::dims,typename T::coord_type > type;
};

#endif /* TOKERNEL_TRANSFORMATION_HPP_ */
