/*
 * tokernel_transformation.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: i-bird
 */

#ifndef TOKERNEL_TRANSFORMATION_HPP_
#define TOKERNEL_TRANSFORMATION_HPP_

#include "data_type/aggregate.hpp"

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

template<template <typename> class layout_base, typename T, bool = is_vector_native<T>::value>
struct toKernel_transform;

template<template <typename> class layout_base, typename ... args>
struct apply_trasform_impl
{
	typedef void type;
};

template<template <typename> class layout_base, typename ... args>
struct apply_trasform_impl<layout_base,boost::fusion::vector<args...>>
{
	typedef aggregate<typename toKernel_transform<layout_base,args>::type ... > type;
};

template<template <typename> class layout_base,typename T>
struct apply_transform
{
	typedef typename apply_trasform_impl<layout_base,typename T::type>::type type;
};

template<template <typename> class layout_base, typename T, bool >
struct toKernel_transform
{
	typedef T type;
};


template<template <typename> class layout_base, typename T>
struct toKernel_transform<layout_base,T,true>
{
	typedef typename apply_transform<layout_base,typename T::value_type>::type aggr;

	typedef openfpm::vector_gpu_ker<aggr,layout_base> type;
};

#endif /* TOKERNEL_TRANSFORMATION_HPP_ */
