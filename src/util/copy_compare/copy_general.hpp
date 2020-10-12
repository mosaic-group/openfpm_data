/*
 * aggregate_copy.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_COPY_GENERAL_HPP_
#define OPENFPM_DATA_SRC_UTIL_COPY_GENERAL_HPP_

#include "util/common.hpp"
#include "util/util_debug.hpp"
#include "copy_compare_aggregates.hpp"
#include "util/for_each_ref.hpp"
#include <boost/mpl/range_c.hpp>
#include <iostream>
#include "util/cuda_util.hpp"
#include "data_type/aggregate.hpp"

/*! \brief This structure define the operation add to use with copy general
 *
 * \tparam Tdst destination object type
 * \tparam Tsrc source object type
 *
 */
template<typename Tdst, typename Tsrc>
struct max_
{
	/*! \brief Defition of the add operation
	 *
	 * \param dst Destination object
	 * \param src Source object
	 *
	 */
	static inline void operation(Tdst & dst, const Tsrc & src)
	{
		dst = (src > dst)?src:dst;
	}
};

/*! \brief This structure define the operation add to use with copy general
 *
 * \tparam Tdst destination object type
 * \tparam Tsrc source object type
 *
 */
template<typename Tdst, typename Tsrc>
struct min_
{
	/*! \brief Defition of the add operation
	 *
	 * \param dst Destination object
	 * \param src Source object
	 *
	 */
	static inline void operation(Tdst & dst, const Tsrc & src)
	{
		dst = (src < dst)?src:dst;
	}
};

/*! \brief This structure define the operation add to use with copy general
 *
 * \tparam Tdst destination object type
 * \tparam Tsrc source object type
 *
 */
template<typename Tdst, typename Tsrc>
struct add_
{
	/*! \brief Defition of the add operation
	 *
	 * \param dst Destination object
	 * \param src Source object
	 *
	 */
	__device__ __host__ static inline void operation(Tdst & dst, const Tsrc & src)
	{
		dst += src;
	}
};

/*! \brief This structure define the operation add to use with copy general
 *
 * \tparam Tdst destination object type
 * \tparam Tsrc source object type
 *
 */
template<typename Tdst, typename Tsrc>
struct add_atomic_
{
	/*! \brief Defition of the add operation
	 *
	 * \param dst Destination object
	 * \param src Source object
	 *
	 */
	__device__ __host__ static inline void operation(Tdst & dst, const Tsrc & src)
	{
		atomicAdd(&dst,src);
	}
};

/*! \brief This structure define the operation add to use with copy general
 *
 * \tparam Tdst destination object type
 * \tparam Tsrc source object type
 *
 */
template<typename Tdst, typename Tsrc>
struct replace_
{
	/*! \brief Defition of the replace operation
	 *
	 * \param dst Destination object
	 * \param src Source object
	 *
	 */
	static inline void operation(Tdst & dst, const Tsrc & src)
	{
		dst = src;
	}
};

/*! \brief This structure define the operation add to use with copy general
 *
 * \tparam Tdst destination object type
 * \tparam Tsrc source object type
 *
 */
template<typename Tdst, typename Tsrc>
struct merge_
{
	/*! \brief Defition of the add operation
	 *
	 * \param dst Destination object
	 * \param src Source object
	 *
	 */
	static inline void operation(Tdst & dst, const Tsrc & src)
	{
		dst.add(src);
	}
};

/*! \brief structure to copy aggregates
 *
 * \tparam T type to copy
 *
 */
template<typename T, unsigned int agg=2 * is_aggregate<T>::value + std::is_copy_assignable<T>::value>
struct copy_general
{
	/*! \brief Specialization when there is unknown copy method
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	inline copy_general(const T & src, T & dst)
	{
#ifndef DISABLE_ALL_RTTI
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << "  " << demangle(typeid(T).name()) << " does not have an operator= and is not an aggregate or an openfpm native structure, copy is not possible" << "\n";
#endif
	}
};

//! Specialization for if dst type is copy assignable from src type
template<typename T>
struct copy_general<T,1>
{
	/*! \brief copy objects that has an operator= (implicit or explicit)
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	__device__ __host__ inline copy_general(const T & src, T & dst)
	{
		dst = src;
	}
};

//! Specialization for aggregate type object
template<typename T>
struct copy_general<T,2>
{
	/*! \brief copy objects that are aggregates
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	inline copy_general(const T & src, T & dst)
	{
		copy_aggregate<T> cp(src,dst);

		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::max_prop>>(cp);
	}
};

//! Specialization for aggregate type object that define an operator=
template<typename T>
struct copy_general<T,3>
{
	/*! \brief copy objects that are aggregates but define an operator=
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	inline copy_general(const T & src, T & dst)
	{
		dst = src;
	}
};

/////////////////// VERSION WITH OPERATIONS ///////////////////

/*! \brief structure to copy aggregates applying an operation
 *
 * \tparam T type to copy
 *
 */
template<template<typename,typename> class op, typename T, unsigned int agg=2 * is_aggregate<T>::value + std::is_copy_assignable<T>::value>
struct copy_general_op
{
	/*! \brief Specialization when there is unknown copy method
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	inline copy_general_op(const T & src, T & dst)
	{
#ifndef DISABLE_ALL_RTTI
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << "  " << demangle(typeid(T).name()) << " does not have an operator " << demangle(typeid(op<T,T>).name()) << "defined" << std::endl;
#endif
	}
};

//! Specialization for object that can be assigned with an operator copy
template<template<typename,typename> class op,typename T>
struct copy_general_op<op,T,1>
{
	/*! \brief copy objects that has an operator= (implicit or explicit)
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	inline copy_general_op(const T & src, T & dst)
	{
		op<T,T>::operation(dst,src);
	}
};

//! Specialization for aggregate type objects
template<template<typename,typename> class op, typename T>
struct copy_general_op<op,T,3>
{
	/*! \brief copy objects that are aggregates
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	inline copy_general_op(const T & src, T & dst)
	{
		copy_aggregate_op<op,T> cp(src,dst);

		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::max_prop>>(cp);
	}
};


#endif /* OPENFPM_DATA_SRC_UTIL_COPY_GENERAL_HPP_ */
