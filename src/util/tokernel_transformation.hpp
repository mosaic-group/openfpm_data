/*
 * tokernel_transformation.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: i-bird
 */

#ifndef TOKERNEL_TRANSFORMATION_HPP_
#define TOKERNEL_TRANSFORMATION_HPP_

#include "data_type/aggregate.hpp"

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to call hostToDevice for each properties
 *
 */
template<typename Tv>
struct host_to_dev_all_prp
{
	Tv & p;

	inline host_to_dev_all_prp(Tv & p)
	:p(p)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		p.template hostToDevice<T::value>();
	}
};

template<typename T, typename T_ker, typename type_prp, template<typename> class layout_base , int is_vector>
struct call_recursive_host_device_if_vector
{
	template<typename mem_type, typename obj_type> static void transform(mem_type * mem, obj_type & obj, size_t start, size_t stop)
	{
		start /= sizeof(type_prp);
		stop /= sizeof(type_prp);

		// The type of device and the type on host does not match (in general)
		// So we have to convert before transfer

		T * ptr = static_cast<T *>(obj.get_pointer());

		mem_type tmp;

		tmp.allocate(mem->size());

		T_ker * ptr_tt = static_cast<T_ker *>(tmp.getPointer());

		for(size_t i = start ; i < stop ; i++)
		{
			new (&ptr_tt[i]) T_ker();
			ptr_tt[i] = ptr[i].toKernel();
		}

		mem->hostToDevice(tmp);
	}

	//! It is a vector recursively call deviceToHost
	template<typename obj_type>
	static void call(obj_type & obj, size_t start, size_t stop)
	{
		T * ptr = static_cast<T *>(obj.get_pointer());

		for(size_t i = start ; i < stop ; i++)
		{
			host_to_dev_all_prp<T> hdap(ptr[i]);

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::value_type::max_prop>>(hdap);
		}
	}
};

template<typename T, typename T_ker, typename type_prp ,template<typename> class layout_base>
struct call_recursive_host_device_if_vector<T,T_ker,type_prp,layout_base,0>
{
	template<typename mem_type,typename obj_type> static void transform(mem_type * mem, obj_type & obj, size_t start, size_t stop)
	{
		mem->hostToDevice(start,stop);
	}

	//! It is not a vector nothing to do
	template<typename obj_type>
	static void call(obj_type & obj, size_t start, size_t stop) {}
};

template<typename T, typename T_ker, typename type_prp ,template<typename> class layout_base>
struct call_recursive_host_device_if_vector<T,T_ker,type_prp,layout_base,3>
{
	template<typename mem_type,typename obj_type> static void transform(mem_type * mem, obj_type & obj, size_t start, size_t stop)
	{
		// calculate the start and stop elements
		start /= std::extent<type_prp,0>::value;
		stop /= std::extent<type_prp,0>::value;
		size_t sz = mem->size() / std::extent<type_prp,0>::value;

		size_t offset = 0;
		for (size_t i = 0 ; i < std::extent<type_prp,0>::value ; i++)
		{
			mem->hostToDevice(offset+start,offset+stop);
			offset += sz;
		}
	}

	//! It is not a vector nothing to do
	template<typename obj_type>
	static void call(obj_type & obj, size_t start, size_t stop) {}
};

template<typename T, typename T_ker, typename type_prp ,template<typename> class layout_base>
struct call_recursive_host_device_if_vector<T,T_ker,type_prp,layout_base,4>
{
	template<typename mem_type,typename obj_type> static void transform(mem_type * mem, obj_type & obj, size_t start, size_t stop)
	{
		// calculate the start and stop elements
		start = start / std::extent<type_prp,0>::value / std::extent<type_prp,1>::value;
		stop = stop / std::extent<type_prp,0>::value / std::extent<type_prp,1>::value;
		size_t sz = mem->size() / std::extent<type_prp,0>::value / std::extent<type_prp,1>::value;

		size_t offset = 0;
		for (size_t i = 0 ; i < std::extent<type_prp,0>::value ; i++)
		{
			for (size_t j = 0 ; j < std::extent<type_prp,1>::value ; j++)
			{
				mem->hostToDevice(offset+start,offset+stop);
				offset += sz;
			}
		}
	}

	//! It is not a vector nothing to do
	template<typename obj_type>
	static void call(obj_type & obj, size_t start, size_t stop) {}
};

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

template<template <typename> class layout_base, typename T, int = is_vector_native<T>::value + 2*is_vector_dist<T>::value + 4*is_gpu_celllist<T>::value >
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

template<unsigned int dim ,typename T> class Point;

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
