/*
 * memory_config.h
 *
 *  Created on: Aug 20, 2014
 *      Author: i-bird
 */

#include <boost/mpl/vector.hpp>

#ifndef MEMORY_CONFIG_H_
#define MEMORY_CONFIG_H_

#include "Point.hpp"
#include "memory_c.hpp"
#include "to_variadic.hpp"
#include <boost/fusion/include/end.hpp>
#include <boost/fusion/include/begin.hpp>
#include "t_to_memory_c.hpp"

template<typename T>
struct memory_cpu_type
{
	typedef typename T::type type;
	typedef typename T::type rtype;
};

#include "Point.hpp"

template<typename T>
class memory_gpu_type
{
	typedef T type;
};

#include <boost/type_traits/remove_all_extents.hpp>

template<typename T>
struct memory_gpu_type< Point<T> >
{
	typedef typename Point<T>::type ptype;

	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<0> >::type ptype_0;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<1> >::type ptype_1;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<2> >::type ptype_2;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<3> >::type ptype_3;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<4> >::type ptype_4;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<5> >::type ptype_5;

	typedef typename boost::remove_reference<ptype_0>::type mt_0;
	typedef typename boost::remove_reference<ptype_1>::type mt_1;
	typedef typename boost::remove_reference<ptype_2>::type mt_2;
	typedef typename boost::remove_reference<ptype_3>::type mt_3;
	typedef typename boost::remove_reference<ptype_4>::type mt_4;
	typedef typename boost::remove_reference<ptype_5>::type mt_5;

	typedef typename boost::remove_all_extents<mt_4>::type mt_4ele;
	typedef typename boost::remove_all_extents<mt_4>::type mt_5ele;

	typedef boost::fusion::vector< memory_c<mt_0>,memory_c<mt_1>,memory_c<mt_2>,memory_c<mt_3>,memory_c<multi_array<boost::mpl::vector<mt_4ele,boost::mpl::int_<3>>>>,memory_c<multi_array<boost::mpl::vector<mt_5ele,boost::mpl::int_<3>,boost::mpl::int_<3>>>> > type;
};




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
	typedef to_variadic<t_to_memory_c,typename boost::fusion::result_of::begin<Seq>::type ,typename boost::fusion::result_of::end<Seq>::type > type;
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
	//! we define the type structure as ptype, eapected to be a boost::mpl/fusion::vector
	typedef typename T::type ptype;

	//! for each element in the vector interleave memory_c
	typedef typename inter_memc<ptype>::type type;
};

#endif /* MEMORY_CONFIG_H_ */
