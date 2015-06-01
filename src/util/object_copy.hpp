/*
 * vector_prop_copy.hpp
 *
 *  Created on: May 31, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_PROP_COPY_HPP_
#define VECTOR_PROP_COPY_HPP_

#include "for_each_ref.hpp"
#include <boost/mpl/range_c.hpp>
#include <boost/fusion/include/size.hpp>

/*! \brief This is a container for the sending buffers
 *
 * It is used in ghost_get to create a particular object with the properties selected
 *
 * \tparam Is a boost::fusion::vector with the properties selected
 *
 *
 */
template<typename v>
struct object
{
	typedef v type;
	typedef typename boost::fusion::result_of::size<v>::type size_tpy;

	type data;

	static const int max_prop = size_tpy::value;
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy the selected properties
 *
 * \tparam v_src source objecy
 * \tparam d_src destination object
 *
 */

template<typename v_src,typename v_dst, int... prp>
struct object_copy_f
{
	// Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	// Source object
	v_src & src;

	// Destination object
	v_dst & dst;

	/*! \brief Constructor
	 *
	 * Create a vertex properties list
	 *
	 * \param v_node std::string that is filled with the graph properties in the GraphML format
	 * \param n_obj object container to access its properties for example encapc<...>
	 *
	 */
	object_copy_f(v_src & src, v_dst & dst)
	:src(src),dst(dst)
	{
	};

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
    	typedef typename boost::mpl::at<typename v_dst::type,typename boost::mpl::int_<T::value>>::type ctype;

    	meta_copy<ctype>(boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(src.data),boost::fusion::at_c<T::value>(dst.data));
    }
};

/*! \brief It copy the properties from one object to another
 *
 * Given a set of properties for the source (0,1,3) it copy that property to the destination properties
 * (0,1,2)
 *
 * [Example]
 *
 * boost::fusion::vector<float,double,int,float[3]> src;
 *
 * boost::fusion::vector<float,float[3]> dst;
 *
 * object_copy<0,1,3>(src,dst)
 *
 */
template<typename v_src, typename v_dst, int... prp>
inline void object_copy(v_src & vs, v_dst & vd)
{
	object_copy_f<v_src,v_dst,prp...> obj(vs,vd);
	boost::mpl::for_each_ref< boost::mpl::range_c<int,0,v_dst::max_prop> >(obj);
};

#endif /* VECTOR_PROP_COPY_HPP_ */
