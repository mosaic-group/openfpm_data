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

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy the selected properties
 *
 * \tparam v_src source object
 * \tparam d_src destination object
 *
 */

template<typename v_src,typename v_dst, int... prp>
struct object_si_d_e
{
	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	//! Source object
	const v_src & src;

	//! Destination object
	v_dst & dst;

	/*! \brief Constructor
	 *
	 * \param src source object
	 * \param dst destination object
	 *
	 */
	__device__ __host__ object_si_d_e(const v_src & src, v_dst & dst)
	:src(src),dst(dst)
	{
	};

	//! It call the functor for each member
    template<typename T>
    __device__ __host__ void operator()(T& t)
    {
    	typedef decltype(src.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>()) stype;

    	// In case of layout switch use this
    	//meta_copy_d<stype,dtype>::meta_copy_d_(src.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(),dst.template get<T::value>());

    	meta_copy<typename std::remove_reference<stype>::type>::meta_copy_(src.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(),dst.template get<T::value>());
    }
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy the selected properties
 *
 * \tparam v_src source object
 * \tparam d_src destination object
 *
 */

template<typename v_src,typename v_dst, int... prp>
struct object_si_d_f
{
	// Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	// Source object
	const v_src & src;

	// Destination object
	v_dst & dst;

	/*! \brief Constructor
	 *
	 * \param src source object
	 * \param dst destination object
	 *
	 */
	object_si_d_f(const v_src & src, v_dst & dst)
	:src(src),dst(dst)
	{
	};

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
    	typedef typename boost::mpl::at<typename v_dst::type,typename boost::mpl::int_<T::value>>::type ctype;

    	meta_copy<ctype>::meta_copy_(boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(src.data),boost::fusion::at_c<T::value>(dst.data));
    }
};

#define OBJ_ENCAP 1
#define OBJ_NORMAL 2

/*! \brief It copy the properties from one object to another
 *
 * Stub object
 *
 * \see object_copy<v_src,v_dst,NORMAL,prp...> object_copy<v_src,v_dst,ENCAP,prp...>
 *
 *
 */
template<typename v_src, typename v_dst,int type_copy, int... prp>
struct object_si_d
{
	inline object_si_d(const v_src & vs, v_dst & vd)
	{
		std::cerr << "Error object_copy: " << __FILE__ << " " << __LINE__ << "\n";
	};
};

/*! \brief It copy the properties from one object to another
 *
 * Given a set of properties for the source (0,1,3) it copy that property to the destination properties
 * (0,1,2)
 *
 * ### Object copy example
 * \snippet util_test.hpp object copy example
 *
 */
template<typename v_src, typename v_dst, int... prp>
struct object_si_d<v_src,v_dst,OBJ_NORMAL,prp...>
{
	inline object_si_d(const v_src && vs, v_dst && vd)
	{
		object_si_d_f<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,v_dst::max_prop> >(obj);
	}

	inline object_si_d(const v_src & vs, v_dst & vd)
	{
		object_si_d_f<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,v_dst::max_prop> >(obj);
	}
};

/*! \brief It copy the properties from one object to another
 *
 * Given a set of properties for the source object (0,1,3) it copy that properties
 * to the destination object properties (0,1,2)
 *
 * For object we mean an object that follow the OpenFPM data structure format, see openFPM_data wiki
 * for more information
 *
 * ## Create a compile-time object and copy *from* the selected properties
 * \snippet util_test.hpp object copy example
 * ## Create a compile-time Encap object and copy *from* the selected properties
 * \snippet util_test.hpp object copy encap example
 *
 */
template<typename v_src, typename v_dst, int... prp>
struct object_si_d<v_src,v_dst,OBJ_ENCAP,prp...>
{
	__device__ __host__ inline object_si_d(const v_src && vs, v_dst && vd)
	{
		object_si_d_e<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,v_dst::max_prop> >(obj);
	}

	__device__ __host__ inline object_si_d(const v_src & vs, v_dst & vd)
	{
		object_si_d_e<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,v_dst::max_prop> >(obj);
	}
};

template<typename v_src, typename v_dst>
struct object_si_d<v_src,v_dst,OBJ_ENCAP>
{
	__device__ __host__  inline object_si_d(const v_src && vs, v_dst && vd)
	{
	}

	__device__ __host__  inline object_si_d(const v_src & vs, v_dst & vd)
	{
	}
};

#endif /* VECTOR_PROP_COPY_HPP_ */
