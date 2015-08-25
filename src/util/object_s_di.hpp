/*
 * object_write.hpp
 *
 *  Created on: Jun 11, 2015
 *      Author: i-bird
 */

#ifndef OBJECT_WRITE_HPP_
#define OBJECT_WRITE_HPP_


#include "for_each_ref.hpp"
#include "util/variadic_to_vmpl.hpp"
#include "util/meta_copy.hpp"
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
struct object_s_di_e
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
	object_s_di_e(const v_src & src, v_dst & dst)
	:src(src),dst(dst)
	{
	};

#ifdef DEBUG
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	object_s_di_e(const v_src && src, v_dst & dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
    	typedef typename boost::mpl::at<typename v_dst::type,typename boost::mpl::int_<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>>::type ctype;

    	meta_copy<ctype>(src.template get<T::value>(),dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>());
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
struct object_s_di_f
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
	object_s_di_f(const v_src & src, v_dst & dst)
	:src(src),dst(dst)
	{
	};

#ifdef DEBUG
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	object_s_di_f(const v_src && src, v_dst & dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
    	typedef typename boost::mpl::at<typename v_dst::type,typename boost::mpl::int_<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>>::type ctype;

    	meta_copy<ctype>(boost::fusion::at_c<T::value>(src.data),boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst.data));
    }
};

#define ENCAP 1
#define NORMAL 2

/*! \brief It copy the properties from one object to another
 *
 * Stub object
 *
 * \see object_copy<v_src,v_dst,NORMAL,prp...> object_copy<v_src,v_dst,ENCAP,prp...>
 *
 *
 */
template<typename v_src, typename v_dst,int type_copy, int... prp>
struct object_s_di
{
	inline object_s_di(const v_src & vs, v_dst & vd)
	{
		std::cerr << "Error object_copy: " << __FILE__ << " " << __LINE__ << "\n";
	};
};

/*! \brief Given a set of properties for the destination (0,1,3,5) it copy the source properties (0,1,2,3)
 *
 *
 * ### Object copy example
 * \snippet util_test.hpp object copy example
 *
 */
template<typename v_src, typename v_dst, int... prp>
struct object_s_di<v_src,v_dst,NORMAL,prp...>
{
	inline object_s_di(const v_src && vs, v_dst && vd)
	{
		object_s_di_f<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}

	inline object_s_di(const v_src & vs, v_dst & vd)
	{
		object_s_di_f<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}
};

/*! \brief It copy the properties from one object to another
 *
 * Given a set of properties for the destination object (0,1,3) it copy that properties
 * to the source object properties (0,1,2)
 *
 * For object we mean an object that follow the OpenFPM data structure format, see openFPM_data wiki
 * for more information
 *
 * ## Create a compile-time object and copy *to* the selected properties
 * \snipper util_test.hpp object write example
 * ## Create a compile-time Encap object and copy *to* the selected properties
 * \snipper util_test.hpp object write encap example
 *
 */
template<typename v_src, typename v_dst, int... prp>
struct object_s_di<v_src,v_dst,ENCAP,prp...>
{
	inline object_s_di(const v_src && vs, v_dst && vd)
	{
		object_s_di_e<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}

	inline object_s_di(const v_src & vs, v_dst & vd)
	{
		object_s_di_e<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}
};




#endif /* OBJECT_WRITE_HPP_ */
