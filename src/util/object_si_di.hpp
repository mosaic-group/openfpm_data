/*
 * object_si_di.hpp
 *
 *  Created on: Jun 11, 2015
 *      Author: i-bird
 */

#ifndef OBJECT_SI_DI_HPP_
#define OBJECT_SI_DI_HPP_


#include "for_each_ref.hpp"
#include "util/variadic_to_vmpl.hpp"
#include "util/copy_compare/meta_copy.hpp"
#include <boost/mpl/range_c.hpp>
#include <boost/fusion/include/size.hpp>



/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy the selected properties applying an operation
 *
 * \tparam op operation
 * \tparam v_src source object
 * \tparam d_src destination object
 * \tparam prp properties
 *
 */

template<template<typename,typename> class op, typename v_src,typename v_dst, int... prp>
struct object_si_di_e_op
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
	object_si_di_e_op(const v_src & src, v_dst & dst)
	:src(src),dst(dst)
	{
	};

#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	object_si_di_e_op(v_src && src, v_dst & dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>())>::type copy_dtype;
		typedef typename std::remove_reference<decltype(src.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>())>::type copy_stype;

    	meta_copy_op_d<op,copy_stype,copy_dtype>::meta_copy_op_d_(src.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(),dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>());
    }
};





#define OBJ_ENCAP 1

////////////////////////// WITH OPERATION VERSION

/*! \brief It copy the properties from one object to another applying an operation
 *
 * Stub object
 *
 * \see object_si_di_op<v_src,v_dst,OBJ_ENCAP,prp...> object_copy<v_src,v_dst,OBJ_ENCAP,prp...>
 *
 *
 */
template<template<typename,typename> class op, typename v_src, typename v_dst,int type_copy, int... prp>
struct object_si_di_op
{
	/*! \brief Stub method
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	inline object_si_di_op(const v_src & vs, v_dst & vd)
	{
		std::cerr << "Error object_copy: " << __FILE__ << " " << __LINE__ << "\n";
	};
};




/*! \brief It copy the properties from one object to another applying an operation
 *
 * Given a set of properties for the destination object (0,1,3) it copy that properties
 * to the source object properties (0,1,2) applying an operation
 *
 * For object we mean an object that follow the OpenFPM data structure format, see openFPM_data wiki
 * for more information
 *
 * ## Create a compile-time object and copy *to* the selected properties applying an operation
 * \snippet util_test.hpp object write example with op
 * ## Create a compile-time Encap object and copy *to* the selected properties applying an operation
 * \snippet util_test.hpp object write encap example with op
 *
 */
template<template<typename,typename> class op, typename v_src, typename v_dst, int... prp>
struct object_si_di_op<op, v_src,v_dst,OBJ_ENCAP,prp...>
{
	/*! \brief Implementation of the copy with operation
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	__device__ __host__ inline object_si_di_op(const v_src & vs, v_dst && vd)
	{
		object_si_di_e_op<op,v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}

	/*! \brief Implementation of the copy with operation
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	__device__ __host__ inline object_si_di_op(const v_src & vs, v_dst & vd)
	{
		object_si_di_e_op<op,v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}
};



#endif /* OBJECT_WRITE_HPP_ */
