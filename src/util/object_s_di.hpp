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
#include "util/copy_compare/meta_copy.hpp"
#include <boost/mpl/range_c.hpp>
#include <boost/fusion/include/size.hpp>

template <typename> struct Debug;

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy the selected properties
 *
 * \tparam v_src source object
 * \tparam d_src destination object
 * \tparam prp properties
 *
 */

template<typename v_src,typename v_dst, int... prp>
struct object_s_di_e_cnk
{
	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	//! Source object
	const v_src & src;

	//! Destination object
	v_dst & dst;

	//! element id
	size_t sub_id;

	/*! \brief Constructor
	 *
	 * \param src source object
	 * \param dst destination object
	 *
	 */
	object_s_di_e_cnk(const v_src & src, v_dst & dst,size_t sub_id)
	:src(src),dst(dst),sub_id(sub_id)
	{
	};

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>())>::type::value_type copy_rtype;

    	meta_copy<copy_rtype>::meta_copy_(src.template get<T::value>(),dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>()[sub_id]);
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
 * \tparam prp properties
 *
 */

template<typename v_src,typename v_dst, int... prp>
struct object_s_di_e
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
	object_s_di_e(const v_src & src, v_dst & dst)
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
	object_s_di_e(const v_src && src, v_dst & dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the functor for each member
    template<typename T>
    __device__ __host__ void operator()(T& t)
    {
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>())>::type copy_dtype;
		typedef typename std::remove_reference<decltype(src.template get<T::value>())>::type copy_stype;

    	meta_copy_d<copy_stype,copy_dtype>::meta_copy_d_(src.template get<T::value>(),dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>());
    }
};


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
struct object_s_di_e_op
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
	object_s_di_e_op(const v_src & src, v_dst & dst)
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
	object_s_di_e_op(v_src && src, v_dst & dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>())>::type copy_dtype;
		typedef typename std::remove_reference<decltype(src.template get<T::value>())>::type copy_stype;

    	meta_copy_op_d<op,copy_stype,copy_dtype>::meta_copy_op_d_(src.template get<T::value>(),dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>());
    }
};

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
struct object_s_di_e_op_cnk
{
	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	//! Source object
	const v_src & src;

	//! Destination object
	v_dst & dst;

	//! element id
	size_t sub_id;

	/*! \brief Constructor
	 *
	 * \param src source object
	 * \param dst destination object
	 *
	 */
	object_s_di_e_op_cnk(const v_src & src, v_dst & dst,size_t sub_id)
	:src(src),dst(dst),sub_id(sub_id)
	{
	};

#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	object_s_di_e_op_cnk(v_src && src, v_dst & dst, size_t sub_id)
	:src(src),dst(dst),sub_id(sub_id)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>()[sub_id])>::type copy_rtype;

    	meta_copy_op<op,copy_rtype>::meta_copy_op_(src.template get<T::value>(),dst.template get<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>()[sub_id]);
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

    	meta_copy<ctype>::meta_copy_(boost::fusion::at_c<T::value>(src.data),boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst.data));
    }
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy the selected properties
 *
 * \tparam op operation to apply
 * \tparam v_src source object
 * \tparam d_src destination object
 *
 */

template<template<typename,typename> class op, typename v_src,typename v_dst, int... prp>
struct object_s_di_f_op
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
	object_s_di_f_op(const v_src & src, v_dst & dst)
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
	object_s_di_f_op(v_src && src, v_dst & dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
    	typedef typename boost::mpl::at<typename v_dst::type,typename boost::mpl::int_<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>>::type ctype;

    	meta_copy_op<op,ctype>::meta_copy_op_(boost::fusion::at_c<T::value>(src.data),boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst.data));
    }
};

#define OBJ_ENCAP 1
#define OBJ_NORMAL 2
#define OBJ_ENCAP_CHUNKING 3

/*! \brief It copy the properties from one object to another
 *
 * Stub object
 *
 * \see object_s_di<v_src,v_dst,OBJ_ENCAP,prp...> object_copy<v_src,v_dst,OBJ_ENCAP,prp...>
 *
 *
 */
template<typename v_src, typename v_dst,int type_copy, int... prp>
struct object_s_di
{
	/*! \brief Stub method
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
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
struct object_s_di<v_src,v_dst,OBJ_NORMAL,prp...>
{
	/*! \brief Implementation of the copy
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	inline object_s_di(const v_src && vs, v_dst && vd)
	{
		object_s_di_f<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}

	/*! \brief Implementation of the copy
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
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
 * \snippet util_test.hpp object write example
 * ## Create a compile-time Encap object and copy *to* the selected properties
 * \snippet util_test.hpp object write encap example
 *
 */
template<typename v_src, typename v_dst, int... prp>
struct object_s_di<v_src,v_dst,OBJ_ENCAP,prp...>
{
	/*! \brief Implementation of the copy
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	__host__ __device__ inline object_s_di(const v_src & vs, v_dst && vd)
	{
		object_s_di_e<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,(int)sizeof...(prp)> >(obj);
	}

	/*! \brief Implementation of the copy
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	__host__ __device__ inline object_s_di(const v_src & vs, v_dst & vd)
	{
		object_s_di_e<v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}
};



template<typename v_src, typename v_dst, int... prp>
struct object_s_di<v_src,v_dst,OBJ_ENCAP_CHUNKING,prp...>
{
	/*! \brief Implementation of the copy
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	inline object_s_di(const v_src & vs, v_dst && vd, size_t sub_id)
	{
		object_s_di_e_cnk<v_src,v_dst,prp...> obj(vs,vd,sub_id);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}

	/*! \brief Implementation of the copy
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	inline object_s_di(const v_src & vs, v_dst & vd, size_t sub_id)
	{
		object_s_di_e_cnk<v_src,v_dst,prp...> obj(vs,vd,sub_id);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}
};

////////////////////////// WITH OPERATION VERSION

/*! \brief It copy the properties from one object to another applying an operation
 *
 * Stub object
 *
 * \see object_s_di_op<v_src,v_dst,OBJ_ENCAP,prp...> object_copy<v_src,v_dst,OBJ_ENCAP,prp...>
 *
 *
 */
template<template<typename,typename> class op, typename v_src, typename v_dst,int type_copy, int... prp>
struct object_s_di_op
{
	/*! \brief Stub method
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	inline object_s_di_op(const v_src & vs, v_dst & vd)
	{
		std::cerr << "Error object_copy: " << __FILE__ << " " << __LINE__ << "\n";
	};
};

/*! \brief Given a set of properties for the destination (0,1,3,5)
 *         it copy the source properties (0,1,2,3) applying an operation
 *
 *
 * ### Object copy example
 * \snippet util_test.hpp object copy example
 *
 */
template<template<typename,typename> class op, typename v_src, typename v_dst, int... prp>
struct object_s_di_op<op,v_src,v_dst,OBJ_NORMAL,prp...>
{
	/*! \brief Implementation of the copy with operation
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	inline object_s_di_op(const v_src && vs, v_dst && vd)
	{
		object_s_di_f_op<op,v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}

	/*! \brief Implementation of the copy with operation
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	inline object_s_di_op(const v_src & vs, v_dst & vd)
	{
		object_s_di_f_op<op,v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}
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
struct object_s_di_op<op, v_src,v_dst,OBJ_ENCAP,prp...>
{
	/*! \brief Implementation of the copy with operation
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	__device__ __host__ inline object_s_di_op(const v_src & vs, v_dst && vd)
	{
		object_s_di_e_op<op,v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}

	/*! \brief Implementation of the copy with operation
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	__device__ __host__ inline object_s_di_op(const v_src & vs, v_dst & vd)
	{
		object_s_di_e_op<op,v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}
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
struct object_s_di_op<op, v_src,v_dst,OBJ_ENCAP_CHUNKING,prp...>
{
	/*! \brief Implementation of the copy with operation
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	inline object_s_di_op(const v_src & vs, v_dst && vd)
	{
		object_s_di_e_op_cnk<op,v_src,v_dst,prp...> obj(vs,vd);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}

	/*! \brief Implementation of the copy with operation
	 *
	 * \param vs source object
	 * \param vd destination object
	 *
	 */
	inline object_s_di_op(const v_src & vs, v_dst & vd, size_t sub_id)
	{
		object_s_di_e_op_cnk<op,v_src,v_dst,prp...> obj(vs,vd,sub_id);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(obj);
	}
};

#endif /* OBJECT_WRITE_HPP_ */
