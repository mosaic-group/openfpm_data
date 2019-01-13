/*
 * Encap.hpp
 *
 *  Created on: Feb 2, 2015
 *      Author: i-bird
 */

#ifndef ENCAP_HPP_
#define ENCAP_HPP_

#include "util/for_each_ref.hpp"
#include "util/copy_compare/meta_copy.hpp"
#include "boost/mpl/range_c.hpp"
#include <boost/fusion/container/vector.hpp>
#ifdef SE_CLASS2
#include "Memleak_check.hpp"
#endif
#include "util/se_util.hpp"
#include "util/copy_compare/copy_fusion_vector.hpp"
#include "util/copy_compare/compare_fusion_vector.hpp"
#include "memory_ly/memory_conf.hpp"

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<typename e_src, typename e_dst, unsigned int ... prp>
struct copy_cpu_encap_encap_prp
{
	//! encapsulated source object
	e_src & src;
	//! encapsulated destination object
	e_dst & dst;

	typedef typename to_boost_vmpl<prp...>::type vprp;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	__device__ __host__ inline copy_cpu_encap_encap_prp(e_src & src, e_dst & dst)
	:src(src),dst(dst)
	{
#if defined(SE_CLASS1) && !defined(__NVCC__)
		// e_src and e_dst must have the same number of properties

		if (e_src::T_type::max_prop != e_dst::T_type::max_prop)
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the number of properties between src and dst must match";
#endif
	};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_cpu_encap_encap_prp(const e_src && src, const e_dst && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<vprp,T>::type prp_cp;

		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get<prp_cp::value>())>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(src.template get<prp_cp::value>(),dst.template get<prp_cp::value>());
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<typename e_src, typename e_dst>
struct copy_cpu_encap_encap
{
	//! encapsulated source object
	const e_src & src;
	//! encapsulated destination object
	e_dst & dst;


	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	__device__ __host__ inline copy_cpu_encap_encap(const e_src & src, e_dst & dst)
	:src(src),dst(dst)
	{
#if defined(SE_CLASS1) && !defined(__NVCC__)
		// e_src and e_dst must have the same number of properties

		if (e_src::T_type::max_prop != e_dst::T_type::max_prop)
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the number of properties between src and dst must match";
#endif
	};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_cpu_encap_encap(const e_src && src, const e_dst && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t) const
	{
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get<T::value>())>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(src.template get<T::value>(),dst.template get<T::value>());
	}
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<typename e_src, typename e_dst>
struct copy_cpu_encap_encap_general
{
	//! encapsulated source object
	const e_src & src;
	//! encapsulated destination object
	e_dst & dst;


	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	__device__ __host__ inline copy_cpu_encap_encap_general(const e_src & src, e_dst & dst)
	:src(src),dst(dst)
	{
#if defined(SE_CLASS1) && !defined(__NVCC__)
		// e_src and e_dst must have the same number of properties

		if (e_src::T_type::max_prop != e_dst::T_type::max_prop)
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the number of properties between src and dst must match";
#endif
	};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_cpu_encap_encap_general(const e_src && src, const e_dst && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t) const
	{
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get<T::value>())>::type copy_dtype;

		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(src.template get<T::value>())>::type copy_stype;

		meta_copy_d<copy_stype,copy_dtype>::meta_copy_d_(src.template get<T::value>(),dst.template get<T::value>());
	}
};



/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<typename e_src>
struct copy_cpu_encap_single
{
	//! encapsulated source object
	const e_src & src;
	//! encapsulated destination object
	e_src & dst;


	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	__device__ __host__ inline copy_cpu_encap_single(const e_src & src, e_src & dst)
	:src(src),dst(dst)
	{
#if defined(SE_CLASS1) && !defined(__NVCC__)
		// e_src and e_dst must have the same number of properties

		if (e_src::T_type::max_prop != e_src::T_type::max_prop)
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the number of properties between src and dst must match";
#endif
	};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_cpu_encap_single(const e_src && src, const e_src && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t) const
	{
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get<T::value>())>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(src.template get<T::value>(),dst.template get<T::value>());
	}
};


///////////////////////////////////

//! It copy two encap object
template<template<typename,typename> class op,typename e_src, typename e_dst, int ... prp>
struct copy_cpu_encap_encap_op_prp
{
	//! encapsulated object source
	const e_src & src;
	//! encapsulated object destination
	e_dst & dst;


	/*! \brief constructor
	 *
	 *
	 * \param src source encapsulated object
	 * \param dst destination encapsulated object
	 *
	 */
	inline copy_cpu_encap_encap_op_prp(const e_src & src, e_dst & dst)
	:src(src),dst(dst)
	{
#ifdef SE_CLASS1
		// e_src and e_dst must have the same number of properties

		if (e_src::T_type::max_prop != e_dst::T_type::max_prop)
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the number of properties between src and dst must match";
#endif
	};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_cpu_encap_encap_op_prp(const e_src && src, const e_dst && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		// Convert variadic to boost::vector
		typedef typename boost::mpl::vector_c<unsigned int,prp...> prpv;

		// element id to copy applying an operation
		typedef typename boost::mpl::at<prpv,T>::type ele_cop;

		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(src.template get< ele_cop::value >())>::type copy_rtype;

		meta_copy_op<op,copy_rtype>::meta_copy_op_(src.template get< ele_cop::value >(),dst.template get< ele_cop::value >());
	}
};



///////////////////////////////

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */

template<typename e_src, typename e_dst>
struct compare_cpu_encap_encap
{
	//! object 1 we have to compare
	const e_src & src;
	//! object 2 we have to compare
	const e_dst & dst;


	/*! \brief constructor
	 *
	 *
	 * \param src encapsulated object1
	 * \param dst encapsulated object2
	 *
	 */
	inline compare_cpu_encap_encap(const e_src & src, const e_dst & dst)
	:src(src),dst(dst)
	{
#ifdef SE_CLASS1
		// e_src and e_dst must have the same number of properties

		if (e_src::max_prop != e_dst::max_prop)
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the number of properties between src and dst must match";
#endif
	};


#ifdef SE_CLASS1

	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 * \param e_src encapsulated object1
	 * \param e_dst encapsulated object2
	 *
	 */
	inline compare_cpu_encap_encap(const e_src && src, const e_dst && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	/*!  \brief It call the copy function for each property
	 *
	 * \param t each member
	 *
	 */
	template<typename T>
	inline void operator()(T& t) const
	{
		//! This is the type of the object we have to copy
		typedef typename boost::fusion::result_of::at_c<typename e_src::type,T::value>::type copy_type;

		//! Remove the reference from the type to copy
		typedef typename boost::remove_reference<copy_type>::type copy_rtype;

		meta_compare<copy_rtype> cp(src.template get<T::value>(),dst.template get<T::value>());
	}
};

/*! Stub encap
 *
 * \tparam dim dimension
 * \tparam T object it encapsulate
 * \tparam layout
 *
 */
template<unsigned int dim,typename T,typename layout>
class encapc
{

};

/*! \brief this structure encapsulate an object of the grid
 *
 * This structure encapsulate an object of the grid
 * It give the possibility to select the property in a secondary moment
 *
 * Can be thought as a reference to an object of the grid. So every time we
 * use the term encapsulated object we mean reference to object
 *
 * \note A vector is a 1D grid
 *
 * \param dim Dimensionality of the grid
 * \param T type of object the grid store
 *
 *
 */

template<unsigned int dim,typename T>
class encapc<dim,T,typename memory_traits_lin<T>::type >
{
public:

	//! object type the encap object encapsulate
	typedef typename T::type type;

private:

	//! reference to the encapsulated object
	type & data_c;

	//! layout of the encapsulated object
	typedef typename memory_traits_lin<T>::type Mem;

	//! layout of the encapsulated object
	typedef typename memory_traits_inte<T>::type Mem2;

public:

	//! indicate the it is an encapsulated object
	typedef int yes_i_am_encap;

	//! original object
	typedef T T_type;

	//! number of properties
	static const int max_prop = T::max_prop;

	//! constructor from a reference object
	inline encapc(type & data_c)
	:data_c(data_c)
	{}

	/*! \brief Return the address of the base
	 *
	 * \return the address of the data encapsulated
	 *
	 */
	inline type * operator&()
	{
		return &data_c;
	}

	/*! \brief access the data
	 *
	 * \return the reference
	 *
	 */
	template <unsigned int p, typename r_type=decltype(boost::fusion::at_c<p>(data_c))>
	inline r_type get()
	{
#ifdef SE_CLASS2
		check_valid(&boost::fusion::at_c<p>(data_c),sizeof(typename boost::mpl::at<type,boost::mpl::int_<p>>::type));
#endif
		return boost::fusion::at_c<p>(data_c);
	}

	/*! \brief access the data but return a constant reference
	 *
	 * \return the reference
	 *
	 */
	template <unsigned int p, typename r_type=decltype(boost::fusion::at_c<p>(data_c))>
	inline const r_type get() const
	{
#ifdef SE_CLASS2
		check_valid(&boost::fusion::at_c<p>(data_c),sizeof(typename boost::mpl::at<type,boost::mpl::int_<p>>::type));
#endif
		return boost::fusion::at_c<p>(data_c);
	}

	/*! \brief Set one property of the encapsulated object
	 *
	 * \tparam p property to set
	 *
	 * \param ele value to set
	 *
	 */
	template <unsigned int p> inline void set(decltype(boost::fusion::at_c<p>(data_c)) & ele)
	{
#ifdef SE_CLASS2
			check_valid(&boost::fusion::at_c<p>(data_c),sizeof(typename boost::mpl::at<type,boost::mpl::int_<p>>::type));
#endif
			return boost::fusion::at_c<p>(data_c) = ele;
	}

	/*! \brief Set one property of the encapsulated object
	 *
	 * \tparam dim2 dimensionality of the multy-array
	 * \tparam p property to set
	 *
	 * \param ec value to set as encapsulated object
	 *
	 * \return itself
	 *
	 */
	template<unsigned int dim2> inline encapc<dim,T,Mem> & set(const encapc<dim2,T,Mem> & ec)
	{
		copy_cpu_encap_encap<encapc<dim2,T,Mem>,encapc<dim,T,Mem>> cp(ec,*this);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

		return *this;
	}

	/*! \brief Assignment
	 *
	 * \param ec object encapsulated to copy
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ inline encapc<dim,T,Mem> & operator=(const encapc<dim,T,Mem> & ec)
	{
		copy_cpu_encap_single<encapc<dim,T,Mem>> cp(ec,*this);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

		return *this;
	}

	/*! \brief Assignment
	 *
	 * \param ec object encapsulated to copy
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ inline encapc<dim,T,Mem> & operator=(const encapc<dim,T,Mem2> & ec)
	{
		copy_cpu_encap_encap_general<encapc<dim,T,Mem2>,encapc<dim,T,Mem>> cp(ec,*this);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

		return *this;
	}

	/*! \brief Assignment
	 *
	 * \param obj object to copy
	 *
	 * \return itself
	 *
	 */
	__device__  __host__ inline encapc<dim,T,Mem> & operator=(const T & obj)
	{
		copy_fusion_vector<typename T::type> cp(obj.data,data_c);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

		return *this;
	}

	/*! \brief Compare
	 *
	 * \param ec encapsulator
	 *
	 * \return true if the two encap store the same information
	 *
	 */
	inline bool operator==(const encapc<dim,T,Mem> & ec) const
	{
		compare_fusion_vector<typename T::type> cp(ec.data_c,data_c);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

		return cp.result();
	}

	/*! \brief Compare
	 *
	 * \param ec encapsulator
	 *
	 * \return true if the two encap store different information
	 *
	 */
	inline bool operator!=(const encapc<dim,T,Mem> & ec) const
	{
		return ! this->operator==(ec);
	}
};

/*! \brief this structure specialize the class for a void object or null
 *
 * \param dim Dimensionality of the grid
 * \param Mem suppose to be a boost::fusion::vector of arrays
 *
 */

template<unsigned int dim,typename Mem>
class encapc<dim,void,Mem>
{
public:

	// constructor require a key and a memory data
	encapc()
	{}

	//! access the data
	template <unsigned int p> void get()
	{}

	//! set the data
	template <unsigned int p, typename S> void set(S & ele)
	{}
};

/*! \brief this structure encapsulate an object of the grid
 *
 * This structure encapsulate an object of the grid
 * It give the possibility to select the property in a secondary moment
 *
 * Can be thought as a reference to an object of the grid. So every time we
 * use the term encapsulated object we mean reference to object
 *
 * \note A vector is a 1D grid
 *
 *	\param dim Dimensionality of the grid
 *	\param T type of object the grid store
 *
 */
template<unsigned int dim,typename T>
class encapc<dim,T,typename memory_traits_inte<T>::type>
{
	//! type of layout
	typedef typename memory_traits_inte<T>::type Mem;

	//! layout of the encapsulated object
	typedef typename memory_traits_lin<T>::type Mem2;

	//! reference to the encapsulated object
	Mem & data;

	//! element id
	size_t k;

public:

	//! Original list if types
	typedef typename T::type type;

	//! indicate it is an encapsulated object
	typedef int yes_i_am_encap;

	//! original object type
	typedef T T_type;

	//! number of properties
	static const int max_prop = T::max_prop;

	//! constructor require a key and a memory data
	encapc(typename memory_traits_inte<T>::type & data, size_t k)
	:data(data),k(k)
	{}

	/*! \brief Access the data
	 *
	 * \tparam p property selected
	 *
	 * \return The reference of the data
	 *
	 */
	template <unsigned int p> __device__ __host__  auto get() -> decltype(boost::fusion::at_c<p>(data).mem_r.operator[](k))
	{
		return boost::fusion::at_c<p>(data).mem_r.operator[](k);
	}

	/*! \brief Access the data
	 *
	 * \tparam p property selected
	 *
	 * \return The reference of the data
	 *
	 */
	template <unsigned int p> __device__ __host__ auto get() const -> decltype(boost::fusion::at_c<p>(data).mem_r.operator[](k))
	{
		return boost::fusion::at_c<p>(data).mem_r.operator[](k);
	}

	/*! \brief Assignment
	 *
	 * \param ec encapsulator
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ inline encapc<dim,T,Mem> & operator=(const encapc<dim,T,Mem> & ec)
	{
		copy_cpu_encap_single<encapc<dim,T,Mem>> cp(ec,*this);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

		return *this;
	}

	/*! \brief Assignment
	 *
	 * \param ec encapsulator
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ inline encapc<dim,T,Mem> & operator=(const encapc<dim,T,Mem2> & ec)
	{
		copy_cpu_encap_encap_general<encapc<dim,T,Mem2>,encapc<dim,T,Mem>> cp(ec,*this);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

		return *this;
	}

	/*! \brief Assignment
	 *
	 * \param obj object to copy
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ inline encapc<dim,T,Mem> & operator=(const T & obj)
	{
		copy_fusion_vector_encap<typename T::type,decltype(*this)> cp(obj.data,*this);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

		return *this;
	}
};

#include "util/common.hpp"

template<typename T, typename Sfinae = void>
struct is_encap: std::false_type {};


/*! \brief is_encap check if the type is an encap type
 *
 * ### Example
 *
 * \snippet util_test.hpp Check is_encap
 *
 * return true if T is an encap
 *
 */
template<typename T>
struct is_encap<T, typename Void< typename T::yes_i_am_encap>::type> : std::true_type
{};

#endif /* ENCAP_HPP_ */
