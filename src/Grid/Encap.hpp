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
	//! object we have to store
	const e_src & src;
	e_dst & dst;


	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param key which element we are modifying
	 * \param grid_dst grid we are updating
	 * \param obj object we have to set in grid_dst (encapsulated)
	 *
	 */
	inline copy_cpu_encap_encap(const e_src & src, e_dst & dst)
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
	 */
	inline copy_cpu_encap_encap(const e_src && src, const e_dst && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		// This is the type of the object we have to copy
		typedef typename boost::fusion::result_of::at_c<typename e_src::type,T::value>::type copy_type;

		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<copy_type>::type copy_rtype;

		meta_copy<copy_rtype> cp(src.template get<T::value>(),dst.template get<T::value>());
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
struct compare_cpu_encap_encap
{
	//! object we have to compare
	const e_src & src;
	const e_dst & dst;


	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param key which element we are modifying
	 * \param grid_dst grid we are updating
	 * \param obj object we have to set in grid_dst (encapsulated)
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
	 */
	inline compare_cpu_encap_encap(const e_src && src, const e_dst && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		// This is the type of the object we have to copy
		typedef typename boost::fusion::result_of::at_c<typename e_src::type,T::value>::type copy_type;

		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<copy_type>::type copy_rtype;

		meta_compare<copy_rtype> cp(src.template get<T::value>(),dst.template get<T::value>());
	}
};

/*! \brief This class is an helper to get the return type for get method for each property
 *
 * This class is an helper to get the return type for get method for each property
 *
 * \param p id of the property
 * \param T original boost fusion vector, T is suppose to be a boost::fusion::vector<memory_c<...>,.....>
 *
 */

template<unsigned int p,typename mem>
struct type_cpu_prop
{
	//! return a boost::fusion::vector<memory_c<....>....>
	typedef typename mem::vtype vtype;
	//! return a memory_c<...>
	typedef typename boost::fusion::result_of::at< vtype,boost::mpl::int_<p> >::type type;
};

/*! \brief This class is an helper to get the return type of get for each property
 *
 * This class is an helper to get the return type of get for each property
 *
 * \param p id of the property
 * \param T original boost fusion vector, T is suppose to be a boost::fusion::vector<memory_c<...>,.....>
 *
 */

template<unsigned int p,typename Mem>
struct type_gpu_prop
{
	//! return a boost::fusion::vector<memory_c<....>....>
	typedef Mem vtype;
	//! return a memory_c<...>
	typedef typename boost::fusion::result_of::at< vtype,boost::mpl::int_<p> >::type rtype;
	//! remove the reference
	typedef typename boost::remove_reference<rtype>::type mtype;
	//! get the base type that the buffer is storing
	typedef typename mtype::type type;
};



/*! \brief this structure encapsulate an object of the grid
 *
 * This structure encapsulate an object of the grid
 * It give the possibility to select the property in a secondary moment
 *
 * \param dim Dimensionality of the grid
 * \param T type of object the grid store
 * \param Mem suppose to be a boost::fusion::vector of arrays
 *
 */

template<unsigned int dim,typename T,typename Mem>
class encapc
{
public:
	typedef typename T::type type;

private:

	type & data_c;

public:

	typedef int yes_i_am_encap;

	typedef T T_type;

	static const int max_prop = T::max_prop;

	// constructor require a key and a memory data
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
	template <unsigned int p> inline typename type_cpu_prop<p,Mem>::type get()
	{
#ifdef SE_CLASS2
		check_valid(&boost::fusion::at_c<p>(data_c),sizeof(typename type_cpu_prop<p,Mem>::type));
#endif
		return boost::fusion::at_c<p>(data_c);
	}

	/*! \brief access the data but return a constant reference
	 *
	 * \return the reference
	 *
	 */
	template <unsigned int p> inline const typename type_cpu_prop<p,Mem>::type get() const
	{
#ifdef SE_CLASS2
		check_valid(&boost::fusion::at_c<p>(data_c),sizeof(typename type_cpu_prop<p,Mem>::type));
#endif
		return boost::fusion::at_c<p>(data_c);
	}

	// access the data
	template <unsigned int p> inline void set(typename type_cpu_prop<p,Mem>::type & ele)
	{
#ifdef SE_CLASS2
			check_valid(&boost::fusion::at_c<p>(data_c),sizeof(typename type_cpu_prop<p,T>::type));
#endif
			return boost::fusion::at_c<p>(data_c) = ele;
	}

	/*! \brief Assignment
	 *
	 * \param ec encapsulator
	 *
	 */
	inline encapc<dim,T,Mem> & operator=(const encapc<dim,T,Mem> & ec)
	{
		copy_cpu_encap_encap<encapc<dim,T,Mem>,encapc<dim,T,Mem>> cp(ec,*this);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

		return *this;
	}

	/*! \brief Assignment
	 *
	 * \param ec encapsulator
	 *
	 */
	inline encapc<dim,T,Mem> & operator=(const T & obj)
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

	// access the data
	template <unsigned int p> void get()
	{}

	// access the data
	template <unsigned int p, typename S> void set(S & ele)
	{}
};

/*! \brief this structure encapsulate an object of the grid
 *
 * This structure encapsulate an object of the grid
 * It give the possibility to select the property in a secondary moment
 *
 *	\param dim Dimensionality of the grid
 *	\param T type of object the grid store
 *	\param Mem interface used to allocate memory
 *
 */

template<unsigned int dim,typename T,typename Mem>
class encapg
{
	// constructor require a key
	Mem & data;
	size_t k;

public:

	typedef int yes_i_am_encap;

	typedef T T_type;

	// constructor require a key and a memory data
	encapg(Mem & data, size_t k)
	:data(data),k(k)
	{}

	// access the data
	template <unsigned int p> typename type_gpu_prop<p,Mem>::type::reference get()
	{
		return boost::fusion::at_c<p>(data).mem_r->operator[](k);
	}

	// access the data
	template <unsigned int p> const typename type_gpu_prop<p,Mem>::type::reference get() const
	{
		return boost::fusion::at_c<p>(data).mem_r->operator[](k);
	}
};

#include "util/common.hpp"

template<typename T, typename Sfinae = void>
struct is_encap: std::false_type {};


/*! \brief is_encap check if the type is an encap type
 *
 * ### Example
 *
 * \snippet util.hpp Check is_encap
 *
 * return true if T is an encap
 *
 */
template<typename T>
struct is_encap<T, typename Void< typename T::yes_i_am_encap>::type> : std::true_type
{};

#endif /* ENCAP_HPP_ */
