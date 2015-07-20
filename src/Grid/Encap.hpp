/*
 * Encap.hpp
 *
 *  Created on: Feb 2, 2015
 *      Author: i-bird
 */

#ifndef ENCAP_HPP_
#define ENCAP_HPP_

//#include "grid_sm.hpp"

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

	type & data;

public:

	typedef int yes_i_am_encap;

	typedef T T_type;

	static const int max_prop = T::max_prop;

	// constructor require a key and a memory data
	encapc(type & data)
	:data(data)
	{}

	/*! \brief Return the address of the base
	 *
	 * \return the address of the data encapsulated
	 *
	 */
	type * operator&()
	{
		return &data;
	}

	/*! \brief access the data
	 *
	 * \return the reference
	 *
	 */
	template <unsigned int p> typename type_cpu_prop<p,Mem>::type get()
	{
#ifdef MEMLEAK_CHECK
		check_valid(&boost::fusion::at_c<p>(data),sizeof(typename type_cpu_prop<p,Mem>::type));
#endif
		return boost::fusion::at_c<p>(data);
	}

	/*! \brief access the data but return a constant reference
	 *
	 * \return the reference
	 *
	 */
	template <unsigned int p> const typename type_cpu_prop<p,Mem>::type get() const
	{
#ifdef MEMLEAK_CHECK
		check_valid(&boost::fusion::at_c<p>(data),sizeof(typename type_cpu_prop<p,Mem>::type));
#endif
		return boost::fusion::at_c<p>(data);
	}

	// access the data
	template <unsigned int p> void set(typename type_cpu_prop<p,Mem>::type & ele)
	{
#ifdef MEMLEAK_CHECK
			check_valid(&boost::fusion::at_c<p>(data),sizeof(typename type_cpu_prop<p,T>::type));
#endif
			return boost::fusion::at_c<p>(data) = ele;
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
