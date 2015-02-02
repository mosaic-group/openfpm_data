/*
 * Encap.hpp
 *
 *  Created on: Feb 2, 2015
 *      Author: i-bird
 */

#ifndef ENCAP_HPP_
#define ENCAP_HPP_

#include "grid.hpp"

/*! \brief This class is an helper to get the return type for get method for each property
 *
 * This class is an helper to get the return type for get method for each property
 *
 * \param p id of the property
 * \param T original boost fusion vector, T is suppose to be a boost::fusion::vector<memory_c<...>,.....>
 *
 */

template<unsigned int p,typename T>
struct type_cpu_prop
{
	//! return a boost::fusion::vector<memory_c<....>....>
	typedef typename T::memory_lin::vtype vtype;
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

template<unsigned int p,typename T>
struct type_gpu_prop
{
	//! return a boost::fusion::vector<memory_c<....>....>
	typedef typename T::memory_int vtype;
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
	typedef typename T::type type;

	type & data;

public:

	typedef T T_type;

	// constructor require a key and a memory data
	encapc(type & data)
	:data(data)
	{}

	/*! \brief access the data
	 *
	 * \return the reference
	 *
	 */
	template <unsigned int p> typename type_cpu_prop<p,T>::type get()
	{
#ifdef MEMLEAK_CHECK
		check_valid(&boost::fusion::at_c<p>(data),sizeof(typename type_cpu_prop<p,T>::type));
#endif
		return boost::fusion::at_c<p>(data);
	}

	/*! \brief access the data but return a constant reference
	 *
	 * \return the reference
	 *
	 */
	template <unsigned int p> const typename type_cpu_prop<p,T>::type get() const
	{
#ifdef MEMLEAK_CHECK
		check_valid(&boost::fusion::at_c<p>(data),sizeof(typename type_cpu_prop<p,T>::type));
#endif
		return boost::fusion::at_c<p>(data);
	}

	/*! \brief Get the data
	 *
	 * \return the data
	 *
	 */
	template <unsigned int p> typename boost::remove_reference<typename type_cpu_prop<p,T>::type>::type get() const
	{
#ifdef MEMLEAK_CHECK
		check_valid(&boost::fusion::at_c<p>(data),sizeof(typename type_cpu_prop<p,T>::type));
#endif
		return boost::fusion::at_c<p>(data);
	}

	// access the data
	template <unsigned int p> void set(typename type_cpu_prop<p,T>::type & ele)
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
	grid_key_dx<dim> & k;
	grid<dim,void> & g1;

public:

	typedef T T_type;

	// constructor require a key and a memory data
	encapg(Mem & data, grid_key_dx<dim> & k, grid<dim,void> & g1)
	:data(data),k(k),g1(g1)
	{}

	// access the data
	template <unsigned int p> typename type_gpu_prop<p,T>::type::reference get()
	{
#ifdef MEMLEAK_CHECK
		check_valid(&boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(k))));
#endif
		return boost::fusion::at_c<p>(data).mem_r->operator[](g1.LinId(k));
	}
};


#endif /* ENCAP_HPP_ */
