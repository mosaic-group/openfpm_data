/*
 * grid_common.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_COMMON_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_COMMON_HPP_

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to call hostToDevice for each properties
 *
 */
template<typename Tv>
struct host_to_dev_all_prp
{
	Tv & p;

	inline host_to_dev_all_prp(Tv & p)
	:p(p)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		p.template hostToDevice<T::value>();
	}
};

template<typename T, bool is_vector>
struct call_recursive_host_device_if_vector
{
	//! It is a vector recursively call deviceToHost
	template<typename obj_type>
	static void call(obj_type & obj, size_t start, size_t stop)
	{
		T * ptr = static_cast<T *>(obj.get_pointer());

		for(size_t i = start ; i < stop ; i++)
		{
			host_to_dev_all_prp<T> hdap(ptr[i]);

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::value_type::max_prop>>(hdap);
		}
	}
};

template<typename T>
struct call_recursive_host_device_if_vector<T,false>
{
	//! It is not a vector nothing to do
	template<typename obj_type>
	static void call(obj_type & obj, size_t start, size_t stop) {}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one object into one target
 * grid  element in a generic way for a
 * generic object T with variable number of property
 *
 * \tparam dim Dimensionality
 * \tparam S type of grid
 * \tparam Memory type of memory needed for encap
 *
 */

template<unsigned int dim, typename S, typename Memory>
struct copy_cpu_encap
{
	//! size to allocate
	grid_key_dx<dim> & key;

	//! grid where we have to store the data
	S & grid_dst;

	//! type of the object we have to set
	typedef typename S::value_type obj_type;

	//! type of the object boost::sequence
	typedef typename S::value_type::type ov_seq;

	//! object we have to store
	const encapc<1,obj_type,Memory> & obj;

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param key which element we are modifying
	 * \param grid_dst grid we are updating
	 * \param obj object we have to set in grid_dst (encapsulated)
	 *
	 */
	inline copy_cpu_encap(grid_key_dx<dim> & key, S & grid_dst, const encapc<1,obj_type,Memory> & obj)
	:key(key),grid_dst(grid_dst),obj(obj){};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_cpu_encap(grid_key_dx<dim> & key, S & grid_dst, const encapc<1,obj_type,Memory> && obj)
	:key(key),grid_dst(grid_dst),obj(obj)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(grid_dst.template get<T::value>(key))>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(obj.template get<T::value>(),grid_dst.template get<T::value>(key));
	}
};


/*! \brief Metafunction take T and return a reference
 *
 * Metafunction take T and return a reference
 *
 * \param T type
 *
 */

template<typename T>
struct mem_reference
{
	typedef T& type;
};



#endif /* OPENFPM_DATA_SRC_GRID_GRID_COMMON_HPP_ */
