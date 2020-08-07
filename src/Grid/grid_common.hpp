/*
 * grid_common.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_COMMON_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_COMMON_HPP_

#include <type_traits>
#include "util/tokernel_transformation.hpp"

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to call hostToDevice for each properties
 *
 */
template<typename aggrT_src, typename local_grids_type>
struct setBackground_impl
{
	aggrT_src & bck;

	local_grids_type & loc_grid;

	inline setBackground_impl(aggrT_src & bck, local_grids_type & loc_grid)
	:bck(bck),loc_grid(loc_grid)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{loc_grid.get(i).template setBackgroundValue<T::value>(bck.template get<T::value>());}
	}
};

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


template<typename T, typename T_ker, typename type_prp, template<typename> class layout_base , int is_vector>
struct call_recursive_host_device_if_vector
{
	template<typename mem_type, typename obj_type> static void transform(mem_type * mem, obj_type & obj, size_t start, size_t stop)
	{
		start /= sizeof(type_prp);
		stop /= sizeof(type_prp);

		// The type of device and the type on host does not match (in general)
		// So we have to convert before transfer

		T * ptr = static_cast<T *>(obj.get_pointer());

		mem_type tmp;

		tmp.allocate(mem->size());

		T_ker * ptr_tt = static_cast<T_ker *>(tmp.getPointer());

		for(size_t i = start ; i < stop ; i++)
		{
			new (&ptr_tt[i]) T_ker();
			ptr_tt[i] = ptr[i].toKernel();
		}

		mem->hostToDevice(tmp);
	}

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

template<typename T, typename T_ker, typename type_prp ,template<typename> class layout_base>
struct call_recursive_host_device_if_vector<T,T_ker,type_prp,layout_base,0>
{
	template<typename mem_type,typename obj_type> static void transform(mem_type * mem, obj_type & obj, size_t start, size_t stop)
	{
		mem->hostToDevice(start,stop);
	}

	//! It is not a vector nothing to do
	template<typename obj_type>
	static void call(obj_type & obj, size_t start, size_t stop) {}
};

template<typename T, typename T_ker, typename type_prp ,template<typename> class layout_base>
struct call_recursive_host_device_if_vector<T,T_ker,type_prp,layout_base,3>
{
	template<typename mem_type,typename obj_type> static void transform(mem_type * mem, obj_type & obj, size_t start, size_t stop)
	{
		// calculate the start and stop elements
		start /= std::extent<type_prp,0>::value;
		stop /= std::extent<type_prp,0>::value;
		size_t sz = mem->size() / std::extent<type_prp,0>::value;

		size_t offset = 0;
		for (size_t i = 0 ; i < std::extent<type_prp,0>::value ; i++)
		{
			mem->hostToDevice(offset+start,offset+stop);
			offset += sz;
		}
	}

	//! It is not a vector nothing to do
	template<typename obj_type>
	static void call(obj_type & obj, size_t start, size_t stop) {}
};

template<typename T, typename T_ker, typename type_prp ,template<typename> class layout_base>
struct call_recursive_host_device_if_vector<T,T_ker,type_prp,layout_base,4>
{
	template<typename mem_type,typename obj_type> static void transform(mem_type * mem, obj_type & obj, size_t start, size_t stop)
	{
		// calculate the start and stop elements
		start = start / std::extent<type_prp,0>::value / std::extent<type_prp,1>::value;
		stop = stop / std::extent<type_prp,0>::value / std::extent<type_prp,1>::value;
		size_t sz = mem->size() / std::extent<type_prp,0>::value / std::extent<type_prp,1>::value;

		size_t offset = 0;
		for (size_t i = 0 ; i < std::extent<type_prp,0>::value ; i++)
		{
			for (size_t j = 0 ; j < std::extent<type_prp,1>::value ; j++)
			{
				mem->hostToDevice(offset+start,offset+stop);
				offset += sz;
			}
		}
	}

	//! It is not a vector nothing to do
	template<typename obj_type>
	static void call(obj_type & obj, size_t start, size_t stop) {}
};

/////////// destructor


template<typename T, typename T_ker, typename type_prp, template<typename> class layout_base , int is_vector>
struct call_recursive_destructor_if_vector
{
	template<typename mem_type, typename obj_type> static void destruct(mem_type * mem, obj_type & obj)
	{
		size_t sz = mem->size() / sizeof(type_prp);
		// The type of device and the type on host does not match (in general)
		// So we have to convert before transfer

		mem_type tmp;

		tmp.allocate(mem->size());

		mem->deviceToHost(tmp);
		T_ker * ptr = static_cast<T_ker *>(tmp.getPointer());

		for(size_t i = 0 ; i < sz ; i++)
		{
			ptr->~T_ker();
			++ptr;
		}
	}
};

template<typename T, typename T_ker, typename type_prp ,template<typename> class layout_base>
struct call_recursive_destructor_if_vector<T,T_ker,type_prp,layout_base,0>
{
	template<typename mem_type,typename obj_type> static void destruct(mem_type * mem, obj_type & obj)
	{}
};

template<typename T, typename T_ker, typename type_prp ,template<typename> class layout_base>
struct call_recursive_destructor_if_vector<T,T_ker,type_prp,layout_base,3>
{
	template<typename mem_type,typename obj_type> static void destruct(mem_type * mem, obj_type & obj)
	{}
};

template<typename T, typename T_ker, typename type_prp ,template<typename> class layout_base>
struct call_recursive_destructor_if_vector<T,T_ker,type_prp,layout_base,4>
{
	template<typename mem_type,typename obj_type> static void destruct(mem_type * mem, obj_type & obj)
	{}
};

///////////////////////

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

/*! \brief Options for remove copy
 *
 *
 */
enum rem_copy_opt
{
	NONE_OPT = 0,
	PHASE1 = 0,
	PHASE2 = 1,
	PHASE3 = 2,
	// This option indicate that the geometrical structure of the sparse-grid has not changed
	KEEP_GEOMETRY = 4,
};

#endif /* OPENFPM_DATA_SRC_GRID_GRID_COMMON_HPP_ */
