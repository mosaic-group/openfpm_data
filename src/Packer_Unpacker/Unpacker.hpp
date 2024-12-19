/*
 * Unpacker.hpp
 *
 *  Created on: Jul 17, 2015
 *      Author: i-bird
 */

#ifndef SRC_UNPACKER_HPP_
#define SRC_UNPACKER_HPP_

#include "util/object_util.hpp"
//#include "Grid/util.hpp"
//#include "Vector/util.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "util/util_debug.hpp"
#include "Pack_selector.hpp"
#include "util/Pack_stat.hpp"
#include "memory/PtrMemory.hpp"
#include "Packer_util.hpp"
#include "util/multi_array_openfpm/multi_array_ref_openfpm.hpp"
#include "has_pack_encap.hpp"
#include "util/object_creator.hpp"
#include "Grid/grid_common.hpp"

/*! \brief Unpacker class
 *
 * \tparam T object type to unpack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 * \tparam Implementation of the unpacker (the Pack_selector choose the correct one)
 *
 */
template<typename T, typename Mem, int pack_type >
class Unpacker
{
public:

	/*! \brief Error, no implementation
	 *
	 */
	static void unpack(ExtPreAlloc<Mem> , T & obj)
	{
#ifndef DISABLE_ALL_RTTI
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " packing for the type " << demangle(typeid(T).name()) << " is not implemented\n";
#endif
	}
};

/*! \brief Unpacker for primitives
 *
 * \tparam T object type to unpack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Unpacker<T,Mem,PACKER_PRIMITIVE>
{
public:


	/*! \brief It unpack C++ primitives
	 *
	 * \param ext preallocated memory from where to unpack the object
	 * \param obj object where to unpack
	 *
	 */
	static void unpack(ExtPreAlloc<Mem> & ext, T & obj,Unpack_stat & ps)
	{
		T * ptr = static_cast<T *>(ext.getPointerOffset(ps.getOffset()));
		obj = *ptr;

		ps.addOffset(sizeof(T));
	}
};

template<typename T, typename Mem>
class Unpacker<T,Mem,PACKER_ARRAY_PRIMITIVE>
{


public:


	/*! \brief It unpack C++ primitives
	 *
	 * \param ext preallocated memory from where to unpack the object
	 * \param obj object where to unpack
	 *
	 */
	static void unpack(ExtPreAlloc<Mem> & ext, T & obj, Unpack_stat & ps)
	{

		//Unpacking a size of a source vector
		size_t u2 = 0;
		Unpacker<size_t, Mem>::unpack(ext,u2,ps);

		//Resize a destination vector
		obj.resize(u2);

		memcpy(obj.getPointer(),ext.getPointerOffset(ps.getOffset()),sizeof(typename T::value_type)*obj.size());

		ps.addOffset(sizeof(typename T::value_type)*obj.size());
	}
};

template <typename T, bool is_array>
struct get_pointer_unpack
{
	static void * getpointer(T & obj)
	{
		return &obj;
	}
};

template <typename T>
struct get_pointer_unpack<T,true>
{
	static void * getpointer(T & obj)
	{
		return obj;
	}
};

template<typename T, typename Mem>
class Unpacker<T,Mem,PACKER_ARRAY_CP_PRIMITIVE>
{
	//! base type of the primitive arrat
	typedef typename std::remove_all_extents<T>::type base_type;

	//! numver of extents in the primitive array
	typedef typename boost::mpl::int_<std::rank<T>::value> n_ext;

	typedef typename to_memory_multi_array_ref_view<T>::vmpl vmpl;

public:

	/*! \brief It packs arrays of C++ primitives
	 *
	 * \param ext preallocated memory where to pack the object
	 * \param obj object to pack
	 * \param sts pack-stat info
	 *
	 */
	inline static void unpack(ExtPreAlloc<Mem> & ext, T & obj, Unpack_stat & ps)
	{
		typedef typename std::remove_extent<T>::type prim_type;

		meta_copy<T>::meta_copy_((prim_type *)ext.getPointerOffset(ps.getOffset()),obj);

		ps.addOffset(sizeof(T));
	}

	/*! \brief It unpack C++ primitives
	 *
	 * \param ext preallocated memory from where to unpack the object
	 * \param obj object where to unpack
	 *
	 */
	static void unpack(ExtPreAlloc<Mem> & ext,
					   typename openfpm::multi_array_ref_openfpm<base_type,n_ext::value+1,vmpl>::reference obj,
					   Unpack_stat & ps)
	{
		meta_copy_d<T,typename openfpm::multi_array_ref_openfpm<base_type,n_ext::value+1,vmpl>::reference>
								::meta_copy_d_(*static_cast<T *>(ext.getPointerOffset(ps.getOffset())),obj);

		ps.addOffset(sizeof(T));
	}
};

/*! \brief Unpacker for objects with no possibility to check for internal pointers
 *
 * \tparam T object type to unpack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Unpacker<T,Mem,PACKER_OBJECTS_WITH_WARNING_POINTERS>
{
public:

	/*! \brief unpack object
	 *
	 * \param ext preallocated memory from where to unpack the object
	 * \param obj object where to unpack
	 * \param ps unpack statistic
	 *
	 */
	static void unpack(ExtPreAlloc<Mem> & ext, T & obj, Unpack_stat & ps)
	{
		memcpy(get_pointer_unpack<T,std::is_array<T>::value>::getpointer(obj),
			   (void *)ext.getPointerOffset(ps.getOffset()),
			   sizeof(T));

		ps.addOffset(sizeof(T));
	}
};

/*! \brief Unpacker class for objects
 *
 * \tparam T object type to unpack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Unpacker<T,Mem,PACKER_OBJECTS_WITH_POINTER_CHECK>
{
public:

	/*! \brief It unpack any object checking that the object does not have pointers inside
	 *
	 * \param ext preallocated memory from where to unpack the object
	 * \param obj object where to unpack
	 *
	 */
	static void unpack(ExtPreAlloc<Mem> & ext, T & obj, Unpack_stat & ps)
	{
		memcpy(&obj,(T *)ext.getPointerOffset(ps.getOffset()),sizeof(T));

		ps.addOffset(sizeof(T));
	}
};

/*! \brief Unpacker for vectors
 *
 * \tparam T object type to unpack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Unpacker<T,Mem,PACKER_GENERAL>
{
public:

	template<unsigned int ... prp> void static unpack(ExtPreAlloc<Mem> & mem, T & obj, Unpack_stat & ps)
	{
		obj.template unpack<prp...>(mem, ps);
	}

	template<unsigned int ... prp> void static unpack(ExtPreAlloc<Mem> & mem, T & obj, Unpack_stat & ps, size_t n)
	{
		if (mem.size() == 0)
			return;

		obj.template unpack<prp...>(mem, ps);
	}
};

/*! \brief Unpacker for grids
 *
 * \tparam T object type to unpack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Unpacker<T,Mem,PACKER_GRID>
{
public:

	template<unsigned int ... prp> static void unpack(ExtPreAlloc<Mem> & mem, T & obj, Unpack_stat & ps)
	{
		obj.template unpack<prp...>(mem, ps);
	}

	template<typename grid_sub_it_type,
			 typename context_type,
			 unsigned int ... prp> static void unpack(ExtPreAlloc<Mem> & mem, grid_sub_it_type & sub_it, T & obj, Unpack_stat & ps, context_type& gpuContext, rem_copy_opt opt)
	{
		obj.template unpack<prp...>(mem, sub_it, ps, gpuContext, opt);
	}

	template<typename grid_sub_it_type,
			 typename header_type,
			 typename context_type,
			 unsigned int ... prp> static void unpack_with_header(ExtPreAlloc<Mem> & mem, 
			 													  grid_sub_it_type & sub_it, 
																  T & obj, 
																  header_type & header, 
																  int ih, 
																  Unpack_stat & ps, 
																  context_type& gpuContext,
																  rem_copy_opt opt)
	{
		obj.template unpack_with_headers<prp...>(mem, sub_it, header, ih, ps, gpuContext, opt);
	}
};


template<unsigned int prop, typename T>
struct std_array_copy_chunks
{
	template<typename Tsrc, typename Tdst>
	static void copy(Tsrc & src, Tdst & dst, size_t pos)
	{
		typedef typename std::remove_reference<decltype(dst.template get<prop>()[0])>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(src.template get<prop>()[0],dst.template get<prop>()[pos]);
	}

	template<typename Tsrc, typename Tdst>
	static void copy_unpack(Tsrc & src, Tdst & dst, size_t pos)
	{
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(src.template get< prop >())>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(src.template get< prop >(),dst.template get< prop >()[pos]);
	}
};

template<unsigned int prop, typename T, unsigned int N1>
struct std_array_copy_chunks<prop,T[N1]>
{
	template<typename Tsrc, typename Tdst>
	static void copy(Tsrc & src, Tdst & dst, size_t pos)
	{
		typedef typename std::remove_reference<decltype(dst.template get<prop>()[0][pos])>::type copy_rtype;

		for (int i = 0 ; i < N1 ; i++)
		{
			meta_copy<copy_rtype>::meta_copy_(dst.template get<prop>()[i][pos],src.template get<prop>()[i][0]);
		}
	}

	template<typename Tsrc, typename Tdst>
	static void copy_unpack(Tsrc & src, Tdst & dst, size_t pos)
	{
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(src.template get< prop >()[0])>::type copy_rtype;

		for (int i = 0 ; i < N1 ; i++)
		{
			meta_copy<copy_rtype>::meta_copy_(src.template get< prop >()[i],dst.template get< prop >()[i][pos]);
		}
	}
};

//! It copy one element of the chunk for each property
template<typename e_src, typename e_dst>
struct copy_unpacker_chunk
{
	//! encapsulated object source
	const e_src & src;
	//! encapsulated object destination
	e_dst & dst;

	//! element to copy
	size_t sub_id;

	/*! \brief constructor
	 *
	 *
	 * \param src source encapsulated object
	 * \param dst destination encapsulated object
	 *
	 */
	inline copy_unpacker_chunk(const e_src & src, e_dst & dst, size_t sub_id)
	:src(src),dst(dst),sub_id(sub_id)
	{
	};



	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename std::remove_reference<decltype(src.template get<T::value>())>::type copy_rtype;

		std_array_copy_chunks<T::value,copy_rtype>::copy_unpack(src,dst,sub_id);
	}
};

/*! \brief Unpacker for encapsulated objects
 *
 * \tparam T object type to unpack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Unpacker<T,Mem,PACKER_ENCAP_OBJECTS_CHUNKING>
{
public:

	/*! \brief
	 *
	 * is this needed
	 *
	 */
	template<typename T_nc, unsigned int ... prp>
	static void unpack(ExtPreAlloc<Mem> & mem,
						T & obj,
						size_t sub_id,
						Unpack_stat & ps)
	{
		if (has_pack_encap<T,prp ...>::result::value == true)
		{call_encapUnpackChunking<T,Mem,prp ...>::call_unpack(obj,sub_id,mem,ps);}
		else
		{
			// get the first element to get the chunking size
			typedef typename boost::mpl::at<typename T::T_type::type,boost::mpl::int_<0>>::type cnk_size;

			if (sizeof...(prp) == 0)
			{
				encapc<1,T_nc,typename memory_traits_lin< T_nc >::type> enc(*static_cast<typename T_nc::type *>((void *)((unsigned char *)mem.getPointerBase() + ps.getOffset())));
				copy_unpacker_chunk<encapc<1,T_nc,typename memory_traits_lin< T_nc >::type>,
									decltype(obj)> cp(enc,obj,sub_id);

				boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::T_type::max_prop>>(cp);
				ps.addOffset(sizeof(T_nc));
			}
			else
			{
				typedef object<typename object_creator_chunking<typename T::type,prp...>::type> prp_object;
				encapc<1,prp_object,typename memory_traits_lin< prp_object >::type> enc(*static_cast<typename prp_object::type *>((void *)((unsigned char *)mem.getPointerBase() + ps.getOffset())));
				object_s_di<decltype(enc),T,OBJ_ENCAP_CHUNKING,prp ... >(enc,obj,sub_id);
				ps.addOffset(sizeof(prp_object));
			}
		}

		// update statistic
	}

	/*! \brief
	 *
	 * is this needed
	 *
	 */
	template<template<typename,typename> class op, unsigned int ... prp>
	static void unpack_op(ExtPreAlloc<Mem> & mem,
						T & obj,
						size_t sub_id,
						Unpack_stat & ps)
	{
		if (has_pack_encap<T,prp ...>::result::value == true)
		{call_encapUnpackChunking<T,Mem,prp ...>::call_unpack(obj,sub_id,mem,ps);}
		else
		{
			if (sizeof...(prp) == 0)
			{
				encapc<1,typename T::T_type,typename memory_traits_lin< typename T::T_type >::type> enc(*static_cast<typename T::T_type::type *>(mem.getPointer()));
				copy_unpacker_chunk<encapc<1,typename T::T_type,typename memory_traits_lin< typename T::T_type >::type>,
									decltype(obj)>(enc,obj,sub_id);
				ps.addOffset(sizeof(typename T::T_type));
			}
			else
			{
				typedef object<typename object_creator_chunking<typename T::type,prp...>::type> prp_object;
				encapc<1,prp_object,typename memory_traits_lin< prp_object >::type> enc(*static_cast<typename prp_object::type *>((void *)((unsigned char *)mem.getPointerBase() + ps.getOffset())));
				object_s_di_op<op,decltype(enc),T,OBJ_ENCAP_CHUNKING,prp ... >(enc,obj,sub_id);
				ps.addOffset(sizeof(prp_object));
			}
		}

		// update statistic
	}
};



#endif /* SRC_UNPACKER_HPP_ */
