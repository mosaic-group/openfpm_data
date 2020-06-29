/*
 * Packer_cls.hpp
 *
 *  Created on: Jul 15, 2015
 *      Author: i-bird
 */

#ifndef SRC_PACKER_HPP_
#define SRC_PACKER_HPP_

#include "util/object_util.hpp"
//#include "Grid/util.hpp"
#include "Vector/util.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "util/util_debug.hpp"


#include "Grid/grid_sm.hpp"
#include "util/Pack_stat.hpp"
#include "Pack_selector.hpp"
#include "has_pack_encap.hpp"
#include "Packer_util.hpp"

template <typename> struct Debug;

/*! \brief Packing class
 *
 * This class pack objects primitives vectors and grids, the general usage is to create a vector of
 * packing request (std::vector<size_t>) that contain the size of the required space needed to pack
 * the information. Calculate the total size, allocating it on HeapMemory (for example), Create an
 * ExtPreAlloc memory object giving the preallocated memory to it and finally Pack all the objects
 * subsequently
 *
 * In order to unpack the information the Unpacker class can be used
 *
 * \see Unpacker
 *
 * \snippet Packer_unit_tests.hpp Pack into a message primitives objects vectors and grids
 *
 * \tparam T object type to pack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 * \tparam Implementation of the packer (the Pack_selector choose the correct one)
 *
 */
template<typename T, typename Mem, int pack_type >
class Packer
{
public:

	/*! \brief Error, no implementation
	 *
	 */
	static void pack(ExtPreAlloc<Mem> , const T & obj)
	{
#ifndef DISABLE_ALL_RTTI
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " packing for the type " << demangle(typeid(T).name()) << " is not implemented\n";
#endif
	}

	/*! \brief Error, no implementation
	 *
	 */
	static size_t packRequest(const T & obj, size_t & req)
	{
#ifndef DISABLE_ALL_RTTI
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " packing for the type " << demangle(typeid(T).name()) << " is not implemented\n";
#endif
		return 0;
	}
};

/*! \brief Packer for primitives
 *
 * \tparam T object type to pack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Packer<T,Mem,PACKER_PRIMITIVE>
{
public:

	/*! \brief It pack any C++ primitives
	 *
	 * \param ext preallocated memory where to pack the object
	 * \param obj object to pack
	 * \param sts pack-stat info
	 *
	 */
	inline static void pack(ExtPreAlloc<Mem> & ext, const T & obj, Pack_stat & sts)
	{
		ext.allocate(sizeof(T));
		*(typename std::remove_const<T>::type *)ext.getPointer() = obj;

		// update statistic
		sts.incReq();
	}

	/*! \brief It add a request to pack a C++ primitive
	 *
	 * \param req requests vector
	 *
	 */
	static void packRequest(const T & obj, size_t & req)
	{
		req += sizeof(T);
	}

	/*! \brief It add a request to pack a C++ primitive
	 *
	 * \param req requests vector
	 *
	 */
	static void packRequest(size_t & req)
	{
		req += sizeof(T);
	}
};

/*! \brief Packer for primitives
 *
 * \tparam T object type to pack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */


template<typename T, typename Mem>
class Packer<T,Mem,PACKER_ARRAY_PRIMITIVE>
{
public:

	/*! \brief It packs arrays of C++ primitives
	 *
	 * \param ext preallocated memory where to pack the object
	 * \param obj object to pack
	 * \param sts pack-stat info
	 *
	 */
	inline static void pack(ExtPreAlloc<Mem> & ext, const T & obj, Pack_stat & sts, size_t n)
	{
		//Pack the size of a vector
		Packer<size_t, Mem>::pack(ext,obj.size(),sts);

		//Pack a vector
		ext.allocate(sizeof(typename T::value_type)*n);
		memcpy(ext.getPointer(),obj.getPointer(),sizeof(typename T::value_type)*n);

		// update statistic
		sts.incReq();
	}

	/*! \brief It add a request to pack a C++ primitive
	 *
	 * \param req requests vector
	 *
	 */
	static void packRequest(T & obj,size_t & req)
	{
		req += sizeof(typename T::value_type)*obj.size();
	}
};


template<typename T, typename Mem>
class Packer<T,Mem,PACKER_ARRAY_CP_PRIMITIVE>
{
public:

	/*! \brief It packs arrays of C++ primitives
	 *
	 * \param ext preallocated memory where to pack the object
	 * \param obj object to pack
	 * \param sts pack-stat info
	 *
	 */
	inline static void pack(ExtPreAlloc<Mem> & ext, const T & obj, Pack_stat & sts)
	{
		typedef typename std::remove_extent<T>::type prim_type;

		//Pack a vector
		ext.allocate(sizeof(T));

		meta_copy<T>::meta_copy_(obj,(prim_type *)ext.getPointer());

		// update statistic
		sts.incReq();
	}

	/*! \brief It packs arrays of C++ primitives
	 *
	 * \param ext preallocated memory where to pack the object
	 * \param obj object to pack
	 * \param sts pack-stat info
	 *
	 */
	template<typename tp, long unsigned int dim, typename vmpl>
	inline static void pack(ExtPreAlloc<Mem> & ext,
		                  	const openfpm::detail::multi_array::sub_array_openfpm<tp,dim,vmpl> & obj,
		                  	Pack_stat & sts)
	{
		typedef typename std::remove_extent<T>::type prim_type;

		//Pack a vector
		ext.allocate(sizeof(T));

		meta_copy<T>::meta_copy_(obj,(prim_type *)ext.getPointer());

		// update statistic
		sts.incReq();
	}

	/*! \brief It add a request to pack a C++ primitive
	 *
	 * \param req requests vector
	 *
	 */
	static void packRequest(T & obj,size_t & req)
	{
		req += sizeof(T);
	}

	/*! \brief It add a request to pack a C++ primitive
	 *
	 * \param req requests vector
	 *
	 */
	template<typename tp, long unsigned int dim, typename vmpl>
	static void packRequest(const openfpm::detail::multi_array::sub_array_openfpm<tp,dim,vmpl> & obj,
			                size_t & req)
	{
		req += sizeof(T);
	}
};

/*! \brief Packer for objects, with impossibility to check for internal pointers
 *
 * \tparam T object type to pack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Packer<T,Mem,PACKER_OBJECTS_WITH_WARNING_POINTERS>
{
public:

	/*! \brief It pack an object
	 *
	 * \param ext preallocated memory where to pack the objects
	 * \param obj object to pack
	 * \param sts pack-stat info
	 *
	 */
	static void pack(ExtPreAlloc<Mem> & ext, const T & obj, Pack_stat & sts)
	{
#ifdef SE_CLASS1
		if (ext.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";

		if (!(std::is_array<T>::value == true && std::is_fundamental<typename std::remove_all_extents<T>::type>::value == true))
			std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " impossible to check the type " << demangle(typeid(T).name()) << " please consider to add a static method like \"static bool noPointers() {return true;}\" \n" ;
#endif
		ext.allocate(sizeof(T));
		memcpy((typename std::remove_reference<T>::type *)ext.getPointer(),&obj,sizeof(T));

		// update statistic
		sts.incReq();
	}

	/*! \brief it add a request to pack an object
	 *
	 * \param req requests vector
	 *
	 */
	static void packRequest(const T & obj,size_t & req)
	{
		req += sizeof(T);
	}

	/*! \brief it add a request to pack an object
	 *
	 * \param req requests vector
	 *
	 */
	static void packRequest(size_t & req)
	{
		req += sizeof(T);
	}
};

/*! \brief Packer class for objects
 *
 * \tparam T object type to pack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Packer<T,Mem,PACKER_OBJECTS_WITH_POINTER_CHECK>
{
public:

	/*! \brief It pack any object checking that the object does not have pointers inside
	 *
	 * \param ext preallocated memory where to pack the objects
	 * \param obj object to pack
	 * \param sts pack-stat info
	 *
	 */
	static void pack(ExtPreAlloc<Mem> & ext, const T & obj, Pack_stat & sts)
	{
#ifdef SE_CLASS1
		if (ext.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";

		if (obj.noPointers() == false)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " the type " << demangle(typeid(T).name()) << " has pointers inside, sending pointers values has no sense\n";
#endif
		ext.allocate(sizeof(T));
		memcpy((typename std::remove_const<T>::type *)ext.getPointer(),&obj,sizeof(T));

		// Update statistic
		sts.incReq();
	}

	/*! \brief it add a request to pack an object
	 *
	 * \param req requests vector
	 *
	 */
	static void packRequest(const T & obj,size_t & req)
	{
		req += sizeof(T);
	}

	/*! \brief it add a request to pack an object
	 *
	 * \param req requests vector
	 *
	 */
	static void packRequest(size_t & req)
	{
		req += sizeof(T);
	}
};

/*! \brief Packer class for vectors
 *
 * \tparam T vector type to pack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Packer<T,Mem,PACKER_GENERAL>
{
public:

	template<int ... prp> static void packRequest(const T & obj, size_t & req)
	{
		obj.template packRequest<prp...>(req);
	}

	template<int ... prp> static void pack(ExtPreAlloc<Mem> & mem, const T & obj, Pack_stat & sts)
	{
		obj.template pack<prp...>(mem, sts);
	}
};

/*! \brief Packer for grids and sub-grids
 *
 * \tparam T grid type to pack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Packer<T,Mem,PACKER_GRID>
{
public:

	template<int ... prp> static void packRequest(const T & obj, size_t & req)
	{
		obj.template packRequest<prp...>(req);
	}

	template<typename grid_sub_it_type, int ... prp> static void packRequest(T & obj, grid_sub_it_type & sub, size_t & req)
	{
		obj.template packRequest<prp...>(sub, req);
	}

	template<int ... prp> static void pack(ExtPreAlloc<Mem> & mem, const T & obj, Pack_stat & sts)
	{
		obj.template pack<prp...>(mem, sts);
	}

	template<typename grid_sub_it_type, int ... prp> static void pack(ExtPreAlloc<Mem> & mem, T & obj, grid_sub_it_type & sub_it, Pack_stat & sts)
	{
		obj.template pack<prp...>(mem, sub_it, sts);
	}
};

template<typename T, typename Mem>
class Packer<T,Mem,PACKER_ENCAP_OBJECTS>
{
public:

	/*! \brief
	 *
	 *
	 */
	template<int ... prp> static void pack(ExtPreAlloc<Mem> & mem, const T & eobj, Pack_stat & sts)
	{
#ifdef SE_CLASS1
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif

		if (has_pack_encap<T,prp ...>::result::value == true)
		{call_encapPack<T,Mem,prp ...>::call_pack(eobj,mem,sts);}
		else
		{
			if (sizeof...(prp) == 0)
			{
				mem.allocate(sizeof(typename T::T_type));
				encapc<1,typename T::T_type,typename memory_traits_lin< typename T::T_type >::type> enc(*static_cast<typename T::T_type::type *>(mem.getPointer()));
				enc = eobj;
			}
			else
			{
				typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;
				mem.allocate(sizeof(prp_object));
				encapc<1,prp_object,typename memory_traits_lin< prp_object >::type> enc(*static_cast<typename prp_object::type *>(mem.getPointer()));
				object_si_d<T,decltype(enc),OBJ_ENCAP,prp ... >(eobj,enc);
			}
		}

		// update statistic
		sts.incReq();
	}

	/*! \brief
	 *
	 *
	 */
	template<int ... prp> void packRequest(T & eobj,size_t & req)
	{
		if (has_pack_encap<T>::value == true)
		{call_encapPackRequest<T,Mem,prp ...>::call_packRequest(eobj,req);}
		else
		{
			if (sizeof...(prp) == 0)
				return;

			typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

			req += sizeof(prp_object);
		}
	}
};

///////////////////////////////////

//! It copy one element of the chunk for each property
template<typename e_src, typename e_dst>
struct copy_packer_chunk
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
	inline copy_packer_chunk(const e_src & src, size_t sub_id, e_dst & dst)
	:src(src),dst(dst),sub_id(sub_id)
	{
	};



	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<decltype(dst.template get< T::value >())>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(src.template get< T::value >()[sub_id],dst.template get< T::value >());
	}
};

template <typename> struct Debug;

template<bool is_not_zero, typename T, typename T_nc, int ... prp>
struct selector_chunking_prp_has_zero_size
{
	template<typename Mem> static inline void select(ExtPreAlloc<Mem> & mem,  const T & eobj, size_t sub_id)
	{
		mem.allocate(sizeof(T_nc));
		encapc<1,T_nc,typename memory_traits_lin< T_nc >::type> enc(*static_cast<typename T_nc::type *>(mem.getPointer()));
		copy_packer_chunk<decltype(eobj),
		                  encapc<1,T_nc,typename memory_traits_lin< T_nc >::type>> cp(eobj,sub_id,enc);
		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::T_type::max_prop>>(cp);
		//enc = eobj[sub_id];
	}
};

template<typename T,  typename T_nc, int ... prp>
struct selector_chunking_prp_has_zero_size<true,T,T_nc,prp...>
{
	template<typename Mem> static inline void select(ExtPreAlloc<Mem> & mem,  const T & eobj, size_t sub_id)
	{
		// Here we create an encap without array extent
		typedef object<typename object_creator_chunking<typename T::type,prp...>::type> prp_object;
		mem.allocate(sizeof(prp_object));
		encapc<1,prp_object,typename memory_traits_lin< prp_object>::type> enc(*static_cast<typename prp_object::type *>(mem.getPointer()));
		object_si_d<T,decltype(enc),OBJ_ENCAP_CHUNKING,prp ... >(eobj,sub_id,enc);
	}
};

template<typename T, typename Mem>
class Packer<T,Mem,PACKER_ENCAP_OBJECTS_CHUNKING>
{
public:

	/*! \brief Pack the encapsulkated object
	 *
	 * \tparam Base object without chunking
	 *
	 */
	template<typename T_nc, int ... prp> static void pack(ExtPreAlloc<Mem> & mem, const T & eobj, size_t sub_id, Pack_stat & sts)
	{
#ifdef SE_CLASS1
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif

		if (has_pack_encap<T,prp ...>::result::value == true)
		{call_encapPackChunking<T,Mem,prp ...>::call_pack(eobj,sub_id,mem,sts);}
		else
		{
			selector_chunking_prp_has_zero_size<sizeof...(prp) != 0,T, T_nc, prp ...>::select(mem,eobj,sub_id);
		}

		// update statistic
		sts.incReq();
	}

	/*! \brief
	 *
	 *
	 */
	template<int ... prp> void packRequest(T & eobj,size_t & req)
	{
		if (has_pack_encap<T>::value == true)
			call_encapPackRequest<T,Mem,prp ...>::call_packRequest(eobj,req);
		else
		{
			if (sizeof...(prp) == 0)
				return;

			typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

			req += sizeof(prp_object);
		}
	}
};

#endif /* SRC_PACKER_HPP_ */
