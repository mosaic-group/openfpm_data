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

/*! \brief Unpacker class
 *
 * \tparam T object type to unpack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 * \tparam Implementation of the unpacker (the Pack_selector choose the correct one)
 *
 */
template<typename T, typename Mem, int pack_type=Pack_selector<T>::value >
class Unpacker
{
public:

	/*! \brief Error, no implementation
	 *
	 */
	static void unpack(ExtPreAlloc<Mem> , T & obj)
	{
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " packing for the type " << demangle(typeid(T).name()) << " is not implemented\n";
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
	 *
	 */
	static void unpack(ExtPreAlloc<Mem> & ext, T & obj, Unpack_stat & ps)
	{
#ifdef DEBUG
		std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " impossible to check the type " << demangle(typeid(T).name()) << " please consider to add a static method like \"static bool noPointers() {return true;}\" \n" ;
#endif
		memcpy(&obj,(T *)ext.getPointerOffset(ps.getOffset()),sizeof(T));

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
#ifdef DEBUG
		if (obj.noPointers() == false)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " the type " << demangle(typeid(T).name()) << " has pointers inside, sending pointers values has no sense\n";
#endif
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
class Unpacker<T,Mem,PACKER_VECTOR>
{
public:

	template<unsigned int ... prp> void static unpack(ExtPreAlloc<Mem> & mem, T & obj, Unpack_stat & ps)
	{
#ifdef DEBUG
		std::cout << "Inside vector unpack() function! (unpacker.hpp)" << std::endl;
#endif
		obj.template unpack<prp...>(mem, ps);
	};
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
	};

	template<unsigned int ... prp> static void unpack(ExtPreAlloc<Mem> & mem, grid_key_dx_iterator_sub<T::dims> & sub_it, T & obj, Unpack_stat & ps)
	{
		obj.template unpack<prp...>(mem, sub_it, ps);
	};
};

/*! \brief Unpacker for encapsulated objects
 *
 * \tparam T object type to unpack
 * \tparam Mem Memory origin HeapMemory CudaMemory ...
 *
 */
template<typename T, typename Mem>
class Unpacker<T,Mem,PACKER_ENCAP_OBJECTS>
{
public:

	// TODO

	/*! \brief
	 *
	 *
	 */
/*	void pack(ExtPreAlloc<Mem> & mem, T & eobj)
	{
		// Create an object out of the encapsulated object and copy
		typename T::type obj = eobj;

		memcpy(mem.getPointer(),&obj,sizeof(T::type));
	}*/

	/*! \brief
	 *
	 *
	 */
/*	void packRequest(std::vector<size_t> & v)
	{
		v.push_back(sizeof(T::type));
	}*/
};


#endif /* SRC_UNPACKER_HPP_ */
