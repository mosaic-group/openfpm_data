/*
 * Pack_selector.hpp
 *
 *  Created on: Jul 15, 2015
 *      Author: Pietro Incardona
 */


#ifndef SRC_PACK_AAA_SELECTOR_HPP_
#define SRC_PACK_AAA_SELECTOR_HPP_

#include <type_traits>
#include "util/common.hpp"
#include "Grid/Encap.hpp"
#include "Grid/util.hpp"
#include "Vector/util.hpp"
#include "Packer_Unpacker/has_pack_agg.hpp"

//! Primitive packing
#define PACKER_PRIMITIVE 1
//! Array of primitives packing
#define PACKER_ARRAY_PRIMITIVE 2
//! Encapsulated Object packing
#define PACKER_ENCAP_OBJECTS 3
//! Vector of objects packing
#define PACKER_GENERAL 4
//! Grid packing
#define PACKER_GRID 5
//! Packer cannot check for pointers
#define PACKER_OBJECTS_WITH_WARNING_POINTERS 6
//! Packer error structure has pointers
#define PACKER_OBJECTS_WITH_POINTER_CHECK 7
//! Compile time array packing
#define PACKER_ARRAY_CP_PRIMITIVE 8

#define IS_ENCAP 4
#define IS_GRID 2
#define HAS_PACKER 1

/*! \brief Pack selector for unknown type
 *
 *
 */
template <typename T, bool has_noPointers>
struct Pack_selector_unknown_type_impl
{
	enum
	{
		value = PACKER_OBJECTS_WITH_POINTER_CHECK
	};
};

template <typename T>
struct Pack_selector_unknown_type_impl<T,false>
{
	enum
	{
		value = PACKER_OBJECTS_WITH_WARNING_POINTERS
	};
};

/*! \brief Pack selector for unknown type
 *
 *
 */
template <typename T, int known_type>
struct Pack_selector_known_type_impl
{
	enum
	{
		value = Pack_selector_unknown_type_impl<T, has_noPointers<T>::value >::value
	};
};

template <typename T>
struct Pack_selector_known_type_impl<T,IS_GRID>
{
	enum
	{
		value = PACKER_GRID
	};
};

template <typename T>
struct Pack_selector_known_type_impl<T,HAS_PACKER>
{
	enum
	{
		value = PACKER_GENERAL
	};
};

template <typename T>
struct Pack_selector_known_type_impl<T,IS_ENCAP>
{
	enum
	{
		value = PACKER_ENCAP_OBJECTS
	};
};

/////////////////////// ---------- CHECKING FOR PRIMITIVES --------------
/*! \brief it is not a fundamental type
 *
 */
template<typename T, bool is_foundamental>
struct Pack_selector_impl
{
	enum
	{
		value = Pack_selector_known_type_impl< T, 4*is_encap<T>::value + is_grid<T>::value * 2 + (has_pack_gen<T>::value)*(!is_grid<T>::value) >::value
	};
};

template<typename T, bool impl = std::is_array<T>::value>
struct Pack_selector_array_prim
{
	enum
	{
		value = PACKER_PRIMITIVE
	};
};

template<typename T>
struct Pack_selector_array_prim<T,true>
{
	enum
	{
		value = PACKER_ARRAY_CP_PRIMITIVE
	};
};

/*! \brief Select the primitive packing
 *
 */
template<typename T>
struct Pack_selector_impl<T,true>
{
	enum
	{
		value = Pack_selector_array_prim<T>::value
	};
};

//////////////////////////////////////////////////////////////////////////

/*! \brief Pack selector
 *
 *
 */
template <typename T>
struct Pack_selector
{
	enum
	{
		value = Pack_selector_impl< T,std::is_fundamental<T>::value  || std::is_fundamental< typename std::remove_all_extents<T>::type >::value >::value
	};
};


#endif /* SRC_PACK_SELECTOR_HPP_ */
