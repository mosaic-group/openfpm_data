/*
 * is_vtk_writable.hpp
 *
 *  Created on: Jul 18, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_VTKWRITER_IS_VTK_WRITABLE_HPP_
#define OPENFPM_IO_SRC_VTKWRITER_IS_VTK_WRITABLE_HPP_

//! the vtk type
template<typename T, bool is_w>
struct vtk_type
{
	//! get the vtk type for the property
	typedef decltype(std::declval<T>().get_vtk(0)) type;
};

//! the vtk type
template<typename T>
struct vtk_type<T,false>
{
	//! non writable vtk property (so void)
	typedef void type;
};

/*! \brief it check if the type is vtk writable
 *
 *
 */
template<typename ObjType, typename Sfinae = void>
struct is_custom_vtk_writable: std::false_type {};

/*! \brief it check if the type is vtk writable
 *
 * \tparam ObjType check the type
 *
 */
template<typename ObjType>
struct is_custom_vtk_writable <ObjType, typename Void< typename ObjType::is_vtk_writable >::type> : std::true_type
{};

/*! \brief it check if the type is vtk writable
 *
 *
 */
template<typename ObjType, typename Sfinae = void>
struct is_vtk_vector_dims: std::false_type {};

//! If it has not dims property defined the object is considered scalar
template<typename ObjType, bool has_dims = is_vtk_vector_dims<ObjType>::value >
struct vtk_dims
{
	//! dimensionality of the vtk property (scalar)
	enum
	{
		value = 1
	};
};

//! return the dimansionality of the object
template<typename ObjType >
struct vtk_dims<ObjType,true>
{
	//! dimansionality of the vtk proverty (vector) in case of an object point
	//! or an object that define dimansionality
	enum
	{
		value = ObjType::dims
	};
};

/*! \brief it check if the type is vtk writable
 *
 * \tparam ObjType check the type
 *
 */
template<typename ObjType>
struct is_vtk_vector_dims<ObjType, typename Void< decltype(ObjType::dims) >::type> : std::true_type
{};

//! check for T to be writable
template<typename T>
struct is_vtk_writable
{
	//! It check if the object is vtk compatible
	enum
	{
		value = is_custom_vtk_writable<T>::value
	};
};

//! check float
template<>
struct is_vtk_writable<float>
{
	//! float is vtk writable
	enum
	{
		value = true
	};
};

//! check double
template<>
struct is_vtk_writable<double>
{
	//! double is vtk writable
	enum
	{
		value = true
	};
};

//! check char
template<>
struct is_vtk_writable<char>
{
	//! char is vtk writable
	enum
	{
		value = true
	};
};

//! check unsigned char
template<>
struct is_vtk_writable<unsigned char>
{
	//! unsigned char is vtk writable
	enum
	{
		value = true
	};
};

//! check short
template<>
struct is_vtk_writable<short>
{
	//! short is vtk writable
	enum
	{
		value = true
	};
};

//! check unsigned short
template<>
struct is_vtk_writable<unsigned short>
{
	//! unsigned short is vtk writable
	enum
	{
		value = true
	};
};

//! check int
template<>
struct is_vtk_writable<int>
{
	//! int is vtk writable
	enum
	{
		value = true
	};
};

//! check unsigned int
template<>
struct is_vtk_writable<unsigned int>
{
	//! unsigned int is vtk writable
	enum
	{
		value = true
	};
};

//! check long int
template<>
struct is_vtk_writable<long int>
{
	//! long int is vtk writable
	enum
	{
		value = false
	};
};

//! check unsigned long int
template<>
struct is_vtk_writable<unsigned long int>
{
	//! unsigned long int is vtk writable
	enum
	{
		value = false
	};
};

//! check bool
template<>
struct is_vtk_writable<bool>
{
	//! bool is vtk writable
	enum
	{
		value = true
	};
};

#endif /* OPENFPM_IO_SRC_VTKWRITER_IS_VTK_WRITABLE_HPP_ */
