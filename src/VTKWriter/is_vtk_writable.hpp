/*
 * is_vtk_writable.hpp
 *
 *  Created on: Jul 18, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_VTKWRITER_IS_VTK_WRITABLE_HPP_
#define OPENFPM_IO_SRC_VTKWRITER_IS_VTK_WRITABLE_HPP_

template<typename T, bool is_w>
struct vtk_type
{
	typedef decltype(std::declval<T>().get_vtk(0)) type;
};

template<typename T>
struct vtk_type<T,false>
{
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

template<typename ObjType, bool has_dims = is_vtk_vector_dims<ObjType>::value >
struct vtk_dims
{
	enum
	{
		value = 1
	};
};


template<typename ObjType >
struct vtk_dims<ObjType,true>
{
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


template<typename T>
struct is_vtk_writable
{
	enum
	{
		value = is_custom_vtk_writable<T>::value
	};
};

template<>
struct is_vtk_writable<float>
{
	enum
	{
		value = true
	};
};

template<>
struct is_vtk_writable<double>
{
	enum
	{
		value = true
	};
};

template<>
struct is_vtk_writable<char>
{
	enum
	{
		value = true
	};
};

template<>
struct is_vtk_writable<unsigned char>
{
	enum
	{
		value = true
	};
};

template<>
struct is_vtk_writable<short>
{
	enum
	{
		value = true
	};
};

template<>
struct is_vtk_writable<unsigned short>
{
	enum
	{
		value = true
	};
};


template<>
struct is_vtk_writable<int>
{
	enum
	{
		value = true
	};
};

template<>
struct is_vtk_writable<unsigned int>
{
	enum
	{
		value = true
	};
};

template<>
struct is_vtk_writable<long int>
{
	enum
	{
		value = true
	};
};

template<>
struct is_vtk_writable<unsigned long int>
{
	enum
	{
		value = true
	};
};


template<>
struct is_vtk_writable<bool>
{
	enum
	{
		value = true
	};
};

#endif /* OPENFPM_IO_SRC_VTKWRITER_IS_VTK_WRITABLE_HPP_ */
