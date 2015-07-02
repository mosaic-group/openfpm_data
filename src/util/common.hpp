#ifndef COMMON_HPP
#define COMMON_HPP

#include <type_traits>
#include <random>
#include "memory/memory.hpp"

namespace std
{
	// We need the definition of std::to_string that work on string

	static std::string to_string(std::string s)
	{
		return s;
	}
}


//! Void structure
template<typename> struct Void
{
	//! define void type
	typedef void type;
};

template<typename T, typename Sfinae = void>
struct has_attributes: std::false_type {};


/*! \brief has_attributes check if a type has defined an
 * internal structure with attributes
 *
 * ### Example
 *
 * \code{.cpp}
 * has_attributes<Test>::value
 * \endcode
 *
 * return true if T::attributes::name[0] is a valid expression
 * and produce a defined type
 *
 */
template<typename T>
struct has_attributes<T, typename Void<decltype( T::attributes::name[0] )>::type> : std::true_type
{};

template<typename T, typename Sfinae = void>
struct has_typedef_type: std::false_type {};

/*! \brief has_typedef_type check if a typedef ... type inside the structure is
 *         defined
 *
 * ### Example
 *
 * \code{.cpp}
 * has_typedef_type<Test>::value
 * \endcode
 *
 * return true if T::type is a valid type
 *
 */
template<typename T>
struct has_typedef_type<T, typename Void< typename T::type>::type> : std::true_type
{};

/*! \brief has_data check if a type has defined a member data
 *
 * ### Example
 *
 * \code{.cpp}
 * has_data<Test>::value
 * \endcode
 *
 * return true if T::type is a valid type
 *
 */

template<typename T, typename Sfinae = void>
struct has_data: std::false_type {};

template<typename T>
struct has_data<T, typename Void<decltype( T::data )>::type> : std::true_type
{};


/*! \brief check if T::type and T.data has the same type
 *
 * \tparam i when different from 0 a check is performed otherwise not, the reason of
 *           this is that the typedef and data could also not exist producing
 *           compilation error, this flag avoid this, it perform the check only if it
 *           is safe
 *
 * \tparam T
 *
 * ### Example
 *
 * \code{.cpp}
 * is_typedef_and_data_same<has_data<T> && has_typedef<T>,T>::value
 * \endcode
 *
 * return true if the type of T::data is the same of T::type, false otherwise
 *
 */
template<bool cond, typename T>
struct is_typedef_and_data_same
{
	enum
	{
		value = std::is_same<decltype(T().data),typename T::type>::value
	};
};


template<typename T>
struct is_typedef_and_data_same<false,T>
{
	enum
	{
		value = false
	};
};


#endif
