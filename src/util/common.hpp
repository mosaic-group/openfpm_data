#ifndef COMMON_HPP
#define COMMON_HPP

#include <type_traits>
#include <random>
#include "memory/memory.hpp"
#include "util/for_each_ref.hpp"
#include "util/variadic_to_vmpl.hpp"

#define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)

namespace std
{
	// We need the definition of std::to_string that work on string
	inline static std::string to_string(std::string s)
	{
		return s;
	}
}

/*! \brief convert a type into constant type
 *
 * \param T type tp convert
 *
 * \return a constant of the type T
 *
 */
namespace openfpm
{
        template <class T>
        constexpr typename std::add_const<T>::type & as_const(T& t) noexcept
        {
                return t;
        }

        template <typename> struct Debug;
}

 //! Compile time array functor needed to generate array at compile-time of type
 // {3,3,3,3,3,3,.....}
 template<size_t index, size_t N> struct Fill_three {
    enum { value = 3 };
 };

 //! {0,0,0,0,....}
 template<size_t index, size_t N> struct Fill_zero {
    enum { value = 0 };
 };

 //! {2,2,2,2,....}
 template<size_t index, size_t N> struct Fill_two {
    enum { value = 2 };
 };

 //! {1,1,1,1,....}
 template<size_t index, size_t N> struct Fill_one {
    enum { value = 1 };
 };

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
 * \snippet util_test.hpp Declaration of struct with attributes and without
 * \snippet util_test.hpp Check has_attributes
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
 * \snippet util_test.hpp Check has_typedef_type
 *
 * return true if T::type is a valid type
 *
 */
template<typename T>
struct has_typedef_type<T, typename Void< typename T::type>::type> : std::true_type
{};

template<typename T, typename Sfinae = void>
struct has_vector_kernel: std::false_type {};

/*! \brief has_vector_kernel check if a type has defined a member data
 *
 * ### Example
 *
 * \snippet util_test.hpp Check has_data
 *
 * return true if T::vector_kernel is a valid type
 *
 */
template<typename T>
struct has_vector_kernel<T, typename Void< typename T::vector_kernel >::type> : std::true_type
{};

template<typename T, typename Sfinae = void>
struct has_data: std::false_type {};

/*! \brief has_data check if a type has defined a member data
 *
 * ### Example
 *
 * \snippet util_test.hpp Check has_data
 *
 * return true if T::type is a valid type
 *
 */
template<typename T>
struct has_data<T, typename Void<decltype( T::data )>::type> : std::true_type
{};

template<typename T, typename Sfinae = void>
struct has_posMask: std::false_type {};

/*! \brief has_posMask check if a type has defined a member stag_mask
 *
 * It is used to indicate a staggered grid
 *
 * return true if T::stag_mask is a valid type
 *
 */
template<typename T>
struct has_posMask<T, typename Void<decltype( T::stag_mask )>::type> : std::true_type
{};

template<typename T, typename Sfinae = void>
struct has_check_device_pointer: std::false_type {};

/*! \brief has_check_device_pointer check if a type has defined a member yes_has_check_device_pointer
 *
 * This mean that the class support a way to check if it is the owner od a particular device pointer
 *
 *
 * return true if T::yes_has_check_device_pointer is a valid type
 *
 */
template<typename T>
struct has_check_device_pointer<T, typename Void< typename T::yes_has_check_device_pointer >::type> : std::true_type
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
 * \snippet util_test.hpp Check is_typedef_and_data_same
 *
 * return true if the type of T::data is the same of T::type, false otherwise
 *
 */
template<bool cond, typename T>
struct is_typedef_and_data_same
{
	enum
	{
		value = std::is_same<decltype(std::declval<T>().data),typename T::type>::value
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


template<typename T, typename Sfinae = void>
struct has_noPointers: std::false_type {};


/*! \brief has_noPointers check if a type has defined a
 * method called noPointers
 *
 * ### Example
 *
 * \snippet util_test.hpp Check no pointers
 *
 * return true if T::noPointers() is a valid expression (function pointers)
 * and produce a defined type
 *
 */
template<typename T>
struct has_noPointers<T, typename Void<decltype( T::noPointers() )>::type> : std::true_type
{};

/*! \brief has_Pack check if a type has defined a
 * method called Pack
 *
 * ### Example
 *
 * \snippet
 *
 * return true if T::pack() is a valid expression (function pointers)
 * and produce a defined type
 *
 */

template<typename ObjType, typename Sfinae = void>
struct has_pack: std::false_type {};

template<typename ObjType>
struct has_pack<ObjType, typename Void<decltype( ObjType::pack() )>::type> : std::true_type
{};

/*! \brief has_toKernel check if a type has defined a
 * method called toKernel
 *
 * ### Example
 *
 * \snippet
 *
 * return true if T.toKernel() is a valid expression
 * and produce a defined type
 *
 */

template<typename ObjType, typename Sfinae = void>
struct has_toKernel: std::false_type {};

template<typename ObjType>
struct has_toKernel<ObjType, typename Void<decltype( std::declval<ObjType>().toKernel() )>::type> : std::true_type
{};

/*! \brief has_packRequest check if a type has defined a
 * method called packRequest
 *
 * ### Example
 *
 * \snippet util_test.hpp Check has pack
 *
 * return true if T::packRequest() is a valid expression (function pointers)
 * and produce a defined type
 *
 */

template<typename ObjType, typename Sfinae = void>
struct has_packRequest: std::false_type {};

template<typename ObjType>
struct has_packRequest<ObjType, typename Void<decltype( ObjType::packRequest() )>::type> : std::true_type
{};


/*! \brief has_calculateMem check if a type has defined a
 * method called calculateMem
 *
 * ### Example
 *
 * \snippet util_test.hpp Check has packRequest
 *
 * return true if T::calculateMem() is a valid expression (function pointers)
 * and produce a defined type
 *
 */

template<typename ObjType, typename Sfinae = void>
struct has_packMem: std::false_type {};

/*! \brief has_PackMem check if a type has packMem() member function
 *
 * ### Example
 * \snippet util_test.hpp Check has packMem
 *
 * return true if ObjType::packMem() is a valid expression
 *
 */
template<typename ObjType>
struct has_packMem<ObjType, typename Void<decltype( ObjType::packMem() )>::type> : std::true_type
{};

/*! \brief is_openfpm_native check if a type is an openfpm native structure type
 *
 * ### Example
 * \snippet util_test.hpp Declaration of an openfpm native structure
 * \snippet util_test.hpp Usage with a non openfpm native structure
 *
 * return true if T::attributes::name[0] is a valid expression
 * and produce a defined type
 *
 */
template<typename T, bool = is_typedef_and_data_same<has_typedef_type<T>::value && has_data<T>::value,T>::value>
struct is_openfpm_native : std::false_type
{};


template<typename T>
struct is_openfpm_native<T,true> : std::true_type
{};

template<typename T, typename Sfinae = void>
struct has_value_type: std::false_type {};

/*! \brief has_value_type check if a type has defined a member value_type
 *
 * ### Example
 *
 * \snippet util_test.hpp Check has_value_type
 *
 * return true if T::value_type is a valid type
 *
 */
template<typename T>
//struct has_value_type<T, typename Void<decltype( typename T::value_type )>::type> : std::true_type
struct has_value_type<T, typename Void< typename T::value_type>::type> : std::true_type
{};


//! [Metafunction definition]
template<size_t index, size_t N> struct MetaFunc {
   enum { value = index + N };
};


template<size_t index, size_t N> struct MetaFuncOrd {
   enum { value = index };
};

///////////// Check if the

template<typename ObjType, typename Sfinae = void>
struct isDynStruct: std::false_type
{
	constexpr static bool value()
	{
		return false;
	}
};

template<typename ObjType>
struct isDynStruct<ObjType, typename Void<decltype( ObjType::isCompressed() )>::type> : std::true_type
{
	constexpr static bool value()
	{
		return ObjType::isCompressed();
	}
};


#endif
