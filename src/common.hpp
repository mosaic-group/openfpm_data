#ifndef COMMON_HPP
#define COMMON_HPP

#include <type_traits>
#include <random>

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

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element the operator() is called. Used
 * to calculate the size of the selected elements
 *
 * \tparam N number of properties
 * \tparam v boost::fusion::vector
 *
 */

/*template<unsigned int N, typename v>
struct el_size
{
	//! total_size
	size_t total_size;

	//! List of properties
	const size_t (& prp)[N];*/

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param key which element we are modifying
	 * \param grid_src grid we are updating
	 * \param obj object we have to set in grid_src
	 *
	 */
/*	el_size(const size_t (& prp)[N])
	:total_size(0),prp(prp)
	{};

	//! It call the copy function for each property
    template<typename T>
    void operator()(T& t) const
    {
    	total_size += boost::fusion::result_of::at< v,boost::mpl::int_<prp[T::value]> >::type;
    }

    size_t size()
    {return total_size;}
};*/

/*template<unsigned int N,typename v> size_t ele_size(const size_t (& prp)[N])
{
	el_size<N,v> sz(prp);

	boost::mpl::for_each_ref< boost::mpl::range_c<int,0,N> >(sz);

	return sz.size();
}*/

/**/

#endif
