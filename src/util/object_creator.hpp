#ifndef VECTOR_CREATOR
#define VECTOR_CREATOR

#include <boost/fusion/container/vector.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/mpl/range_c.hpp>
#include <type_traits>
#include "util_debug.hpp"
#include "check_no_pointers.hpp"
#include "util/for_each_ref.hpp"
#include <iostream>

/*! \brief functor for for_each algorithm
 *
 * It print a warnings for each properties that does not have noPointers method
 *
 * \tparam v original boost::fusion::vector with the properties
 *
 */
template<typename v>
struct check_types
{
	size_t & ret;

	check_types(size_t & ret)
	:ret(ret)
	{
		ret = PNP::NO_POINTERS;
	}

	template<typename T>
	void operator()(T& t)
	{


		typedef typename std::remove_all_extents< typename boost::mpl::at<v,boost::mpl::int_<T::value> >::type>::type tpy;

		// if it is a pointer make no sense
		if (std::is_pointer<tpy>::value == true)
		{ret = PNP::POINTERS;return;}

		// if it is an l-value reference make no send
		if (std::is_lvalue_reference<tpy>::value == true)
		{ret = PNP::POINTERS;return;}

		// if it is an r-value reference make no send
		if (std::is_rvalue_reference<tpy>::value == true)
		{ret = PNP::POINTERS;return;}

		if (std::is_fundamental<tpy>::value == true)
			return;

		// check that T has a method called noPointers
		switch (check_no_pointers<tpy>::value())
		{
			case PNP::UNKNOWN:
			{
				std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " impossible to check the type " << demangle(typeid(tpy).name()) << " please consider to add a static method \"static bool noPointers()\" \n" ;
				ret = PNP::UNKNOWN;
				break;
			}
			case PNP::POINTERS:
			{
				ret = PNP::POINTERS;
				break;
			}
			default:
			{

			}
		}
	}
};

/*! \brief push p_ele into v only of to_push is true
 *
 * \tparam v vector where to push
 * \tparam p_ele object to push
 * \tparam to_push it push the element only if this flag is true
 *
 */
template<typename v,typename p_ele,bool to_push>
struct conditional_push
{
	typedef typename boost::mpl::push_front<v,p_ele>::type type;
};

/*! \brief push p_ele into v only of to_push is true
 *
 * \tparam v vector where to push
 * \tparam p_ele object to push
 * \tparam to_push it push the element only if this flag is true
 *
 */
template<typename v,typename p_ele>
struct conditional_push<v,p_ele,false>
{
	typedef v type;
};

/*! \brief Implementation of noPointer_sequence_impl
 *
 * Case of a property that has not noPointers method
 *
 * \tparam v original boost::fusion::vector
 * \tparam p1 property we are considering
 * \tparam remaining properties to push
 *
 */
template<typename v, int p1, int... prp>
struct noPointers_sequence_impl
{
	// Object we are analyzing
	typedef typename std::remove_reference<typename std::remove_pointer<typename boost::mpl::at< v,boost::mpl::int_<p1> >::type>::type>::type obj_type;

	// analyze the other properties, returning the output vector
	typedef typename noPointers_sequence_impl<v ,prp...>::type v_step;

	// push on the vector the element p1
	typedef typename conditional_push<v_step, boost::mpl::int_<p1>, !has_noPointers<obj_type>::value && !std::is_fundamental<obj_type>::value >::type type;
};

/*! \brief Implementation of noPointer_sequence_impl
 *
 * specialization for last properties
 *
 * \tparam v original boost::fusion::vector
 * \tparam p1 property we are considering
 *
 */
template<typename v, int p1>
struct noPointers_sequence_impl<v,p1>
{
	// Object we are analyzing
	typedef typename std::remove_reference< typename std::remove_pointer<typename boost::mpl::at< v,boost::mpl::int_<p1> >::type>::type>::type obj_type;

	// push on the vector the element p1
	typedef typename conditional_push<boost::mpl::vector<>, boost::mpl::int_<p1>, !has_noPointers<obj_type>::value && !std::is_fundamental<obj_type>::value >::type type;
};

/*! \brief it return a boost::mpl::vector of integers where each integer identify one object without the method "noPointers"
 *
 * ## Example
 *
 * \snippet util.hpp object creator check for no pointers
 *
 * \tparam v boost::fusion::vector
 * \tparam prp selected properties
 *
 */
template<typename v, int... prp>
struct noPointers_sequence
{
	typedef typename noPointers_sequence_impl<v,prp...>::type type;
};


/*! \brief This is a container to create a general object
 *
 * \note It is used in ghost_get to create a particular object with the properties selected
 *
 * For object we mean an object that follow the OpenFPM data structure format, see openFPM_data wiki
 * for more information
 *
 * \tparam Is a boost::fusion::vector with the properties selected
 *
 *
 */
template<typename v>
struct object
{
	typedef v type;
	typedef typename boost::fusion::result_of::size<v>::type size_tpy;

	type data;

	static bool noPointers()
	{
		size_t ret;
		check_types<v> ct(ret);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0, boost::fusion::result_of::size<v>::type::value > > (ct);

		return ret;
	}

	static const int max_prop = size_tpy::value;
};

/*! \brief Implementation of object creator
 *
 * \tparam v boost::fusion::vector
 * \tparam vc basic boost::fusion::vector object from where start to push
 * \tparam prp properties to push
 *
 */
template<typename v, typename vc, int... prp>
struct object_creator_impl
{
};

/*! \brief Implementation of object creator
 *
 * \tparam v original boost::fusion::vector
 * \tparam vc basic boost::fusion::vector object from where start to push
 * \tparam p1 actual property
 * \tparam remaining properties to push
 *
 */
template<typename v, typename vc, int p1, int... prp>
struct object_creator_impl<v,vc,p1,prp...>
{
	typedef typename object_creator_impl<v,vc,prp... >::type vc_step;

	typedef typename boost::remove_reference< typename boost::mpl::at< v,boost::mpl::int_<p1> >::type>::type ele;

	// push on the vector the element p1
	typedef typename boost::mpl::push_front<vc_step, ele >::type type;
};

/*! \brief Implementation of object creator
 *
 * \tparam v original boost::fusion::vector
 * \tparam vc basic boost::fusion::vector object from where start to push
 * \tparam prp property to push
 */
template<typename v, typename vc, int prp>
struct object_creator_impl<v,vc,prp>
{
	typedef typename boost::remove_reference< typename boost::mpl::at< v,boost::mpl::int_<prp> >::type>::type ele;

	// push on the vector the element p1
	typedef typename boost::mpl::push_front<vc, ele >::type type;
};

/*! \brief It create a boost::fusion vector with the selected properties
 *
 * \tparam v boost::fusion::vector
 * \tparam prp selected properties
 *
 * ## Create a compile-time object and copy *from* the selected properties
 * \snippet util_test.hpp object copy example
 * ## Create a compile-time Encap object and copy *from* the selected properties
 * \snippet util_test.hpp object copy encap example
 * ## Create a compile-time object and copy *to* the selected properties
 * \snipper util_test.hpp object write example
 * ## Create a compile-time Encap object and copy *to* the selected properties
 * \snipper util_test.hpp object write encap example
 *
 *
 */
template<typename v, int... prp>
struct object_creator
{
	typedef typename boost::fusion::result_of::as_vector<typename object_creator_impl<v,boost::mpl::vector<>,prp... >::type>::type type;
};

//! specialization when no properties are passed
template<typename v>
struct object_creator<v>
{
	typedef v type;
};

#endif
