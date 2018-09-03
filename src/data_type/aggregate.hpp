/*
 * aggregate.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_AGGREGATE_HPP_
#define OPENFPM_DATA_SRC_UTIL_AGGREGATE_HPP_

#include <boost/fusion/container/vector.hpp>
#include <Packer_Unpacker/has_pack_agg.hpp>

#ifdef SE_CLASS3

#define SE3_MAX_PROP(i) i+2
#define SE3_ADD_PROP(i) size_t[i+1],size_t
#define SE3_SUB_MAX_PROP -2

/*! \brief aggregate of properties, from a list of object if create a struct that follow the OPENFPM native structure
 *
 * see the Wiki for more information about the OPENFPM native structure format
 *
 * \tparam list of properties
 *
 */
template<typename ... list>
struct aggregate
{
	typedef boost::fusion::vector<list... , SE3_ADD_PROP(sizeof...(list))> type;
	typedef boost::fusion::vector<list... > type_real;

	typedef int yes_is_aggregate;

	type data;

	/*! \brief get the properties i
	 *
	 * \return the property i
	 *
	 */
	template<unsigned int i> typename boost::mpl::at<type,boost::mpl::int_<i>>::type & get()
	{
		return boost::fusion::at_c<i>(data);
	}

	/*! \brief it return false if this aggregate has no pointers
	 *
	 *
	 */
	static bool noPointers()
	{
		return !has_pack_gen<aggregate<list ...>>::value;
	}

	/*! \brief get the properties i
	 *
	 * \return the property i
	 *
	 */
	template<unsigned int i> const typename boost::mpl::at<type,boost::mpl::int_<i>>::type & get() const
	{
		return boost::fusion::at_c<i>(data);
	}

	aggregate<list...> & operator=(aggregate<list...> & ag) = delete;

	static const unsigned int max_prop = boost::mpl::size<type>::type::value;
	static const unsigned int max_prop_real = boost::mpl::size<type>::type::value + SE3_SUB_MAX_PROP;
};



#else


/*! \brief aggregate of properties, from a list of object if create a struct that follow the OPENFPM native structure
 *
 * see the Wiki for more information about the OPENFPM native structure format
 *
 * \tparam list of properties
 *
 */
template<typename ... list>
struct aggregate
{
	//! internal type containing the data
	typedef boost::fusion::vector<list...> type;

	//! real internal type containing the data
	typedef boost::fusion::vector<list...> type_real;

	typedef int yes_is_aggregate;

	//! the data
	type data;

	/*! \brief get the properties i
	 *
	 * \return the property i
	 *
	 */
	template<unsigned int i> typename boost::mpl::at<type,boost::mpl::int_<i>>::type & get()
	{
		return boost::fusion::at_c<i>(data);
	}

	/*! \brief it return false if this aggregate has no pointers
	 *
	 *
	 */
	static bool noPointers()
	{
		return !has_pack_gen<aggregate<list ...>>::value;
	}

	/*! \brief get the properties i
	 *
	 * \return the property i
	 *
	 */
	template<unsigned int i> const typename boost::mpl::at<type,boost::mpl::int_<i>>::type & get() const
	{
		return boost::fusion::at_c<i>(data);
	}

	aggregate<list...> & operator=(aggregate<list...> & ag) = delete;

	static const unsigned int max_prop = boost::mpl::size<type>::type::value;
	static const unsigned int max_prop_real = boost::mpl::size<type>::type::value;
};

#endif

template<typename T, typename Sfinae = void>
struct is_aggregate: std::false_type {};


/*! \brief Check if a type T is an aggregate
 *
 * return true if T is an aggregate
 *
 */
template<typename T>
struct is_aggregate<T, typename Void< typename T::yes_is_aggregate>::type> : std::true_type
{};


#endif /* OPENFPM_DATA_SRC_UTIL_AGGREGATE_HPP_ */
