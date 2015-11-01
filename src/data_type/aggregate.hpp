/*
 * aggregate.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_AGGREGATE_HPP_
#define OPENFPM_DATA_SRC_UTIL_AGGREGATE_HPP_

#include <boost/fusion/container/vector.hpp>

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
	typedef boost::fusion::vector<list...> type;

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
};


#endif /* OPENFPM_DATA_SRC_UTIL_AGGREGATE_HPP_ */
