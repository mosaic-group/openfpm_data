/*
 * copy_aggregates.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_COPY_COMPARE_AGGREGATES_HPP_
#define OPENFPM_DATA_SRC_UTIL_COPY_COMPARE_AGGREGATES_HPP_

template<typename T> struct meta_copy;
template<typename T> struct meta_compare;

/*! \brief Structure to copy aggregates
 *
 * \tparam aggregate to copy
 *
 */
template<typename S>
struct copy_aggregate
{
	//! src
	const S & src;

	//! Destination grid
	S & dst;

	//! copy_aggregate
	inline copy_aggregate(const S & src, S & dst)
	:src(src),dst(dst){};

	//! It call the copy function for each member
	template<typename T>
	inline void operator()(T& t) const
	{
		// This is the type of the object we have to copy
		typedef typename boost::fusion::result_of::at_c<typename S::type,T::value>::type copy_type;

		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<copy_type>::type copy_rtype;

		meta_copy<copy_rtype> cp(src.template get<T::value>(),dst.template get<T::value>());
	}
};


/*! \brief Structure to compare aggregates
 *
 * \tparam aggregate to compare
 *
 */
template<typename S>
struct compare_aggregate
{
	bool eq;

	//! src
	const S & src;

	//! Destination grid
	const S & dst;

	//! copy_aggregate
	inline compare_aggregate(const S & src, const S & dst)
	:eq(true),src(src),dst(dst){};

	//! It call the copy function for each member
	template<typename T>
	inline void operator()(T& t)
	{
		// This is the type of the object we have to copy
		typedef typename boost::fusion::result_of::at_c<typename S::type,T::value>::type compare_type;

		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<compare_type>::type compare_rtype;

		if (eq == false)
			return;

		eq = meta_compare<compare_rtype>::meta_compare_f(src.template get<T::value>(),dst.template get<T::value>());
	}

	/*! \brief Returh the result of the comparison
	 *
	 * \return true if aggregates match
	 *
	 */
	inline bool result()
	{
		return eq;
	}
};

#endif /* OPENFPM_DATA_SRC_UTIL_COPY_COMPARE_AGGREGATES_HPP_ */
