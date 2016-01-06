/*
 * compare_fusion_vector.hpp
 *
 *  Created on: Nov 1, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_COMPARE_FUSION_VECTOR_HPP_
#define OPENFPM_DATA_SRC_UTIL_COMPARE_FUSION_VECTOR_HPP_

/*! \brief this class is a functor for "for_each" algorithm
 *
 * It compare a boost::fusion::vector with another boost::fusion::vector
 *
 */
template<typename bfv>
struct compare_fusion_vector
{
	bool eq;

	const bfv & src;
	const bfv & dst;

	/*! \brief constructor
	 *
	 * It define the elements to compare.
	 *
	 * \param obj object we have to set in grid_dst
	 *
	 */
	inline compare_fusion_vector(const bfv & src, const bfv & dst)
	:eq(true),src(src),dst(dst){};

#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline compare_fusion_vector(const bfv && src, const bfv && dst)
	:eq(true),src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		// This is the type of the object we have to copy
		typedef typename boost::fusion::result_of::at_c<bfv,T::value>::type copy_type;

		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<copy_type>::type copy_rtype;

		if (eq == false)
			return;

		eq = meta_compare<copy_rtype>::meta_compare_f(boost::fusion::at_c<T::value>(src),boost::fusion::at_c<T::value>(dst));
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



#endif /* OPENFPM_DATA_SRC_UTIL_COMPARE_FUSION_VECTOR_HPP_ */
