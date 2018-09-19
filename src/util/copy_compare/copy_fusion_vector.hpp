/*
 * copy_fusion_vector.hpp
 *
 *  Created on: Nov 1, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_COPY_FUSION_VECTOR_HPP_
#define OPENFPM_DATA_SRC_GRID_COPY_FUSION_VECTOR_HPP_

/*! \brief this class is a functor for "for_each" algorithm
 *
 * It copy a boost::fusion::vector into another boost::fusion::vector
 *
 */
template<typename bfv>
struct copy_fusion_vector
{
	//! source fusion vector
	const bfv & src;

	//! destination fusion vector
	bfv & dst;

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param src source fusion vector
	 * \param dst destination fusion vector
	 *
	 */
	__device__ __host__ inline copy_fusion_vector(const bfv & src, bfv & dst)
	:src(src),dst(dst){};

#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_fusion_vector(const bfv && src, bfv && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t)
	{
		// This is the type of the object we have to copy
		typedef typename boost::fusion::result_of::at_c<bfv,T::value>::type copy_type;

		// Remove the reference from the type to copy
		typedef typename boost::remove_reference<copy_type>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(boost::fusion::at_c<T::value>(src),boost::fusion::at_c<T::value>(dst));
	}
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * It copy a boost::fusion::vector into an encap
 *
 */
template<typename bfv, typename enc>
struct copy_fusion_vector_encap
{
	//! source fusion vector
	const bfv & src;

	//! destination fusion vector
	enc & dst;

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param src source fusion vector
	 * \param dst destination fusion vector
	 *
	 */
	inline copy_fusion_vector_encap(const bfv & src, enc & dst)
	:src(src),dst(dst){};

#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_fusion_vector_encap(const bfv && src, enc && dst)
	:src(src),dst(dst)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object\n";};
#endif

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<bfv,boost::mpl::int_<T::value> >::type copy_src;
		typedef typename std::remove_reference<decltype(dst.template get<T::value>())>::type copy_dst;

		meta_copy_d<copy_src,copy_dst>::meta_copy_d_(boost::fusion::at_c<T::value>(src),dst.template get<T::value>());
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_COPY_FUSION_VECTOR_HPP_ */
