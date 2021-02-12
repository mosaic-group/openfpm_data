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
#include "util/copy_compare/copy_compare_aggregates.hpp"

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

#ifdef SE_CLASS3

#define SE3_MAX_PROP(i) i+2
#define SE3_ADD_PROP(i) size_t[i+1],size_t
#define SE3_SUB_MAX_PROP -2

/*! \brief An aggregate that accept a boost fusion vector as type
 *
 *
 *
 */
template<typename T>
struct aggregate_bfv
{
	//! type the object store
	typedef T type;

	//! real type the object store
	typedef T type_real;

	//! data to store
	type data;

	aggregate_bfv()	{};

	static const unsigned int max_prop = boost::mpl::size<type>::type::value;
	static const unsigned int max_prop_real = boost::mpl::size<type>::type::value + SE3_SUB_MAX_PROP;
};

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

	__device__ __host__ inline aggregate()
	{}

	__device__ __host__ inline aggregate(const aggregate<list ...> & aggr)
	{
		this->operator=(aggr);
	}

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

	/*! \brief it return false if this aggregate has no pointers
	 *
	 *
	 */
	static bool noPointers()
	{
		return !has_pack_gen<aggregate<list ...>>::value;
	}

	aggregate<list...> & operator=(const aggregate<list...> & ag)
	{
		copy_fusion_vector<aggregate<list...>::type> ca(ag.data,this->data);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(list)> >(ca);

		return *this;
	}

	static const unsigned int max_prop = boost::mpl::size<type>::type::value;
	static const unsigned int max_prop_real = boost::mpl::size<type>::type::value + SE3_SUB_MAX_PROP;
};



#else

/*! \brief An aggregate that accept a boost fusion vector as type
 *
 *
 *
 */
template<typename T>
struct aggregate_bfv
{
	//! type the object store
	typedef T type;

	//! real type the object store
	typedef T type_real;

	//! data to store
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

	static const unsigned int max_prop = boost::mpl::size<type>::type::value;
};

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

	__device__ __host__ inline aggregate()
	{}

	__device__ __host__ inline aggregate(const aggregate<list ...> & aggr)
	{
		this->operator=(aggr);
	}

	/*! \brief get the properties i
	 *
	 * \return the property i
	 *
	 */
	template<unsigned int i> __device__ __host__ typename boost::mpl::at<type,boost::mpl::int_<i>>::type & get()
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
	template<unsigned int i> __device__ __host__ const typename boost::mpl::at<type,boost::mpl::int_<i>>::type & get() const
	{
		return boost::fusion::at_c<i>(data);
	}

	inline aggregate<list...> & operator=(const aggregate<list...> & ag)
	{
		copy_fusion_vector<aggregate<list...>::type> ca(ag.data,this->data);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(list)> >(ca);

		return *this;
	}

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

namespace openfpm
{
	template<unsigned int p, typename aggr>
	auto at_c(aggr & agg) -> decltype(boost::fusion::at_c<p>(agg.data))
	{
		return boost::fusion::at_c<p>(agg.data);
	}

	template<unsigned int p, typename aggr>
	auto get(aggr & agg) -> decltype(boost::fusion::at_c<p>(agg.data))
	{
		return boost::fusion::at_c<p>(agg.data);
	}
}

template<typename BlockT, typename T>
struct AggregateAppend
{
};

template<typename BlockT, typename ... list>
struct AggregateAppend<BlockT, aggregate<list ...>>
{
    typedef aggregate<list..., BlockT> type;
};

#endif /* OPENFPM_DATA_SRC_UTIL_AGGREGATE_HPP_ */
