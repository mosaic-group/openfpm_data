/*
 * map_vector_std_util.hpp
 *
 *  Created on: May 14, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_VECTOR_MAP_VECTOR_STD_UTIL_HPP_
#define OPENFPM_DATA_SRC_VECTOR_MAP_VECTOR_STD_UTIL_HPP_

#include "util/common.hpp"

/*! \brief pack/add function selector
 *
 * This in case of error
 *
 */
template<bool,typename T, typename S>
struct push_back_op_neste
{
	/*! \brief push_back
	 *
	 * Print an error
	 *
	 * \param base vector where to push_back
	 * \param obj object to push
	 *
	 */
	static inline void push_back(std::vector<T> & base, const S & obj)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " error cannot push " << demangle(typeid(S).name()) << " into a vector of " << demangle(typeid(T).name()) << std::endl;
	}
};

/*! \brief pack/add function selector
 *
 * This in case normally add an object
 *
 */
template<typename T, typename S>
struct push_back_op_neste<true,T,S>
{

	/*! \brief Push_back on a vector
	 *
	 * Push on a vector an object
	 *
	 * \param base vector where to push_back
	 * \param obj object to push
	 *
	 */
	static inline void push_back(std::vector<T> & base, const S & obj)
	{
		base.push_back(T());
		base[base.size()-1] = obj;
	}
};

/*! \brief pack/add function selector
 *
 * In case of error
 *
 */
template<bool,typename T, typename S>
struct push_back_std_op_neste
{

	/*! \brief push_back
	 *
	 * Print an error
	 *
	 * \param base vector where to push_back
	 * \param obj object to push
	 *
	 */
	static inline void push_back(std::vector<T> & base, const S & obj)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " error cannot push " << demangle(typeid(S).name()) << " into a vector of " << demangle(typeid(T).name()) << std::endl;
	}
};

template<typename T, typename S>
struct push_back_std_op_neste<true,T,S>
{
	static inline void push_back(std::vector<T> & base, const S & obj)
	{
		base.reserve(base.size() + obj.size());
		for (size_t i = 0 ; i < obj.size() ; i++)
		{
			base.push_back(obj.get(i));
		}
	}
};

template<bool is_t, bool is_s,typename T, typename S>
struct push_back_op
{
	static inline void push_back(std::vector<T> & base, const S & obj)
	{std::cerr << __FILE__ << ":" << __LINE__ << " error cannot push " << demangle(typeid(S).name()) << " into a vector of " << demangle(typeid(T).name()) << std::endl;}
};

template<typename T, typename S>
struct push_back_op<false,false,T,S>
{
	static inline void push_back(std::vector<T> & base, const S & obj)
	{
		base.push_back(obj);
	}
};

template<typename T, typename S>
struct push_back_op<false,true,T,S>
{
	static inline void push_back(std::vector<T> & base, const S & obj)
	{
		push_back_std_op_neste<std::is_same<T,typename S::value_type>::value,T,S>::push_back(base,obj);
	}
};



template<typename T, typename S>
struct push_back_op<true,true,T,S>
{
	static inline void push_back(std::vector<T> & base, const S & obj)
	{push_back_op_neste<std::is_assignable<typename T::value_type,typename S::value_type>::value ||
		                std::is_same<typename T::value_type,typename S::value_type>::value,T,S>::push_back(base,obj);}
};


template<bool has_base, typename base_obj, typename v_obj>
struct base_copy
{
	static inline void copy(base_obj & base, const v_obj & obj)
	{
		base.clear();
		base.resize(obj.size());

		for (size_t i = 0 ; i <  obj.size() ; i++)
		{
			base.get(i) = obj.get(i);
		}
	}
};

template<typename base_obj, typename v_obj>
struct base_copy<true,base_obj,v_obj>
{
	static inline void copy(base_obj & base, const v_obj & obj)
	{
		base = obj.base;
	}
};

template<typename T, typename Sfinae = void>
struct has_base_to_copy: std::false_type {};

/*! \brief has_data check if a type has defined a member data
 *
 * ### Example
 *
 * \snippet map_vector_std_util_unit_test.hpp Check has_base_to_copy
 *
 * return true if T::type is a valid type
 *
 */
template<typename T>
struct has_base_to_copy<T, typename Void< typename T::base_to_copy >::type > : std::true_type
{};

#endif /* OPENFPM_DATA_SRC_VECTOR_MAP_VECTOR_STD_UTIL_HPP_ */
