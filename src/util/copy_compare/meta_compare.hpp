/*
 * meta_compare.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_META_COMPARE_HPP_
#define OPENFPM_DATA_SRC_UTIL_META_COMPARE_HPP_

#include "compare_general.hpp"

/*! \brief This class compare general objects
 *
 * this function is a general function to compare
 *
 * * primitives
 * * array of primitives
 * * complex objects
 * * aggregates
 *
 * ### Usage of meta copy and compare for primitives
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for primitives
 * ### Usage of meta copy and compare for array of primitives
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for array of primitives
 * ### Usage of meta copy and compare for openfpm aggregates
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for openfpm aggregates
 * ### Usage of meta copy and compare for complex object
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for complex object
 * ### Usage of meta copy and compare for complex aggregates object
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for complex aggregates object
 * ### Usage of meta copy and compare for Point_test
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for Point_test
 *
 */
template<typename T>
struct meta_compare
{
	static inline bool meta_compare_f(const T & src, const T & dst)
	{
		return compare_general<T>::compare_general_f(src,dst);
	}
};

//! Partial specialization for N=1 1D-Array
template<typename T,size_t N1>
struct meta_compare<T[N1]>
{
	static inline bool meta_compare_f(const T src[N1], const T dst[N1])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			if (compare_general<T>::compare_general_f(src[i1],dst[i1]) == false)
				return false;
		}

		return true;
	}
};

//! Partial specialization for N=2 2D-Array
template<typename T,size_t N1,size_t N2>
struct meta_compare<T[N1][N2]>
{
	static inline bool meta_compare_f(const T src[N1][N2], const T dst[N1][N2])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				if (compare_general<T>::compare_general_f(src[i1][i2],dst[i1][i2]) == false)
					return false;
			}
		}

		return true;
	}
};

//! Partial specialization for N=3
template<typename T,size_t N1,size_t N2,size_t N3>
struct meta_compare<T[N1][N2][N3]>
{
	static inline bool meta_compare_f(const T src[N1][N2][N3], const T dst[N1][N2][N3])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					if (compare_general<T>::compare_general_f(src[i1][i2][i3],dst[i1][i2][i3]) == false)
						return false;
				}
			}
		}

		return true;
	}
};




#endif /* OPENFPM_DATA_SRC_UTIL_META_COMPARE_HPP_ */
