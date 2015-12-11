/*
 * mul_array_extents.hpp
 *
 *  Created on: Dec 2, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_MUL_ARRAY_EXTENTS_HPP_
#define OPENFPM_DATA_SRC_UTIL_MUL_ARRAY_EXTENTS_HPP_

/*! \brief Struct that give functionalities on array extensions
 *
 * # Multiplication of all the extents
 * \snippet util_test.hpp Usage mul_array_extents
 *
 */
template<typename T>
struct array_extents
{
	static inline size_t mul()
	{
		return 1;
	}
};

/*! \brief Struct that give functionalities on array extensions
 *
 * # Multiplication of all the extents
 * \snippet util_test.hpp Usage mul_array_extents
 *
 */
template<typename T,size_t N1>
struct array_extents<T[N1]>
{
	static inline size_t mul()
	{
		return N1;
	}
};

/*! \brief Struct that give functionalities on array extensions
 *
 * # Multiplication of all the extents
 * \snippet util_test.hpp Usage mul_array_extents
 *
 */
template<typename T,size_t N1,size_t N2>
struct array_extents<T[N1][N2]>
{
	static inline bool mul()
	{
		return N1 * N2;
	}
};

/*! \brief Struct that give functionalities on array extensions
 *
 * # Multiplication of all the extents
 * \snippet util_test.hpp Usage mul_array_extents
 *
 */
template<typename T,size_t N1,size_t N2,size_t N3>
struct array_extents<T[N1][N2][N3]>
{
	static inline bool mul()
	{
		return N1 * N2 * N3;
	}
};


#endif /* OPENFPM_DATA_SRC_UTIL_MUL_ARRAY_EXTENTS_HPP_ */
