/*
 * util.hpp
 *
 *  Created on: Feb 23, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_PLOT_UTIL_HPP_
#define OPENFPM_IO_SRC_PLOT_UTIL_HPP_

/*! \brief It fill the vector x with function values
 *
 * ### Define a function
 * \snippet Plot_unit_tests.hpp Definition of a function
 * ### Example vector with points on a specified range
 * \snippet Plot_unit_tests.hpp fill a vector with a function
 *
 */
template<typename T> static inline void Fill1D(T start, T stop, size_t np, openfpm::vector<T> & x,T f(T x))
{
	x.resize(np);

	T spacing = (stop - start) / (np - 1);

	for (size_t i = 0 ; i < np ; i++)
		x.get(i) = f(start + i*spacing);
}

/*! \brief It fill the vector x with uniformly distributed set of points
 *
 * ### Example vector with points on a specified range
 * \snippet Plot_unit_tests.hpp fill a vector
 *
 */
template<typename T> static inline void Fill1D(T start, T stop, size_t np, openfpm::vector<T> & x)
{
	x.resize(np);

	T spacing = (stop - start) / (np - 1);

	for (size_t i = 0 ; i < np ; i++)
		x.get(i) = start + i*spacing;
}

#endif /* OPENFPM_IO_SRC_PLOT_UTIL_HPP_ */
