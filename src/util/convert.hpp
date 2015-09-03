/*
 * convert.hpp
 *
 *  Created on: Aug 26, 2015
 *      Author: i-bird
 */

#ifndef SRC_UTIL_CONVERT_HPP_
#define SRC_UTIL_CONVERT_HPP_

/*! Set of utility functions to convert to a point from other data types
 *
 * \tparam dimensionality of the Point
 * \tparam St type of point
 *
 */

template<unsigned int dim, typename St>
class toPoint
{
public:
	/*! \brief Return the combination converted to point
	 *
	 */
	static inline Point<dim,St> convert(const comb<dim> & c)
	{
		// Point
		Point<dim,St> ret;

		// set the point
		for (size_t i = 0; i < dim ; i++)
			ret.get(i) = c.c[i];

		return ret;
	}
};



#endif /* SRC_UTIL_CONVERT_HPP_ */
