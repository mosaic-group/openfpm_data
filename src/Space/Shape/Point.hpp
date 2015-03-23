#ifndef POINT_HPP
#define POINT_HPP

#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include "boost/multi_array.hpp"
#include "base_type.hpp"
#include "memory_conf.hpp"

/*! \brief This class implement the point shape in an N-dimensional space
 *
 * This class implement the point shape in an N-dimensional space
 *
 * \param T type of the space
 * \param dim dimensionality
 *
 */

template<unsigned int dim ,typename T> class Point
{
	public:

	//! boost fusion that store the point
	typedef boost::fusion::vector<T[dim]> type;
	//! layout that interleave the properties
	typedef typename memory_traits_inte<type>::type memory_int;
	//! layout with linear properties
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! structure that store the data of the point
	type data;

	//! Property id of the point
	static const unsigned int x = 0;

	/*! \brief Get coordinate
	 *
	 * \param i dimension
	 * \return the i-coordinate of the point
	 *
	 */

	T get(int i)
	{
		return boost::fusion::at_c<x>(data)[i];
	}

	static const unsigned int max_prop = 1;
	static const unsigned int dims = dim;
};


#endif
