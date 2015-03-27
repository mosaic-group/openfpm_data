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
#include "Grid/grid_key.hpp"

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

	inline T get(int i) const
	{
		return boost::fusion::at_c<x>(data)[i];
	}

	/*! \brief Get coordinate
	 *
	 * \param i dimension
	 * \return the i-coordinate of the point
	 *
	 */

	inline T& get(int i)
	{
		return boost::fusion::at_c<x>(data)[i];
	}

	/*! \brief Get the component i
	 *
	 * \return the i-component
	 *
	 */

	inline T& operator[](size_t i)
	{
		return get(i);
	}

	/*! \brief operator= between points
	 *
	 * \param p Point
	 *
	 */
	inline Point<dim,T> & operator=(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			get(i) = p.get(i);
		}

		return *this;
	}

	/*! \brief Multiply each components
	 *
	 * \param p Point
	 *
	 */
	template<typename aT> inline Point<dim,T> operator*(Point<dim,aT> & p)
	{
		Point<dim,T> result;

		for (size_t i = 0 ; i < dim ; i++)
		{
			result.get(i) = get(i) * p.get(i);
		}

		return result;
	}

	/*! \brief Sum each components
	 *
	 * \param p Point
	 *
	 */
	template<typename aT> inline Point<dim,T> & operator+=(Point<dim,aT> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			get(i) += p.get(i);
		}

		return *this;
	}

	/*! \brief Sum each components
	 *
	 * \param p Point
	 *
	 */
	template<typename aT> inline Point<dim,T> operator+(Point<dim,aT> & p)
	{
		Point<dim,T> result;

		for (size_t i = 0 ; i < dim ; i++)
		{
			result.get(i) = get(i) + p.get(i);
		}

		return result;
	}

	/*! \brief divide each component
	 *
	 * \param ar Component wise division
	 *
	 */
	template<typename aT> inline Point<dim,T> operator/(aT (&ar)[dim])
	{
		Point<dim,T> result;

		for (size_t i = 0 ; i < dim ; i++)
		{
			result.get(i) = get(i) / ar[i];
		}

		return result;
	}

	/*! \brief divide each component
	 *
	 * \param c Component wise division
	 *
	 */
	template<typename aT> inline Point<dim,T> operator/(aT c)
	{
		Point<dim,T> result;

		for (size_t i = 0 ; i < dim ; i++)
		{
			result.get(i) = get(i) / c;
		}

		return result;
	}

	/*! \brief Point constructor from point
	 *
	 * \param p the point
	 *
	 */
	Point(const Point<dim,T> && p)
	{
	    for(size_t i = 0; i < dim ; i++)
	    {get(i) = p.get(i);}
	}

	/*! \brief Point constructor from point
	 *
	 * \param p the point
	 *
	 */
	Point(const Point<dim,T> & p)
	{
	    for(size_t i = 0; i < dim ; i++)
	    {get(i) = p.get(i);}
	}

	/*! \brief Constructor from an array
	 *
	 * \param p array with the coordinate of the point
	 *
	 */
	Point(const T (&p)[dim])
	{
	    for(size_t i = 0; i < dim ; i++)
	    {get(i) = p[i];}
	}

	/*! \brief Constructor from a grid_key_dx<dim>
	 *
	 * \param key from where to initialize
	 *
	 */
	Point(grid_key_dx<dim> key)
	{
	    for(size_t i = 0 ; i < dim ; i++)
	    {get(i) = key.k[i];}
	}

	/*! \brief Constructor from a list
	 *
	 * [Example] Point<3,float> p({0.0,0.0,1.0})
	 *
	 */
	Point(std::initializer_list<T> p1)
	{
		size_t i = 0;
	    for(T x : p1)
	    {get(i) = x;i++;}
	}

	//! Default contructor
	Point()
	{}

	static const unsigned int max_prop = 1;
	static const unsigned int dims = dim;
};


#endif
