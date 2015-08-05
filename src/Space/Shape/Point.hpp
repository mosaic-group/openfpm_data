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
#include "Grid/Encap.hpp"

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

	typedef T coord_type;

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

	/*! \brief Get the component i
	 *
	 * \return the i-component
	 *
	 */

	inline const T operator[](size_t i) const
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

	/*! \brief Compare all components
	 *
	 * \param p Point
	 *
	 */
	template<typename aT> inline bool operator==(const Point<dim,aT> & p)
	{
		bool equal = true;

		for (size_t i = 0 ; i < dim ; i++)
		{
			if(get(i) != p.get(i))
				equal = false;
		}

		return equal;
	}

	/*! \brief Multiply each components
	 *
	 * \param p Point
	 *
	 */
	template<typename aT> inline Point<dim,T> operator*(const Point<dim,aT> & p)
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
	template<typename aT> inline Point<dim,T> & operator+=(const Point<dim,aT> & p)
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
	template<typename aT> inline Point<dim,T> operator+(const Point<dim,aT> & p)
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
	template<typename aT> inline Point<dim,T> operator/(const aT (&ar)[dim])
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
	template<typename aT> inline Point<dim,T> operator/(const aT c)
	{
		Point<dim,T> result;

		for (size_t i = 0 ; i < dim ; i++)
		{
			result.get(i) = get(i) / c;
		}

		return result;
	}

	/*! \brief Operator subtraction
	 *
	 *  it produce a point that is the subtraction of two points
	 *
	 * \param p Point
	 *
	 */
	inline Point<dim,T> operator-(const Point<dim,T> & p) const
	{
		Point<dim,T> result;

		for (size_t i = 0 ; i < dim ; i++)
		{
			result.get(i) = get(i) - p.get(i);
		}

		return result;
	}

	/*! \brief Set to zero the point coordinate
	 *
	 *
	 */
	inline void zero()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			get(i) = 0;
		}
	}

	/*! \brief Create a point set to zero
	 *
	 * \return a point with all coorfinate set to 0
	 *
	 */
	inline static Point<dim,T> zero_p()
	{
		Point<dim,T> p;

		for (size_t i = 0 ; i < dim ; i++)
		{
			p.get(i) = 0;
		}

		return p;
	}

	/*! \brief Convert the point into a string
	 *
	 * \return the string
	 *
	 */
	std::string toPointString() const
	{
		std::stringstream ps;

		for (size_t i = 0 ; i < dim ; i++)
			ps << "x[" << i << "]=" << get(i) << " ";

		ps << "\n";

		return ps.str();
	}

	/*! \brief Return the string with the point coordinate
	 *
	 * \return the string
	 *
	 */
	std::string toString() const
	{
		std::string str;

		for (size_t i = 0 ; i < dim - 1 ; i++)
		{
			str += std::to_string(get(i)) + " ";
		}
		str += std::to_string(get(dim-1));

		return str;
	}

	/*! \brief Point constructor from point
	 *
	 * \param p the point
	 *
	 */
	inline Point(const Point<dim,T> && p)
	{
	    for(size_t i = 0; i < dim ; i++)
	    {get(i) = p.get(i);}
	}

	/*! \brief Point constructor from point
	 *
	 * \param p the point
	 *
	 */
	inline Point(const Point<dim,T> & p)
	{
	    for(size_t i = 0; i < dim ; i++)
	    {get(i) = p.get(i);}
	}

	/*! \brief Constructor from an array
	 *
	 * \param p array with the coordinate of the point
	 *
	 */
	inline Point(const T (&p)[dim])
	{
	    for(size_t i = 0; i < dim ; i++)
	    {get(i) = p[i];}
	}

	/*! \brief Constructor from a grid_key_dx<dim>
	 *
	 * \param key from where to initialize
	 *
	 */
	inline Point(grid_key_dx<dim> key)
	{
	    for(size_t i = 0 ; i < dim ; i++)
	    {get(i) = key.k[i];}
	}

	/*! \brief Constructor from a list
	 *
	 * [Example] Point<3,float> p({0.0,0.0,1.0})
	 *
	 */
	inline Point(std::initializer_list<T> p1)
	{
		size_t i = 0;
	    for(T x : p1)
	    {get(i) = x;i++;}
	}

	//! Default contructor
	inline Point()
	{}

	static const unsigned int max_prop = 1;
	static const unsigned int dims = dim;
};

/*! \brief Convert an array of point coordinate into string
 *
 * \param p coordinate on each dimension
 *
 * \return the string
 *
 */
template <unsigned int N, typename T> std::string toPointString(const T (&p)[N] )
{
	std::stringstream ps;

	for (size_t i = 0 ; i < N ; i++)
		ps << "x[" << i << "]=" << p[i] << " ";

	ps << "\n";

	return ps.str();
}

/*! \brief Convert an encapsulated point into string
 *
 * \param p coordinate on each dimension
 *
 * \return the string
 *
 */
template <unsigned int N, typename T, typename Mem> std::string toPointString(const encapc<1,Point<N,T>,Mem> & p )
{
	std::stringstream ps;

	for (size_t i = 0 ; i < N ; i++)
		ps << "x[" << i << "]=" << p.template get<Point<N,T>::x>()[i] << " ";

	ps << "\n";

	return ps.str();
}

#endif
