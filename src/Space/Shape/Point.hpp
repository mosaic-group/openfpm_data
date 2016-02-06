#ifndef POINT_HPP
#define POINT_HPP

#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include "boost/multi_array.hpp"
#include "Grid/Encap.hpp"

/*! \brief This class implement the point shape in an N-dimensional space
 *
 * \param T type of the space
 * \param dim dimensionality
 *
 */

template<unsigned int dim ,typename T> class Point
{
	public:

	//! coordinate type
	typedef T coord_type;

	//! boost fusion that store the point
	typedef boost::fusion::vector<T[dim]> type;

	//! structure that store the data of the point
	type data;

	//! Property id of the point
	static const unsigned int x = 0;

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

	/*! \brief Point constructor
	 *
	 * \param p Point
	 *
	 */
	template <typename S> inline Point(const Point<dim,S> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
			get(i) = static_cast<S>(p.get(i));
	}

	/*! \brief Point constructor
	 *
	 * \param p encapc Point
	 *
	 */
	template <unsigned int d, typename M> inline Point(const encapc<d,Point<dim,T>,M> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
			get(i) = p.template get<0>()[i];
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

	/*! \brief Multiply each components by a constant
	 *
	 * \param c constanr
	 *
	 */
	inline Point<dim,T> operator*(T c)
	{
		Point<dim,T> result;

		for (size_t i = 0 ; i < dim ; i++)
		{
			result.get(i) = get(i) * c;
		}

		return result;
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

	/*! \brief Subtract each components
	 *
	 * \param p Point
	 *
	 */
	inline Point<dim,T> & operator-=(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			get(i) -= p.get(i);
		}

		return *this;
	}

	/*! \brief Sum each components
	 *
	 * \param p Point
	 *
	 */
	inline Point<dim,T> & operator+=(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			get(i) += p.get(i);
		}

		return *this;
	}

	/*! \Brief norm of the vector
	 *
	 * \return the norm of the vector
	 *
	 */
	T norm()
	{
		T n = 0.0;

		for (size_t i = 0 ; i < dim ; i++)
		{
			n+=get(i) * get(i);
		}

		return sqrt(n);
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

	/*! \brief  It calculate the distance between 2 points
	 *
	 * Itself (p) and the other point (q)
	 *
	 * \parameter q target point
	 *
	 * \return the distance
	 *
	 */
	T distance(const Point<dim,T> & q)
	{
		T tot = 0.0;

		for (size_t i = 0 ; i < dim ; i++)
			tot += (this->get(i)  - q.get(i)) * (this->get(i)  - q.get(i));

		return sqrt(tot);
	}

	/*! \brief  It calculate the square distance between 2 points
	 *
	 * Itself (p) and the other point (q)
	 *
	 * \parameter q target point
	 *
	 * \return the square of the distance
	 *
	 */
	T distance2(const Point<dim,T> & q)
	{
		T tot = 0.0;

		for (size_t i = 0 ; i < dim ; i++)
			tot += (this->get(i)  - q.get(i)) * (this->get(i)  - q.get(i));

		return tot;
	}

	/*! \brief Operator subtraction
	 *
	 *  it produce a point that is the subtraction of two points
	 *
	 * \param p Point
	 *
	 */
	inline Point<dim,T> operator-(const Point<dim,T> & p)
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

	/*! \brief Set to one the point coordinate
	 *
	 *
	 */
	inline void one()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			get(i) = 1;
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

	/*! \brief Return the reference to the value at coordinate i
	 *
	 * \return the reference
	 *
	 */
	T & value(size_t i)
	{
		return get(i);
	}

	/*! \brief Return the value at coordinate i
	 *
	 * \return the value
	 *
	 */
	T value(size_t i) const
	{
		return get(i);
	}

	/*! \brief Return the coordinated of the point as reference array
	 *
	 * \return the reference array
	 *
	 */
	T (&asArray())[dim]
	{
		return boost::fusion::at_c<x>(data);
	}

	/*! Convert the point from Point<dim,T> to Point<dim,A>
	 *
	 * \return the converted point
	 *
	 */
	template<typename A> Point<dim,A> convertPoint() const
	{
		Point<dim,A> p;

		for (size_t i = 0; i < dim ; i++)
			p.get(i) = static_cast<A>(get(i));

		return p;
	}

	//! This structure has no internal pointers
	static bool noPointers()
	{
		return true;
	}

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
