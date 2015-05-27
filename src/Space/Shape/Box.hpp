
#ifndef BOX_HPP_
#define BOX_HPP_

#include "Space/Shape/Sphere.hpp"
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include "Grid/grid_key.hpp"
#include "Grid/Encap.hpp"

/*! \brief It define if we want the upper base or the down base (Lower or upper)
 * extreme of the interval
 *
 */

enum Base
{
	UP,
	DOWN
};

/** \brief This class represent an N-dimensional box
 *
 * This class represent an N-dimensional box
 *
 * \tparam dim dimansionality of the Box
 * \tparam T type of space ... Real Complex Integer
 *
 */

template<unsigned int dim , typename T>
class Box
{
public:

	//! boost fusion that store the point
	typedef boost::fusion::vector<T[dim],T[dim]> type;
	//! type of the box
	typedef T btype;
	//! layout that interleave the properties
	typedef typename memory_traits_inte<type>::type memory_int;
	//! layout with linear properties
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! It store the two point bounding the box
	type data;

	//! Low point
	static const unsigned int p1 = 0;
	//! High point
	static const unsigned int p2 = 1;
	//! Maximum number of properties
	static const unsigned int max_prop = 2;

	static const unsigned int dims = dim;

	/*! \brief Intersect
	 *
	 * Intersect two boxes and return the result boxes, if the boxes does not intersect, return false
	 *
	 * \param b box to intersect with
	 * \param b_out box result of the intersection
	 *
	 * \return true if they intersect
	 *
	 */
	bool Intersect(const Box<dim,T> & b, Box<dim,T> & b_out)
	{
		// check if p1 of b is smaller than

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (getLow(i) <= b.getLow(i))
				b_out.setLow(i,b.getLow(i));
			else if (getLow(i) <= b.getHigh(i))
				b_out.setLow(i,getLow(i));
			else
				return false;

			if (getHigh(i) >= b.getHigh(i))
				b_out.setHigh(i,b.getHigh(i));
			else if (getHigh(i) >= b.getLow(i))
				b_out.setHigh(i,getHigh(i));
			else
				return false;
		}
		return true;
	}

	/*! \brief Intersect
	 *
	 * Intersect two boxes and return the result boxes, if the boxes does not intersect, return false
	 *
	 * \param e_b encapsulator box to intersect with
	 * \param b_out box result of the intersection
	 *
	 * \return true if they intersect
	 *
	 */
	template<typename Mem> bool Intersect(const encapc<1,Box<dim,T>,Mem> & e_b, Box<dim,T> & b_out)
	{
		return Intersect(e_b,b_out);
	}

	/*!
	 *
	 * \brief Check if the sphere intersect the box
	 *
	 *
	 *
	 *   p1 _____
	 *      |    |
	 * p.   |    |
	 *      |____|
	 *            p2
	 *
	 * Given a point p and we search for the distance of the nearest point of the box
	 * from p. In this case the distance is p1.x-p0.x the Y dimension is alligned with p
	 * In general for each dimension we check if the point is in the interval if it is
	 * we do not accumulate, otherwise we accumulate the smallest between (p1-p0) (p2-p0).
	 *
	 * if the distance of the nearest point is smaller than the radius we have an intersection
	 *
	 *
	 *
	 * \tparam distance functor
	 * \param sphere to check the intersection
	 * \return true if intersect false otherwise
	 *
	 */
	template <typename distance> bool Intersect(Sphere<dim,T> & sphere)
	{
		// distance functor
		distance dist;

		// Get the nearest point of the box from the center of the sphere
		typename distance::ResultType distance_r = 0;

		for (size_t i = 0 ; i < dim ; i++)
		{

			// if the center of the sphere on dimension i is not in the i box interval
			// do not accumulate, otherwise accumulate from the nearest point on that
			// dimension
			if (boost::fusion::at_c<p1>(data)[i] < sphere.center(i))
			{
				// accumulate the distance from p1
				distance_r += dist.accum_dist(sphere.center(i),boost::fusion::at_c<p1>(data)[i],i);
			}
			else if ( boost::fusion::at_c<p2>(data)[i] <= sphere.center(i))
			{
				// accumulate the distance from p2
				distance_r += dist.accum_dist(sphere.center(i),boost::fusion::at_c<p2>(data)[i],i);
			}
		}

		// return if there is intersection
		return distance_r < sphere.radius();
	}

	/*! \brief Get the coordinate of the bounding point
	 *
	 * \tparam b integer define (lower or upper interval)
	 * \param i i-coordinate of the point
	 * \return the i-coordinate of the bounding point
	 *
	 */

	template<unsigned int b> T getBase(const unsigned int i)
	{
		return boost::fusion::at_c<b>(data)[i];
	}

	/*! \brief Get the the box of dimensionality dim-1 (it eliminate the last dimension)
	 *
	 *  \return Return the sub-box of dimension dim-1
	 *
	 */

	Box<dim-1,T> getSubBox()
	{
		return Box<dim-1,T>(data);
	}

	/*! \brief Operator= between boxes
	 *
	 * Operator= between boxes of the same size
	 *
	 * \param box is the box that store the interval
	 * \return itself
	 *
	 */

	Box<dim,T> & operator=(const Box<dim,T> & box)
	{
	    for(size_t i = 0 ; i < dim ; i++)
	    {setLow(i,box.getLow(i));}

	    for(size_t i = 0 ; i < dim ; i++)
	    {setHigh(i,box.getHigh(i));}

		// return itself
		return *this;
	}

	public:

	//! default constructor
	Box()
	{}

	/*! \brief Constructor from initializer list
	 *
	 * Constructor from initializer list
	 *
	 * \param p1 Low point, initialize as a list example {0.0,0.0,0.0}
	 * \param p2 High point, initialized as a list example {1.0,1.0,1.0}
	 *
	 */

	Box(std::initializer_list<T> p1, std::initializer_list<T> p2)
	{
		set(p1,p2);
	}

	/*! \brief Box constructor from a box
	 *
	 * Box constructor from a box
	 *
	 * \param high array indicating the coordinates of the low point
	 * \param low array indicating the coordinates of the high point
	 *
	 */
	inline Box(T * high, T * low)
	{
		// copy all the data

		for (int i = 0 ; i < dim ; i++)
		{
			// p1 is low p2 is high

			boost::fusion::at_c<Box::p1>(data)[i] = low[i];
			boost::fusion::at_c<Box::p2>(data)[i] = high[i];
		}
	}

	/*! \brief Box constructor from a box
	 *
	 * Box constructor from a box
	 *
	 * \param box from which to construct
	 *
	 */
	inline Box(const Box<dim,T> & box)
	{
		// we copy the data

		for (size_t i = 0 ; i < dim ; i++)
		{
			boost::fusion::at_c<p1>(data)[i] = boost::fusion::at_c<p1>(box.data)[i];
			boost::fusion::at_c<p2>(data)[i] = boost::fusion::at_c<p2>(box.data)[i];
		}
	}

	/*! \brief Box constructor from vector::fusion
	 *
	 * Box constructor from vector::fusion
	 *
	 * \param box_data from which to construct
	 *
	 */
	inline Box(type box_data)
	{
		// we copy the data

		for (size_t i = 0 ; i < dim ; i++)
		{
			boost::fusion::at_c<p1>(data)[i] = boost::fusion::at_c<p1>(box_data)[i];
			boost::fusion::at_c<p2>(data)[i] = boost::fusion::at_c<p2>(box_data)[i];
		}
	}

	/*! \brief Box constructor from an array reference
	 *
	 * \param array from which to construct the box
	 *
	 */
	inline Box(T (& box_data)[dim])
	{
		// we copy the data

		for (size_t i = 0 ; i < dim ; i++)
		{
			boost::fusion::at_c<p1>(data)[i] = 0;
			boost::fusion::at_c<p2>(data)[i] = box_data[i];
		}
	}

	/*! \brief Box constructor from vector::fusion of higher dimension
	 *
	 * \param box_data fusion vector from which to construct the vector
	 *
	 */

	template<unsigned int dimS> inline Box(boost::fusion::vector<T[dimS],T[dimS]> & box_data)
	{
		// we copy the data

		for (size_t i = 0 ; i < dim ; i++)
		{
			boost::fusion::at_c<p1>(data)[i] = boost::fusion::at_c<p1>(box_data)[i];
			boost::fusion::at_c<p2>(data)[i] = boost::fusion::at_c<p2>(box_data)[i];
		}
	}

	/*! \brief Box constructor from ecapsulated box
	 *
	 * \param box_data fusion vector from which to construct the vector
	 *
	 */

	template<typename Mem> inline Box(const encapc<1,Box<dim,T>,Mem> & b)
	{
		// we copy the data

		for (size_t i = 0 ; i < dim ; i++)
		{
			boost::fusion::at_c<p1>(data)[i] = b.template get<p1>()[i];
			boost::fusion::at_c<p2>(data)[i] = b.template get<p2>()[i];
		}
	}

	/*! \brief Constructor from initializer list
	 *
	 * Constructor from initializer list
	 *
	 * \param p1 Low point, initialize as a list example {0.0,0.0,0.0}
	 * \param p2 High point, initialized as a list example {1.0,1.0,1.0}
	 *
	 */

	inline void set(std::initializer_list<T> p1, std::initializer_list<T> p2)
	{
		size_t i = 0;
	    for(T x : p1)
	    {setLow(i,x);i++;}

	    i = 0;
	    for(T x : p2)
	    {setHigh(i,x);i++;}
	}

	/*! \brief set the low interval of the box
	 *
	 * set the low interval of the box
	 *
	 * \param i dimension
	 * \param val value to set
	 *
	 */
	inline void setLow(int i, T val)
	{
		boost::fusion::at_c<p1>(data)[i] = val;
	}

	/*! \brief set the high interval of the box
	 *
	 * set the high interval of the box
	 *
	 * \param i dimension
	 * \param val value to set
	 *
	 */
	inline void setHigh(int i, T val)
	{
		boost::fusion::at_c<p2>(data)[i] = val;
	}

	/*! \brief get the i-coordinate of the low bound interval of the box
	 *
	 * \param i dimension
	 *
	 * \return i-coordinate
	 *
	 */
	inline T getLow(int i) const
	{
		return boost::fusion::at_c<p1>(data)[i];
	}

	/*! \brief get the high interval of the box
	 *
	 * \param i dimension
	 * \return i coordinate of the high interval
	 *
	 */
	inline T getHigh(int i) const
	{
		return boost::fusion::at_c<p2>(data)[i];
	}

	/*! \brief Get the box enclosing this Box
	 *
	 * basically return itself
	 *
	 * \return itself
	 *
	 */
	Box<dim,T> & getBox()
	{
		return *this;
	}

	/*! \brief Get the internal boost::fusion::vector that store the data
	 *
	 * \return the internal boost::fusion::vector that store the data
	 *
	 */

	type & getVector()
	{
		return data;
	}

	/* \brief Get the key to the point 1
	 *
	 * \return the key to the point 1
	 *
	 */

	grid_key_dx<dim> getKP1() const
	{
		// grid key to return
		grid_key_dx<dim> ret(boost::fusion::at_c<p1>(data));

		return ret;
	}

	/* \brief Get the key to point 2
	 *
	 * \return the key to the point 2
	 *
	 */

	grid_key_dx<dim> getKP2() const
	{
		// grid key to return
		grid_key_dx<dim> ret(boost::fusion::at_c<p2>(data));

		return ret;
	}

	/* \brief Get the point 1
	 *
	 * \return the point 1
	 *
	 */

	inline Point<dim,T> getP1() const
	{
		// grid key to return
		Point<dim,T> ret(boost::fusion::at_c<p1>(data));

		return ret;
	}

	/* \brief Get the point 2
	 *
	 * \return the point 2
	 *
	 */

	inline Point<dim,T> getP2() const
	{
		// grid key to return
		Point<dim,T> ret(boost::fusion::at_c<p2>(data));

		return ret;
	}

	/* \brief expand expand the box by a vector
	 *
	 * \param vector
	 *
	 */
	inline void expand(T (& exp)[dim])
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			boost::fusion::at_c<p2>(data)[i] = boost::fusion::at_c<p2>(data)[i] + exp[i];
		}
	}

	/*! \brief Enlarge the box with ghost margin
	 *
	 * \param the box
	 * \param gh spacing of the margin to enlarge
	 *
	 */
	void enlarge(Box<dim,T> & gh)
	{
		typedef ::Box<dim,T> g;

		for (size_t j = 0 ; j < dim ; j++)
		{
			this->setLow(j,this->template getBase<g::p1>(j) + gh.template getBase<g::p1>(j));
			this->setHigh(j,this->template getBase<g::p2>(j) + gh.template getBase<g::p2>(j));
		}
	}

	/*! \brief Refine the box to enclose the given box and itself
	 *
	 * \param en Box to enclose
	 *
	 */
	inline void enclose(Box<dim,T> & en)
	{
		for (size_t j = 0 ; j < dim ; j++)
		{
			if (getLow(j) > en.getLow(j))
				this->setLow(j,en.getLow(j));

			if (getHigh(j) < en.getHigh(j))
				this->setHigh(j,en.getHigh(j));
		}
	}

	/*! \brief Refine the box to be contained in the given box and itself
	 *
	 * All the boxes are considered centered at p1, so it only count its relative size
	 *
	 * \param en Box to be contained
	 *
	 */
	inline void contained(Box<dim,T> & en, const bool reset_p1 = true)
	{
		for (size_t j = 0 ; j < dim ; j++)
		{
			if (getHigh(j) > (en.getHigh(j) - en.getLow(j)))
				setHigh(j,en.getHigh(j) - en.getLow(j));

			if (reset_p1 == true)
				setLow(j,0);
		}
	}
};

#endif
