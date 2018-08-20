
#ifndef BOX_HPP_
#define BOX_HPP_

#include "Space/Shape/Sphere.hpp"
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include "Grid/grid_key.hpp"
#include "Grid/Encap.hpp"
#include <sstream>

/*! \brief It define if we want the upper base or the down base (Lower or upper)
 * extreme of the interval
 *
 */

enum Base
{
	UP,
	DOWN
};

/*! \brief This class represent an N-dimensional box
 *
 * The box is defined by two points p2 and p1
 *
  \verbatim

                       +---------+ p2
                       |         |
                       |         |
                       |         |
                       |         |
                       |         |
                  p1   +---------+

  \endverbatim
 *
 * \tparam dim dimensionality of the Box
 * \tparam T type of space ... double float int size_t
 *
 * ### Expand the box with some spacing
 * \snippet Box_unit_tests.hpp expand the box with some spacing
 * ### Create an enclosing box
 * \snippet Box_unit_tests.hpp create an enclosing box
 * ### Create the smallest boxes between several boxes
 * \snippet Box_unit_tests.hpp Create the smallest boxes between several boxes
 * ### Enlarge the box
 * \snippet Box_unit_tests.hpp Enlarge the box
 * ### Enlarge the box with fixed P1
 * \snippet Box_unit_tests.hpp Enlarge the box with fixed P1
 *
 * \see SpaceBox
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

	//! It store the two point bounding the box
	type data;

	//! Low point
	static const unsigned int p1 = 0;
	//! High point
	static const unsigned int p2 = 1;
	//! Maximum number of properties
	static const unsigned int max_prop = 2;

	//! dimensionality of the box
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
	bool Intersect(const Box<dim,T> & b, Box<dim,T> & b_out) const
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
	template<typename Mem> bool Intersect(const encapc<1,Box<dim,T>,Mem> & e_b, Box<dim,T> & b_out) const
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

	template<unsigned int b> T getBase(const unsigned int i) const
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

	__device__ __host__ Box<dim,T> & operator=(const Box<dim,T> & box)
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

	/*! \brief Constructor from two points
	 *
	 * \param p1 Low point, initialize as a list example {0.0,0.0,0.0}
	 * \param p2 High point, initialized as a list example {1.0,1.0,1.0}
	 *
	 */
	Box(const Point<dim,T> & p1, const Point<dim,T> & p2)
	{
		setP1(p1);
		setP2(p2);
	}

	/*! \brief Constructor from initializer list
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
	 * \param high array indicating the coordinates of the low point
	 * \param low array indicating the coordinates of the high point
	 *
	 */
	inline Box(T * low, T * high)
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
	 * \param box from which to construct
	 *
	 */
	__device__ __host__ inline Box(const Box<dim,T> & box)
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
	 * \param box_data array from which to construct the box
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

	/*! \brief constructor from 2 grid_key_dx
	 *
	 * \param key1 start point
	 * \param key2 stop point
	 *
	 */
	inline Box(const grid_key_dx<dim> & key1, const grid_key_dx<dim> & key2)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			setLow(i,key1.get(i));
			setHigh(i,key2.get(i));
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

	/*! \brief Box constructor from encapsulated box
	 *
	 * \param b box from which to construct the vector (encapsulated)
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

	/*! \brief constructor from a Box of different type
	 *
	 * \param b box
	 *
	 */
	template <typename S> inline Box(const Box<dim,S> & b)
	{
		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,b.getLow(d));}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,b.getHigh(d));}
	}

	/*! \brief Divide component wise each box points with a point
	 *
	 * \param p point
	 *
	 * \return itself
	 *
	 */
	inline Box<dim,T> & operator/=(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			setLow(i, getLow(i)/p.get(i));
			setHigh(i, getHigh(i)/p.get(i));
		}
		return *this;
	}

	/*! \brief Multiply component wise each box points with a point
	 *
	 * \param p point
	 *
	 * \return the result box
	 *
	 */
	inline Box<dim,T> operator*(const Point<dim,T> & p)
	{
		Box<dim,T> ret;

		for (size_t i = 0 ; i < dim ; i++)
		{
			ret.setLow(i, getLow(i)*p.get(i));
			ret.setHigh(i, getHigh(i)*p.get(i));
		}
		return ret;
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
		{
			setLow(i,x);
			i++;
			if (i >= dim)
				break;
		}

		i = 0;
		for(T x : p2)
		{
			setHigh(i,x);
			i++;
			if (i >= dim)
				break;
		}
	}

	/*! \brief set the low interval of the box
	 *
	 * \param i dimension
	 * \param val value to set
	 *
	 */
	__device__ __host__  inline void setLow(int i, T val)
	{
		boost::fusion::at_c<p1>(data)[i] = val;
	}

	/*! \brief set the high interval of the box
	 *
	 * \param i dimension
	 * \param val value to set
	 *
	 */
	__device__ __host__  inline void setHigh(int i, T val)
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
	__device__ __host__  inline T getLow(int i) const
	{
		return boost::fusion::at_c<p1>(data)[i];
	}

	/*! \brief get the high interval of the box
	 *
	 * \param i dimension
	 * \return i coordinate of the high interval
	 *
	 */
	__device__ __host__ inline T getHigh(int i) const
	{
		return boost::fusion::at_c<p2>(data)[i];
	}

	/*! \brief Set the point P1 of the box
	 *
	 * \param p1 point
	 *
	 */
	inline void setP1(const grid_key_dx<dim> & p1)
	{
		for (size_t i = 0 ; i < dim ; i++)
			setLow(i,p1.get(i));
	}

	/*! \brief Set the point P2 of the box
	 *
	 * \param p2 point
	 *
	 */
	inline void setP2(const grid_key_dx<dim> & p2)
	{
		for (size_t i = 0 ; i < dim ; i++)
			setHigh(i,p2.get(i));
	}

	/*! \brief Set the point P1 of the box
	 *
	 * \param p1 point
	 *
	 */
	inline void setP1(const Point<dim,T> & p1)
	{
		for (size_t i = 0 ; i < dim ; i++)
			setLow(i,p1.get(i));
	}

	/*! \brief Set the point P2 of the box
	 *
	 * \param p2 point
	 *
	 */
	inline void setP2(const Point<dim,T> & p2)
	{
		for (size_t i = 0 ; i < dim ; i++)
			setHigh(i,p2.get(i));
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

	/*! \brief Get the box enclosing this Box
	 *
	 * basically return itself
	 *
	 * \return itself
	 *
	 */
	const Box<dim,T> & getBox() const
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

	/*! \brief Get the point p1 as grid_key_dx
	 *
	 * \return the key
	 *
	 */
	grid_key_dx<dim> getKP1() const
	{
		// grid key to return
		grid_key_dx<dim> ret(boost::fusion::at_c<p1>(data));

		return ret;
	}

	/*! \brief Get the point p12 as grid_key_dx
	 *
	 * \return the key
	 *
	 */
	grid_key_dx<dim> getKP2() const
	{
		// grid key to return
		grid_key_dx<dim> ret(boost::fusion::at_c<p2>(data));

		return ret;
	}

	/*! \brief Get the point p1
	 *
	 * \return the point p1
	 *
	 */
	inline Point<dim,T> getP1() const
	{
		// grid key to return
		Point<dim,T> ret(boost::fusion::at_c<p1>(data));

		return ret;
	}

	/*! \brief Get the point p2
	 *
	 * \return the point p2
	 *
	 */

	inline Point<dim,T> getP2() const
	{
		// grid key to return
		Point<dim,T> ret(boost::fusion::at_c<p2>(data));

		return ret;
	}

	/*! \brief Translate the box
	 *
	 * \param p Point translation vector
	 *
	 * \return itself
	 *
	 */
	inline Box<dim,T> & operator-=(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			boost::fusion::at_c<p2>(data)[i] -= p.get(i);
			boost::fusion::at_c<p1>(data)[i] -= p.get(i);
		}

		return *this;
	}

	/*! \brief Translate the box
	 *
	 * \param p Point translation vector
	 *
	 * \return itself
	 *
	 */
	inline Box<dim,T> & operator+=(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			boost::fusion::at_c<p2>(data)[i] += p.get(i);
			boost::fusion::at_c<p1>(data)[i] += p.get(i);
		}

		return *this;
	}

	/*! \brief expand the box by a vector
	 *
	 * only P2 is expanded
	 *
	 * \param exp expand vector
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
	 * \param gh spacing of the margin to enlarge
	 *
	 *
	 *\verbatim
                            ^ gh.p2[1]
                            |
                            |
                       +----+----+
                       |         |
                       |         |
         gh.p1[0]<-----+         +----> gh.p2[0]
                       |         |
                       |         |
                       +----+----+
                            |
                            v  gh.p1[1]

       \endverbatim
	 *
	 */
	void enlarge(const Box<dim,T> & gh)
	{
		typedef ::Box<dim,T> g;

		for (size_t j = 0 ; j < dim ; j++)
		{
			this->setLow(j,this->template getBase<g::p1>(j) + gh.template getBase<g::p1>(j));
			this->setHigh(j,this->template getBase<g::p2>(j) + gh.template getBase<g::p2>(j));
		}
	}

	/*! \brief Enlarge the box with ghost margin keeping fix the point P1
	 *
	 * \param gh spacing of the margin to enlarge
	 *
	 *
	 *\verbatim
                            ^ gh.p2[1]
                            |
                            |
                       +----+----+
                       |         |
                       |         |
         gh.p1[0]<-----+         +----> gh.p2[0]
                       |         |
                       |         |
                       +----+----+
                            |
                            v  gh.p1[1]

       \endverbatim
	 *
	 */
	template<typename S> inline void enlarge_fix_P1(Box<dim,S> & gh)
	{
		typedef ::Box<dim,T> g;

		for (size_t j = 0 ; j < dim ; j++)
		{
			this->setHigh(j,this->template getBase<g::p2>(j) + gh.template getBase<g::p2>(j) - gh.template getBase<g::p1>(j));
		}
	}

	/*! \brief Invalidate the box
	 *
	 * Bring the state of this box in a way that isValid return false
	 *
	 */
	void invalidate()
	{
		for (size_t j = 0 ; j < dim ; j++)
		{
			this->setLow(j,this->getHigh(j)+1);
		}
	}


	/*! \brief Magnify the box
	 *
	 * For example 1.001 enlarge the box of 0.1% on each direction
	 *
	 * \warning P1 is mooved if not zero
	 *
	 * \param mg Magnification factor
	 *
	 */
	void magnify(T mg)
	{
		typedef ::Box<dim,T> g;

		for (size_t j = 0 ; j < dim ; j++)
		{
			this->setLow(j,mg * this->template getBase<g::p1>(j));
			this->setHigh(j,mg * this->template getBase<g::p2>(j));
		}
	}

	/*! \brief Magnify the box by a factor keeping fix the point P1
	 *
	 * For example 1.001 enlarge the box of 0.1% on each direction
	 *
	 * \param mg Magnification factor
	 *
	 */
	inline void magnify_fix_P1(T mg)
	{
		typedef ::Box<dim,T> g;

		for (size_t j = 0 ; j < dim ; j++)
		{
			this->setHigh(j,this->template getBase<g::p1>(j) + mg * (this->template getBase<g::p2>(j) - this->template getBase<g::p1>(j)));
		}
	}

	/*! \brief Shrink moving p2 of sh quantity (on each direction)
	 *
	 * \param sh
	 *
	 */
	inline void shrinkP2(T sh)
	{
		for (size_t j = 0 ; j < dim ; j++)
		{
			this->setHigh(j,this->getHigh(j) - sh);
		}
	}

	/*! \brief Refine the box to enclose the given box and itself
	 *
	 * \param en Box to enclose
	 *
	 */
	inline void enclose(const Box<dim,T> & en)
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
	 * \param reset_p1 if true set p1 to 0
	 *
	 */
	inline void contained(const Box<dim,T> & en, const bool reset_p1 = true)
	{
		for (size_t j = 0 ; j < dim ; j++)
		{
			if (getHigh(j) > (en.getHigh(j) - en.getLow(j)))
				setHigh(j,en.getHigh(j) - en.getLow(j));

			if (reset_p1 == true)
				setLow(j,0);
		}
	}

	/*! \brief Set p1 and p2 to 0
	 *
	 */
	inline void zero()
	{
		for (size_t j = 0 ; j < dim ; j++)
		{
			setHigh(j,0);
			setLow(j,0);
		}
	}

	/*! \brief Check if the box is contained
	 *
	 * \param b Box
	 *
	 * \return true if the box is contained
	 *
	 */
	inline bool isContained(const Box<dim,T> & b) const
	{
		bool isc = true;

		isc &= isInside(b.getP1());
		isc &= isInside(b.getP2());

		return isc;
	}

	/*! \brief Check if the point is inside the box
	 *
	 * \param p point to check
	 * \return true if the point is inside the space
	 *
	 */

	inline bool isInside(const Point<dim,T> & p) const
	{
		// check if bound

		for (size_t i = 0 ; i < dim ; i++)
		{
			// if outside the region return false
			if (   boost::fusion::at_c<Point<dim,T>::x>(p.data)[i] < boost::fusion::at_c<Box<dim,T>::p1>(this->data)[i]
			    || boost::fusion::at_c<Point<dim,T>::x>(p.data)[i] > boost::fusion::at_c<Box<dim,T>::p2>(this->data)[i])
			{
				// Out of bound

				return false;
			}

		}

		// In bound

		return true;
	}

	/*! \brief Check if the point is inside the region excluding the positive part
	 *
	 * In periodic boundary conditions the positive border is not included, but match the beginning
	 *
	 * \param p point to check
	 * \return true if the point is inside the space
	 *
	 */
	__device__ __host__ inline bool isInsideNP(const Point<dim,T> & p) const
	{
		// check if bound

		for (size_t i = 0 ; i < dim ; i++)
		{
			// if outside the region return false
			if (   boost::fusion::at_c<Point<dim,T>::x>(p.data)[i] < boost::fusion::at_c<Box<dim,T>::p1>(this->data)[i]
			    || boost::fusion::at_c<Point<dim,T>::x>(p.data)[i] >= boost::fusion::at_c<Box<dim,T>::p2>(this->data)[i])
			{
				// Out of bound



				return false;
			}

		}

		// In bound

		return true;
	}

	/*! \brief Check if the point is inside the region excluding the borders
	 *
	 * \param p point to check
	 * \return true if the point is inside the space
	 *
	 */
	inline bool isInsideNB(const Point<dim,T> & p) const
	{
		// check if bound

		for (size_t i = 0 ; i < dim ; i++)
		{
			// if outside the region return false
			if (   boost::fusion::at_c<Point<dim,T>::x>(p.data)[i] <= boost::fusion::at_c<Box<dim,T>::p1>(this->data)[i]
			    || boost::fusion::at_c<Point<dim,T>::x>(p.data)[i] >= boost::fusion::at_c<Box<dim,T>::p2>(this->data)[i])
			{
				// Out of bound

				return false;
			}

		}

		// In bound

		return true;
	}

	/*! \brief Check if the point is inside the region
	 *
	 * \param p point to check
	 * \return true if the point is inside the space
	 *
	 */

	inline bool isInside(const T (&p)[dim]) const
	{
		// check if bound

		for (size_t i = 0 ; i < dim ; i++)
		{
			// if outside the region return false
			if (   p[i] < boost::fusion::at_c<Box<dim,T>::p1>(this->data)[i]
			    || p[i] > boost::fusion::at_c<Box<dim,T>::p2>(this->data)[i])
			{
				// Out of bound

				return false;
			}

		}

		// In bound

		return true;
	}


	/*! \brief Check if the Box is a valid box P2 >= P1
	 *
	 * \return true if it is valid
	 *
	 */
	inline bool isValid() const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (getLow(i) > getHigh(i))
				return false;
		}

		return true;
	}

	/*! \brief Check if the Box is a valid box P2 > P1
	 *
	 * \return true if it is valid
	 *
	 */
	inline bool isValidN() const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (getLow(i) >= getHigh(i))
				return false;
		}

		return true;
	}

	/*! \brief Apply the ceil operation to the point P1
	 *
	 *
	 */
	inline void floorP1()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			setLow(i,std::floor(getLow(i)));
		}
	}

	/*! \brief Apply the ceil operation to the point P2
	 *
	 *
	 */
	inline void floorP2()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			setHigh(i,std::floor(getHigh(i)));
		}
	}

	/*! \brief Apply the ceil operation to the point P1
	 *
	 *
	 */
	inline void ceilP1()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			setLow(i,std::ceil(getLow(i)));
		}
	}

	/*! \brief Apply the ceil operation to the point P2
	 *
	 *
	 */
	inline void ceilP2()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			setHigh(i,std::ceil(getHigh(i)));
		}
	}

	/*! \brief Shrink the point P2 by one vector
	 *
	 * \param p vector
	 *
	 */
	inline void shrinkP2(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			setHigh(i,getHigh(i) - p.get(i));
		}
	}

	/*! \brief exchange the data of two boxes
	 *
	 * \param b box to switch
	 *
	 */
	void swap(Box<dim,T> & b)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			T tmp_l = getLow(i);
			T tmp_h = getHigh(i);

			setLow(i,b.getLow(i));
			setHigh(i,b.getHigh(i));

			b.setLow(i,tmp_l);
			b.setHigh(i,tmp_h);
		}
	}

	/*! \brief Check if the point is inside the region
	 *
	 * \param p point to check
	 * \return true if the point is inside the space
	 *
	 */

	template <typename Mem> inline bool isInside(const encapc<1,Point<dim,T>,Mem> & p)
	{
		// check if bound

		for (size_t i = 0 ; i < dim ; i++)
		{
			// if outside the region return false
			if (   p.template get<Point<dim,T>::x>()[i] < boost::fusion::at_c<Box<dim,T>::p1>(this->data)[i]
			    || p.template get<Point<dim,T>::x>()[i] > boost::fusion::at_c<Box<dim,T>::p2>(this->data)[i])
			{
				// Out of bound

				return false;
			}

		}

		// In bound

		return true;
	}

	/*! \brief Get the volume of the box
	 *
	 * \return the box volume
	 *
	 */
	inline T getVolume() const
	{
		T vol = 1.0;

		for (size_t i = 0 ; i < dim ; i++)
			vol *= (getHigh(i) - getLow(i));

		return vol;
	}

	/*! \brief Get the volume spanned by the Box P1 and P2 interpreted as grid key
	 *
	 * \return The volume
	 *
	 */
	inline T getVolumeKey() const
	{
		T vol = 1.0;

		for (size_t i = 0 ; i < dim ; i++)
			vol *= (getHigh(i) - getLow(i) + 1.0);

		return vol;
	}

	/*! \brief Get the volume spanned by the Box as grid_key_dx_iterator_sub
	 *
	 * \warning Careful it is not the simple volume calculation there is a +1 on each dimension, consider the case of a subgrid iterator with
	 *          P1 = {5,7} ; P2  = {5,7}, the sub-grid iterator has one point {5,7}, that mean Volume=1, so
	 *          the volume formula is (5 - 5 + 1) * (7 - 7 + 1)
	 *
	 * \param p1 point p1
	 * \param p2 point p2
	 *
	 * \return The volume
	 *
	 */
	inline static T getVolumeKey(const T (&p1)[dim], const T(&p2)[dim])
	{
		T vol = 1.0;

		for (size_t i = 0 ; i < dim ; i++)
			vol *= (p2[i] - p1[i] + 1.0);

		return vol;
	}

	//! This structure has no internal pointers
	static bool noPointers()
	{
		return true;
	}

	/*! \brief Return the middle point of the box
	 *
	 * \return the middle point of the box
	 *
	 */
	inline Point<dim,T> middle() const
	{
		Point<dim,T> p;

		for (size_t i = 0 ; i < dim ; i++)
			p.get(i) = (getLow(i) + getHigh(i))/2;

		return p;
	}

	/*! \brief Produce a string from the object
	 *
	 * \return string
	 *
	 */
	std::string toString() const
	{
		std::stringstream str;

		for (size_t i = 0 ; i < dim ; i++)
			str << "x[" << i << "]=" << getLow(i) << " ";

		str << "   |  ";

		for (size_t i = 0 ; i < dim ; i++)
			str << "x[" << i << "]=" << getHigh(i) << " ";

		return str.str();
	}

	/*! \brief Compare two boxes
	 *
	 * \param b
	 *
	 * \return true if the boxes are equal
	 *
	 */
	bool operator==(const Box<dim,T> & b) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (getLow(i) != b.getLow(i))
				return false;

			if (getHigh(i) != b.getHigh(i))
				return false;
		}

		return true;
	}


	/*! \brief Compare two boxes
	 *
	 * \param b
	 *
	 * \return true if the boxes are equal
	 *
	 */
	bool operator!=(const Box<dim,T> & b) const
	{
		return ! this->operator==(b);
	}
};

#endif
