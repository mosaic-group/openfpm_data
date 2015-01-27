#ifndef ADAPTIVE_CONE_HPP_
#define ADAPTIVE_CONE_HPP_

#include <boost/fusion/include/mpl.hpp>
#include "Space/Shape/Box.hpp"
#include "Space/Shape/Point.hpp"
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>

/*! \brief this class is a functor for "for_each" algorithm
 *
 * For each element of the boost::vector the operator() is called to copy
 * the buffer containing the component of the point to the center
 *  point of the adaptive component
 *
 * \param S type of the data to copy (if is a float point is float, int point is
 *          int ... )
 *
 */

template<typename S>
struct copy_acc
{
	//! Pointer storing the data point
	const S * ptr;

	//! constructor it fix the size
	copy_acc(const S * ptr)
	:ptr(ptr){};

	//! It call the copy function for each member
    template<typename T>
    void operator()(T& t) const
    {
    	t = ptr[t];
    }
};

/** \brief This class represent an Adaptive cylinder cone
 *
 * This class represent an N-dimensional adaptive cylinder cone
 *
 * \param T type of space ... Real Complex Integer
 * \param N dimensionality of the adaptive cylinder cone
 *
 * This shape is useful for adaptive particle search, in order to explain how is defined we take a
 * 1D space where particles can interact each other with different radius
 *
 * Given one particle we want to know all the neighborhood particles, all the neighborhood particles
 * are given by all the particles included in the radius of interaction of the particle, plus all
 * the other particles with radius big enough to interact with this particle or in other word all
 * the particles that has this particle as neighborhood
 *
 * If we see the problem to find the neighborhood of the 1D particles in a 2D space (x,r)
 * with x=position and r=interaction radius of the particles.
 *
 * all the particles neighborhood of one particle (x_p,r_p) are all the particles in the rectangle
 * (x-r_p,0) and (x+r_p,r_p)
 *
 * in the area r < r_p
 *
 * and the cone defined by the two bases
 *
 * (x-r_p,r_p) - (x+r_p,r_p) and (x-r_max,r_max) - (x+r_max,r_max)
 *
 * This mean that one point (x_p,r_p) and r_max define the shape
 *
 */

template<unsigned int dim , typename T>
class AdaptiveCylinderCone
{
	//! Base structure that define a Point
	typedef Point<dim,T> PointA;

	public:

	//! boost fusion that store the point
	typedef boost::fusion::vector<T[dim]> type;

	//! Center point of the cylynder cone (Where the cylinder finish), where
	//! the cone start
	type data;

	//! Point where cone and cylinder base match
	static const unsigned int p1 = 0;

	/*! \brief This function check if a box intersect this shape
	 *
	 * \param b is the box to check
	 * \return true if the box intersect the cylinder cone shape
	 *
	 */
	template <typename distance>  bool Intersect(Box<dim,T> b)
	{
//		std::cout << "Intersection between" << "\n";
/*		for (int i = 0 ; i < dim ; i++)
		{
				std::cout << "Box: " << boost::fusion::at_c<0>(b.data)[i] << "     " << boost::fusion::at_c<1>(b.data)[i] << "\n";
		}*/

		// We have to check if the upper-base of the box intersect the shape, if not we have done

		// we get the upper base that is a box of dimension dim-1
		Box<dim-1,T> upBase = b.getSubBox();

		// Get the up level
		T up = b.template getBase<Base::UP>(dim-1);

		// we get the circle of dimension n-1 at the same level of the upper base of the box

		// this is the radius of the circle
		T r_p = 0;

		// if the plane is at lower level than the point r_p is the cylinder radius
		// otherwise is the cone cut

		if (b.template getBase<Base::UP>(dim-1) < boost::fusion::at_c<PointA::x>(data)[dim-1] )
		{
			r_p = boost::fusion::at_c<PointA::x>(data)[dim-1];
		}
		else
		{
			// r_p is the cone cut (plane of base box cutting the cone)

			r_p = b.template getBase<Base::UP>(dim-1);
		}

		// Create the cone-cylinder cut with the plane of the base
		Sphere<dim-1,T> acon(data,up);

		return upBase.Intersect<distance>(acon);
	}

	/*! \brief This function check if a box intersect this shape
	 *
	 * It differer from the previous one because it try to optimize the
	 * intersection calculation from the previous box intersection,
	 * usefull in KDTree, where the child are divisions of the root
	 *
	 * \tparam distance functor to calculate the distance between two points
	 * \param b the bounding Box
	 * \param[in] div where we are cutting from the previous box
	 * \param[in] dists array of the previous distances calculated on each dimension
	 * \param[in] minsqdistance minimum square distance
	 *
	 * \return true if the Box intersect the adaptive cylinder cone
	 *
	 */
	template <typename distance>  bool Intersect(Box<dim,T> b, int div, T * dists, T & minsqdistance)
	{
		// we simply ignore all the optimization parameters we just check that
		// the box intersect the adaptive cone

		return Intersect<distance>(b);
	}


	/*! \brief Create an adaptive cylinder Cone for the point pt
	 *
	 * \param pt array storing the multidimensional point
	 *
	 */
	AdaptiveCylinderCone(const T * pt)
	{
		// for all components from the pt and fill the center point of the
		// Adaptive Cylinder Cone

		for (int i = 0 ; i < dim ; i++)
		{
			boost::fusion::at_c<0>(data)[i] = pt[i];
		}
	}

	/*! \brief Is a point inside the adaptive cone
	 *
	 *  \param pnt Point to check
	 *
	 *  \return true if the point is inside, false otherwise
	 *
	 */

	template<typename Distance> bool isInside(Point<dim,T> pnt)
	{
		// Check if the point is inside the Adaptive cone

		// get the radius

		T r = pnt.get(dim-1);

		// create a sphere

		Sphere<dim-1,T> sphere(data,r);

		sphere.isInside<Distance>(pnt.getSubPoint(dim-1));
	}

	/*! \brief Is a point inside the adaptive cone
	 *
	 *  \param pnt Point to check
	 *
	 *  \return true if the point is inside false otherwise
	 *
	 */
	template<typename Distance> bool isInside(float * pnt)
	{
		// Check if the point is inside the Adaptive cone

		// get the radius

		T r = pnt[dim-1];

		// create a sphere

		Sphere<dim-1,T> sphere(data,r);

		sphere.isInside<Distance>(pnt);
	}
};

#endif
