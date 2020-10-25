#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include "boost/multi_array.hpp"
#include "Point.hpp"
#include <Space/Shape/Point.hpp>

/*! \brief This class implement the Sphere concept in an N-dimensional space
 *
 * This class implement the sphere shape in an N-dimensional space
 *
 * \tparam T type of the space
 * \tparam dim dimensionality
 *
 */

template<unsigned int dim ,typename T> class Sphere
{
	public:

	//! boost fusion that store the point
	typedef boost::fusion::vector<T[dim],T> type;

	//! Structure that store the data
	type data;

	//! property id of the center position of the sphere
	static const unsigned int x = 0;
	//! property id of the radius of the sphere
	static const unsigned int r = 1;

	/*! \brief Get the component i of the center
	 *
	 * \param i coordinate
	 * \return the coordinate i of the center
	 *
	 */
	__device__ __host__ T center(unsigned int i)
	{
		return boost::fusion::at_c<x>(data)[i];
	}

	/*! \brief Sphere constructor
	 *
	 * Sphere constructor
	 *
	 * \tparam k dimensionality of the center point (and the sphere)
	 * \param c center point
	 * \param radius
	 *
	 */
	__device__ __host__ Sphere(const Sphere<dim,T> & sph)
	{
		// Copy the center
		for (int i = 0 ;  i < dim ; i++)
		{
			boost::fusion::at_c<x>(data)[i] = boost::fusion::at_c<x>(sph.data)[i];
		}

		boost::fusion::at_c<r>(data) = boost::fusion::at_c<r>(sph.data);
	}

	/*! \brief Sphere constructor
	 *
	 * Sphere constructor
	 *
	 * \tparam k dimensionality of the center point (and the sphere)
	 * \param c center point
	 * \param radius
	 *
	 */
	template<unsigned int k>
	Sphere(boost::fusion::vector<T[k]> & c, T radius)
	{
		// Copy the center
		for (int i = 0 ;  i < dim ; i++)
		{
			boost::fusion::at_c<x>(data)[i] = boost::fusion::at_c<x>(c)[i];
		}

		boost::fusion::at_c<r>(data) = radius;
	}

	/*! \brief Sphere constructor
	 *
	 * Sphere constructor
	 *
	 * \tparam k dimensionality of the center point (and the sphere)
	 * \param c center point
	 * \param radius
	 *
	 */
	Sphere(Point<dim,double> & c, T radius)
	{
		// Copy the center
		for (int i = 0 ;  i < dim ; i++)
		{
			boost::fusion::at_c<x>(data)[i] = c.get(i);
		}

		boost::fusion::at_c<r>(data) = radius;
	}

	/*! \brief Get the radius of the sphere
	 *
	 * \return the radius of the sphere
	 *
	 */
	__device__ __host__ T radius() const
	{
		return boost::fusion::at_c<r>(data);
	}

	/*! \brief Check if a point is inside
	 *
	 * \tparam distance distance functor used to calculate the distance two points
	 * \param p Point to check if it is inside
	 *
	 * \return true if the point is inside
	 *
	 */
	__device__ __host__ bool isInside(Point<dim,T> p) const
	{
		T dist;

		// calculate the distance of the center from the point

		for (int i = 0; i < dim ; i++)
		{
			dist += (boost::fusion::at_c<x>(data)[i] - p.get(i) )*(boost::fusion::at_c<x>(data)[i] - p.get(i) );
		}

		// Check if the distance is smaller than the radius

		if (dist <= boost::fusion::at_c<r>(data)*boost::fusion::at_c<r>(data))
		{return true;}


		return false;
	}

	/*! \brief Check if a point is inside
	 *
	 * \tparam distance distance functor used to calculate the distance two points
	 * \param p Point to check if it is inside
	 *
	 * \return true if the point is inside
	 *
	 */
	template<typename Distance>
	__device__ __host__ bool isInside(Point<dim,T> p) const
	{
		T dist;

		// Object to calculate distances
		Distance d;

		// calculate the distance of the center from the point

		for (int i = 0; i < dim ; i++)
		{
			dist += d.accum_dist(boost::fusion::at_c<x>(data)[i],p.get(i) );
		}

		// Check if the distance is smaller than the radius

		if (dist <= boost::fusion::at_c<r>(data))
		{return true;}


		return false;
	}

	/*! \brief Is the point inside the sphere
	 *
	 * \param pnt check if the point is inside
	 * \return true if it is inside
	 *
	 */

	template<typename Distance> bool
	__device__ __host__ isInside(float * pnt) const
	{
		T dist;

		// Object to calculate distances
		Distance d;

		// calculate the distance of the center from the point

		for (int i = 0; i < dim ; i++)
		{
			dist += d.accum_dist(boost::fusion::at_c<x>(data)[i],pnt[i],i);
		}

		// Check if the distance is smaller than the radius

		if (dist <= boost::fusion::at_c<r>(data))
		{return true;}

		return false;
	}

	/*! \brief Return the distance from the surface
	 *
	 * \param p point
	 * \return the distance from the surface. The sign indicate if is outside or inside the shape
	 *
	 */
	__device__ __host__ T distance(Point<dim,T> & p) const
	{
		T dist = 0.0;

		// calculate the distance of the center from the point

		for (int i = 0; i < dim ; i++)
		{
			dist += (boost::fusion::at_c<x>(data)[i] - p.get(i))*(boost::fusion::at_c<x>(data)[i] - p.get(i));
		}

		return sqrt(dist) - radius();
	}
};

#endif
