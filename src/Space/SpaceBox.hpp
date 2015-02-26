
#ifndef SPACEBOX_HPP_
#define SPACEBOX_HPP_

#include "Shape/Point.hpp"
#include "Shape/Box.hpp"
#include <boost/fusion/include/vector.hpp>
#include "Grid/Encap.hpp"

/** \brief This class represent an N-dimensional box
 *
 * This class represent an N-dimensional box embedded in an N dimensional space
 *
 * \param T type of space ... Real Complex Integer
 * \param N dimensionality of the Box
 *
 */

template<unsigned int dim, typename T>
class SpaceBox : public Box<dim,T>
{
	public:

	//! layout that interleave the properties
	typedef typename Box<dim,T>::memory_int memory_int;
	//! layout with linear properties
	typedef typename Box<dim,T>::memory_lin memory_lin;

	/*! \brief Check if the point is inside the region
	 *
	 * \param p point to check
	 * \return true if the point is inside the space
	 *
	 */

	bool isBound(Point<dim,T> p)
	{
		// check if bound

		for (int i = 0 ; i < dim ; i++)
		{
			// if outside the region return false
			if (   boost::fusion::at_c<Point<dim,T>::x>(p.data)[i] < boost::fusion::at_c<Box<dim,T>::p1>(this->data)[i]
			    && boost::fusion::at_c<Point<dim,T>::x>(p.data)[i] < boost::fusion::at_c<Box<dim,T>::p2>(this->data)[i])
			{
				// Out of bound

				return false;
			}

		}

		// In bound

		return true;
	}

	/*! \brief Define the box from a box shape
	 *
	 * Define the box from a box shape
	 *
	 * \param b is the box
	 * \return itself
	 *
	 */

	SpaceBox<dim,T> & operator=(Box<dim,T> & b)
	{
		// Copy the element of the box to this box

		this->data = b.data;

		return *this;
	}

	/*! \brief constructor from a box
	 *
	 * constructor from a box
	 *
	 * \param b is the box
	 *
	 */

	SpaceBox(Box<dim,T> & b)
	{
		// for each dimension set high and low

		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,b.getLow(d));}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,b.getHigh(d));}
	}

	/*! \brief Constructor from a Box
	 *
	 * \param Box
	 *
	 */

	template<unsigned int dim_s,typename Mem, typename S>SpaceBox(const encapc<dim_s,Box<dim,S>,Mem> & box)
	{
		// for each dimension set high and low

		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,box.template get<Box<dim,S>::p1>()[d]);}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,box.template get<Box<dim,S>::p2>()[d]);}
	}

	/*! \brief Re-scale the space box with the coefficient defined in sp
	 *
	 * \param sp
	 *
	 */

	void mul(float (& sp)[dim])
	{
		for (int i = 0  ; i < dim ; i++)
		{
			for (size_t d = 0 ; d < dim ; d++)
			{this->setLow(d,this->getLow(d) * sp[i]);}

			for (size_t d = 0 ; d < dim ; d++)
			{this->setHigh(d,this->getHigh(d) * sp[i]);}
		}
	}

	/*! \brief Re-scale the space box with the coefficient defined in sp
	 *
	 * \param sp
	 *
	 */

	void mul(size_t (& sp)[dim])
	{
		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,this->getLow(d) * sp[d]);}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,this->getHigh(d) * sp[d]);}
	}

	/*! \brief Re-scale the space box with the spacing defined in sp
	 *
	 * \param sp spacing
	 * \param s_sp sub-domain size
	 *
	 */

	void spacing(float (& sp)[dim])
	{
		mul(sp);
	}

	/*! \brief Re-scale the space box with the spacing defined in sp
	 *
	 * \param sp spacing
	 *
	 */

	void spacing(size_t (& sp)[dim])
	{
		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,this->getLow(d) * sp[d]);}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,((this->getHigh(d)+1) * sp[d])-1);}
	}

	//! Default constructor
	SpaceBox<dim,T>()	{}
};

#endif
