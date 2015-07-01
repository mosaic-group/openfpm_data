
#ifndef SPACEBOX_HPP_
#define SPACEBOX_HPP_

#include "Shape/Point.hpp"
#include "Shape/Box.hpp"
#include <boost/fusion/include/vector.hpp>
#include "Grid/Encap.hpp"
#include "Ghost.hpp"
#include <stdlib.h>     /* srand, rand */

/** \brief This class represent an N-dimensional box
 *
 * \param T type of space ... double float int size_t
 * \param N dimensionality of the Box
 *
 * ### Definition of a spacebox and rescale
 * \snippet SpaceBox_unit_tests.hpp Definition of a spacebox and rescale
 * ### Definition of a spaceboxes and intersection between them
 * \snippet SpaceBox_unit_tests.hpp Definition of a spacebox and intersection between them
 * ### Create random points inside the SpaceBox
 * \snippet SpaceBox_unit_tests.hpp Create random points inside the SpaceBox
 *
 */
template<unsigned int dim, typename T>
class SpaceBox : public Box<dim,T>
{
	public:

	/*! \brief Define the box from a box shape
	 *
	 * Define the box from a box shape
	 *
	 * \param b is the box
	 * \return itself
	 *
	 */

	inline SpaceBox<dim,T> & operator=(const Box<dim,T> & b)
	{
		// for each dimension set high and low

		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,b.getLow(d));}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,b.getHigh(d));}

		return *this;
	}

	/*! \brief constructor from a Box of different type
	 *
	 * \param b box
	 *
	 */
	template <typename S> inline SpaceBox(const Box<dim,S> & b)
	{
		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,b.getLow(d));}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,b.getHigh(d));}
	}

	/*! \brief constructor from a SpaceBox
	 *
	 * constructor from a SpaceBox
	 *
	 * \param b is the SpaceBox
	 *
	 */
	SpaceBox(const SpaceBox<dim,T> & b)
	:Box<dim,T>(b)
	{
	}

	/*! \brief constructor from a box
	 *
	 * constructor from a box
	 *
	 * \param b is the box
	 *
	 */

	SpaceBox(const Box<dim,T> & b)
	:Box<dim,T>(b)
	{
	}

	/*! \brief Constructor from a Box
	 *
	 * \param box Box (Encapsulated)
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

	/*! \brief Constructor from a Box
	 *
	 * \param box box (Encapsulated)
	 *
	 */

	template<unsigned int dim_s,typename Mem, typename S>SpaceBox(const encapc<dim_s,SpaceBox<dim,S>,Mem> & box)
	{
		// for each dimension set high and low

		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,box.template get<Box<dim,S>::p1>()[d]);}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,box.template get<Box<dim,S>::p2>()[d]);}
	}

	/*! \brief Constructor from initializer list
	 *
	 * Constructor from initializer list
	 *
	 * \param p1 Low point, initialize as a list example {0.0,0.0,0.0}
	 * \param p2 High point, initialized as a list example {1.0,1.0,1.0}
	 *
	 */

	SpaceBox(std::initializer_list<T> p1, std::initializer_list<T> p2)
	{
		// for each dimension set high and low

		size_t i = 0;
	    for(T x : p1)
	    {this->setLow(i,x);i++;}

	    i = 0;
	    for(T x : p2)
	    {this->setHigh(i,x);i++;}
	}

	/*! \brief Re-scale the space box with the coefficient defined in sp
	 *
	 * \param sp
	 *
	 */

	void rescale(float (& sp)[dim])
	{
		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,this->getLow(d) + (this->getHigh(d) -this->getLow(d)) * sp[d]);}
	}

	/*! \brief Re-scale the space box with the coefficient defined in sp
	 *
	 * \param sp
	 *
	 */

	void rescale(size_t (& sp)[dim])
	{
		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,this->getLow(d) + (this->getHigh(d) -this->getLow(d)) * sp[d]);}
	}

	/*! \brief multiply the space box with the coefficient defined in sp
	 *
	 * It rescale the domain where the space box live
	 *
	 * \param sp coefficents
	 *
	 */

	void mul(float (& sp)[dim])
	{
		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,this->getLow(d) * sp[d]);}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,this->getHigh(d) * sp[d]);}
	}

	/*! \brief multiply the space box with the coefficient defined in sp
	 *
	 * It rescale the domain where the space box live
	 *
	 * \param sp coefficents
	 *
	 */

	void mul(size_t (& sp)[dim])
	{
		for (size_t d = 0 ; d < dim ; d++)
		{this->setLow(d,this->getLow(d) * sp[d]);}

		for (size_t d = 0 ; d < dim ; d++)
		{this->setHigh(d,this->getHigh(d) * sp[d]);}
	}

	/*! \brief Generate a random point inside the box
	 *
	 * \return a random point inside the box
	 *
	 */
	Point<dim,T> rnd()
	{
		Point<dim,T> p;

		for (size_t i = 0 ; i < dim ; i++)
			p.get(i) = ((T)rand())/RAND_MAX * (this->getHigh(i) - this->getLow(i)) + this->getLow(i);

		return p;
	}

	//! Default constructor
	SpaceBox<dim,T>()	{}
};

#include "memory_c.hpp"

/*! \brief It make explicit the inheritance of SpaceBox to Box
 * for encap
 *
 * \param dim Dimensionality of the grid
 * \param T type of object the grid store
 * \param Mem suppose to be a boost::fusion::vector of arrays
 *
 */

template<unsigned int dim,typename T,typename Mem>
class encapc<dim,SpaceBox<dim,T>,Mem> : encapc<dim,Box<dim,T>,Mem>
{
};

#endif
