/*
 * Ghost.hpp
 *
 *  Created on: Apr 28, 2015
 *      Author: Pietro Incardona
 */

#ifndef GHOST_HPP_
#define GHOST_HPP_

#include "SpaceBox.hpp"

/*! Ghost
 *
 * it indicate the ghost extension
 *
 * Ghost margins for each dimensions (p1 negative part) (p2 positive part)
                ^ p2[1]
                |
                |
           +----+----+
           |         |
           |         |
p1[0]<-----+         +----> p2[0]
           |         |
           |         |
           +----+----+
                |
                v  p1[1]


    \warning p1[0] and p1[1] must be negative hile p2[1] and p2[0] must be positive
 *
 */

template<unsigned int dim, typename T>
class Ghost : public Box<dim,T>
{
public:

	/*! constructor from another Ghost
	 *
	 * \param g ghost
	 *
	 */
	template <typename S> inline Ghost(const Ghost<dim,S> & g)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			this->setLow(i,g.getLow(i));
			this->setHigh(i,g.getHigh(i));
		}
	}

	// construct a ghost based on interaction radius
	inline Ghost(T r)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			this->setLow(i,-r);
			this->setHigh(i,r);
		}
	}

	// Basic constructor
	inline Ghost()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			this->setLow(i,0);
			this->setHigh(i,0);
		}
	}

	/*! \brief Divide component wise the ghost box with a point
	 *
	 * \param p point
	 *
	 * \return itself
	 *
	 */
	inline Ghost<dim,T> & operator/=(const Point<dim,T> & p)
	{
		Box<dim,T>::operator/=(p);

		return *this;
	}

};

/*! \brief Class that contain Padding information on each direction positive and Negative direction
 *
 * It is equivalent to a Ghost<dim,size_t>
 *
 * \see Ghost
 *
 */
template<unsigned int dim>
class Padding : public Ghost<dim,long int>
{
public:
	/*! \brief Constructor from initializer list
	 *
	 * \param p1 Padding left, initialize as a list example {0.0,0.0,0.0}
	 * \param p2 Padding right, initialized as a list example {1.0,1.0,1.0}
	 *
	 */
	Padding(std::initializer_list<long int> p1, std::initializer_list<long int> p2)
	{
		Box<dim,long int>::set(p1,p2);
	}
};

#endif /* GHOST_HPP_ */
