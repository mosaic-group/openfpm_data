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


#endif /* GHOST_HPP_ */
