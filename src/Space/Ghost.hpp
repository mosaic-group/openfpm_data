/*
 * Ghost.hpp
 *
 *  Created on: Apr 28, 2015
 *      Author: Pietro Incardona
 */

#ifndef GHOST_HPP_
#define GHOST_HPP_

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

	// construct a ghost based on interaction radius
	Ghost(T r)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			this->setLow(i,-r);
			this->setHigh(i,r);
		}
	}
};


#endif /* GHOST_HPP_ */
