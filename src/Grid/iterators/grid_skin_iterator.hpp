/*
 * grid_skin_iterator.hpp
 *
 *  Created on: Jun 24, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_SKIN_ITERATOR_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_SKIN_ITERATOR_HPP_

#include "grid_key_dx_iterator_sub_bc.hpp"

/**
 *
 * Grid key class iterator, iterate through the grid elements on the skin on the grid
 *
 * In particular given 2 boxes A and B
 *
 *
 \verbatim

+----------------------------+
|       skin    B            |
|   +----------------------+ |
|   |                      | |
| s |                      | |
| k |                      | |
| i |                      | |
| n |          A           | |
|   |                      | |
|   |                      | |
|   |                      | |
|   |                      | |
|   |                      | |
|   +----------------------+ |
|            skin            |
+----------------------------+

\endverbatim

 *
 * The skin is the part in between the two boxes A and B
 *
 * More in general is the Box B removed of A without border. This iterator respect the boundary
 * conditions
 *
 * \param dim dimensionality of the grid
 *
 */

template<unsigned int dim>
class grid_skin_iterator_bc
{
	grid_key_dx_iterator_sub_bc<dim> sub_it[2*dim];

	//! Actual iterator
	size_t act;

public:

	/*! \brief Constructor require a grid_sm<dim,T>
	 *
	 * \param g_sm grid information
	 * \param A box A
	 * \param B box B
	 * \param bc boundary conditions
	 *
	 */
	template <typename T> grid_skin_iterator_bc(grid_sm<dim,T> & g_sm, const Box<dim,size_t> & A, const Box<dim,size_t> & B, const size_t (& bc)[dim])
	:act(0)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			grid_key_dx<dim> k1;
			grid_key_dx<dim> k2;

			for (size_t j = 0 ; j < dim ; j++)
			{
				if (j == i)
				{
					k1.set_d(j,A.getHigh(j));
					k2.set_d(j,B.getHigh(j));
				}
				else
				{
					if ( j < i )
					{
						k1.set_d(j,A.getLow(j)+1);
						k2.set_d(j,A.getHigh(j)-1);
					}
					else
					{
						k1.set_d(j,B.getLow(j));
						k2.set_d(j,B.getHigh(j));
					}
				}
			}

			// Initialize a box from the keys and check if it is a valid box
			Box<dim,long int> br(k1,k2);
			if (br.isValid() == true)
				sub_it[2*i].Initialize(g_sm,k1,k2,bc);

			for (size_t j = 0 ; j < dim ; j++)
			{
				if (j == i)
				{
					k1.set_d(j,B.getLow(j));
					k2.set_d(j,A.getLow(j));
				}
				else
				{
					if ( j < i )
					{
						k1.set_d(j,A.getLow(j)+1);
						k2.set_d(j,A.getHigh(j)-1);
					}
					else
					{
						k1.set_d(j,B.getLow(j));
						k2.set_d(j,B.getHigh(j));
					}
				}
			}

			// Initialize a box from the keys and check if it is a valid box
			Box<dim,long int> bl(k1,k2);
			if (bl.isValid() == true && bl != br)
				sub_it[2*i+1].Initialize(g_sm,k1,k2,bc);
		}

		while (sub_it[act].isNext() == false)
			act++;
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	grid_skin_iterator_bc<dim> & operator++()
	{
		++sub_it[act];

		while (act < 2*dim && sub_it[act].isNext() == false)
			act++;

		return *this;
	}

	/*! \brief Check if there is the next element
	 *
	 * Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */
	bool isNext()
	{
		if (act < 2*dim)
		{
			//! we did not reach the end of the iterator

			return true;
		}

		//! we reach the end of the iterator
		return false;
	}

	/*! \brief Get the actual key
	 *
	 * Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	inline grid_key_dx<dim> get() const
	{
		return sub_it[act].get();
	}

	/*! \brief Reset the iterator (it restart from the beginning)
	 *
	 */
	void reset()
	{
		act = 0;
		//! Initialize to 0 the index

		for (size_t i = 0 ; i < 2*dim ; i++)
		{sub_it[i].reset();}
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_GRID_SKIN_ITERATOR_HPP_ */
