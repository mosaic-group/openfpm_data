/*
 * grid_key_dx_iterator_sp.hpp
 *
 *  Created on: Dec 15, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SP_HPP_
#define OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SP_HPP_



/**
 *
 * Grid key class iterator, iterate through a starting linearized grid element
 * to a stop linearized grid element in particular if linearize is the function
 *  that linearize all the grid_key_dx, it create an iterator that pass through
 *  Linearize^(-1)(start) Linearize^(-1)(start+1) ....... Linearize^(-1)(stop)
 *
 * \param dim dimensionality of the grid
 *
 * Usage: In general you never create object directly, but you get it from a grid_cpu or grid_gpu with
 *        getIteratorLinStartStop()
 *
 */

template<unsigned int dim>
class grid_key_dx_iterator_sp : public grid_key_dx_iterator<dim>
{
	//! stop point
	grid_key_dx<dim> gk_stop;

public:

	/*! \brief Constructor require a grid grid<dim,T>
	 *
	 * It construct an iterator from one index to another, in particular
	 * if linearize is the function that linearize all the grid_key, it
	 * create an iterator that pass through Linearize^(-1)(start)
	 * Linearize^(-1)(start+1) ....... Linearize^(-1)(stop)
	 *
	 * For example for start (1,1) and stop (3,3) the point indicated with # are
	 * explored by the iterator
	 *
	 * \verbatim
	 *
                      +-----+-----+-----+-----+-----+-----+ (6,5)
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      #-----#-----#-----#-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      #-----#-----#-----#-----#-----#-----#
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----#-----#-----#-----#-----#-----#
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                    (0,0)
	 *
	 *
	 * \endverbatim
	 *
	 *
	 * \tparam T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 * \param from starting point
	 * \param to end point
	 *
	 */
	template<typename T> grid_key_dx_iterator_sp(grid_sm<dim,T> & g, mem_id from, mem_id to)
	:grid_key_dx_iterator<dim>(g)
	 {
		//! Convert to a grid_key
		this->gk = g.InvLinId(from);

		//! Convert to a grid_key
		gk_stop = g.InvLinId(to);
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
		//! for all dimensions
		for (int i = dim-1 ; i >= 0 ; i-- )
		{
			//! check if we still have points
			if (this->gk.get(i) < gk_stop.get(i))
				return true;
			else if (this->gk.get(i) > gk_stop.get(i))
				return false;
		}

		//! (Final point) we we still have one point
		return true;
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SP_HPP_ */
