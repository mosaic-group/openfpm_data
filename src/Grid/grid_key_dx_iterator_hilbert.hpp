/*
 * grid_key_dx_iterator_hilbert.hpp
 *
 *  Created on: Feb 24, 2016
 *      Author: yaroslav
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_KEY_DX_ITERATOR_HILBERT_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_KEY_DX_ITERATOR_HILBERT_HPP_

extern "C"
{
#include "hilbertKey.h"
}


/*
 *
 * Grid key class iterator, iterate through the grid elements following an
 * hilbert space filling curve
 *
 * \param dim dimensionality of the grid
 *
 * ### Grid iterator declaration and usage
 * \snippet grid_unit_tests.hpp Grid iterator test usage
 *
 */

template<unsigned int dim>
class grid_key_dx_iterator_hilbert
{
#ifdef SE_CLASS1
	//! Actual status of the iterator, when the iterator is not initialized cannot be used
	//! and reinitialize must be called
	bool initialized = false;
#endif

	//! Actual position
	uint64_t hkey = 0;

	//! Order of a hilbert curve
	size_t m;

	//! Size of the hilbert grid in each dimension
	grid_sm<dim,void> grid_base;

protected:

	//! Actual position in the grid
	grid_key_dx<dim> gk;

public:

	/*! \brief Constructor
	 *
	 * m is the order of the hilber curve m=2 produce an hilber curve
	 * passing 2^(2) point in each direction so 4x4 points in 2D 4x4x4 in 3D
	 *
	 * \param m order of the hilber curve
	 *
	 */
	grid_key_dx_iterator_hilbert(int32_t m)
	:m(m)
	{
		// create from m the correct grid_sm

		size_t dims[dim];

		//Set the 2^m value to the each elements of dims[]
		for (size_t i = 0 ; i < dim ; i++)
			dims[i] = 1 << m;

		//Set grid_sm dimensions
		grid_base.setDimensions(dims);

		reset();

#ifdef SE_CLASS1
		initialized = true;
#endif
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	grid_key_dx_iterator_hilbert<dim> & operator++()
	{
		//An integer to handle errors
		int err;

		hkey++;

		//Array to handle output
		uint64_t nextCoord[dim];

		//Get the coordinates of the next cell
		getIntCoordFromHKey(nextCoord, m, dim, hkey, &err);

		//Set the new coordinates
		for (size_t i = 0; i < dim; i++)
			gk.set_d(i, nextCoord[i]);

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
		if ( hkey < (size_t)1 << (m*dim))
		{
			//! we did not reach the end of the grid

			return true;
		}

		//! we reach the end of the grid
		return false;
	}

	/*! \brief Get the actual key
	 *
	 * Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	const grid_key_dx<dim> & get()
	{
		return gk;
	}


	/*! \brief Reset the iterator (it restart from the beginning)
	 *
	 */
	void reset()
	{
		//! Initialize to 0 the index
		for (size_t i = 0 ; i < dim ; i++)
			gk.set_d(i,0);
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_GRID_KEY_DX_ITERATOR_HILBERT_HPP_ */
