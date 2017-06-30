/*
 * grid_key_dx_iterator_stencil.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_STENCIL_HPP_
#define OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_STENCIL_HPP_


#include "Grid/grid_sm.hpp"

/**
 *
 * Grid key class iterator, iterate through the grid element
 *
 * \param dim dimensionality of the grid
 *
 * \note if you have a grid you can get this object from getIterator()
 *
 * ### Grid iterator declaration and usage
 * \snippet grid_unit_tests.hpp Grid iterator test usage
 *
 */

template<unsigned int dim, unsigned int Np>
class grid_key_dx_iterator_stencil
{
#ifdef SE_CLASS1
	// Actual status of the iterator, when the iterator is not initialized cannot be used
	// and reinitialize must be called
	bool initialized = false;
#endif

	//! information about the grid
	grid_sm<dim,void> grid_base;

	//! set of offsets for the stencil
	long int stencil_offset[Np];

	/*! \brief return the index i of the gk key
	 *
	 * \param i index to get
	 *
	 * \return index value
	 *
	 */

	inline size_t get_gk(size_t i) const
	{
		return gk.get(i);
	}

protected:

	grid_key_dx<dim> gk;

public:

	/*! \brief Default constructor
	 *
	 * \warning entremly unsafe
	 * Before use the iterator you have call reinitialize
	 *
	 */
	inline grid_key_dx_iterator_stencil()
	{
#ifdef SE_CLASS1
		initialized = false;
#endif
	}

	/*! \brief Constructor from a grid_key_dx_iterator<dim>
	 *
	 * \param g_it grid_key_dx_iterator<dim>
	 */
	inline grid_key_dx_iterator_stencil(const grid_key_dx_iterator_stencil<dim,Np> & g_it)
	: grid_base(g_it.grid_base)
	{
		//! Initialize to 0 the index

		for (size_t i = 0 ; i < dim ; i++)
		{gk.set_d(i,g_it.get_gk(i));}

		// Copy the offset part

		for (size_t i = 0 ; i < Np ; i++)
		{stencil_offset[i] = g_it.stencil_offset[i];}

#ifdef SE_CLASS1
		initialized = true;
#endif
	}

	/*! \brief Constructor require a grid_sm<dim,T>
	 *
	 * \param g info of the grid on which iterate
	 */
	template<typename T> grid_key_dx_iterator_stencil(const grid_sm<dim,T> & g, const grid_key_dx<dim> (& stencil)[Np])
	: grid_base(g)
	{
		reset();

		// calculate the offsets

		for (size_t i = 0 ; i < Np ; i++)
		{
			grid_key_dx<dim> zero;

			for (size_t k = 0 ; k < dim ; k++)	{zero.set_d(k,0);}

			zero = zero + stencil[i];

			stencil_offset[i] = g.LinId(zero);
		}

#ifdef SE_CLASS1
		initialized = true;
#endif


	}

	/*! \brief Constructor from another grid_key_dx_iterator
	 *
	 * \param key_it grid_key_dx_iterator
	 */
	inline grid_key_dx_iterator_stencil<dim,Np> operator=(const grid_key_dx_iterator_stencil<dim,Np> & key_it)
	{
		grid_base = key_it.grid_base;

		//! Initialize the index using key_it

		for (size_t i = 0 ; i < dim ; i++)
		{gk.set_d(i,key_it.get_gk(i));}

		for (size_t i = 0 ; i < Np ; i++)
		{stencil_offset[i] = key_it.stencil_offset[i];}

		return *this;
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	inline grid_key_dx_iterator_stencil<dim,Np> & operator++()
	{
		//! increment the first index

		size_t id = gk.get(0);
		gk.set_d(0,id+1);

		// update the offsets
		for (size_t i = 0 ; i < Np ; i++)
			stencil_offset[i] += 1;

		//! check the overflow of all the index with exception of the last dimensionality

		size_t i = 0;
		for ( ; i < dim-1 ; i++)
		{
			size_t id = gk.get(i);
			if (id >= grid_base.size(i))
			{
				// ! overflow, increment the next index

				size_t idr = gk.get(i);
				gk.set_d(i,0);
				id = gk.get(i+1);
				gk.set_d(i+1,id+1);

				size_t str_dw = (i == 0)?1:grid_base.size_s(i-1);

				// update the offsets
				for (size_t k = 0 ; k < Np ; k++)
				{stencil_offset[k] += -str_dw*idr + grid_base.size_s(i);}
			}
			else
			{
				break;
			}
		}

		return *this;
	}

	/*! \brief Set the dimension
	 *
	 * Set the dimension
	 *
	 * \param d is the dimension
	 * \param sz set the counter to sz
	 *
	 */
	inline void set(int d, size_t sz)
	{
		// set the counter dim to sz

		gk.set_d(d,sz);
	}

	/*! \brief Check if there is the next element
	 *
	 * Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */

	inline bool isNext()
	{
		if (gk.get(dim-1) < (long int)grid_base.size(dim-1))
		{
			//! we did not reach the end of the grid

			return true;
		}

		//! we reach the end of the grid
		return false;
	}

	/*! \brief Get the actual position
	 *
	 * Get the actual position
	 *
	 * \return the actual key
	 *
	 */
	inline const grid_key_dx<dim> & getLoc() const
	{
		return gk;
	}


	/*! \brief Get the actual position
	 *
	 * Get the actual position
	 *
	 * \return the actual key
	 *
	 */
	template<unsigned int id> inline size_t get() const
	{
		return stencil_offset[id];
	}

	/*! \brief Reinitialize the grid_key_dx_iterator
	 *
	 * \param key form
	 *
	 */
	inline void reinitialize(const grid_key_dx_iterator_stencil<dim,Np> & key)
	{
		grid_base = key.grid_base;
		reset();
	}

	/*! \brief Reset the iterator (it restart from the beginning)
	 *
	 */
	inline void reset()
	{
		// Initialize to 0 the index

		for (size_t i = 0 ; i < dim ; i++)
		{gk.set_d(i,0);}

		// here we check if grid have a size equal to zero or negative
		// in this case the grid has no points

		for (size_t i = 0 ; i < dim ; i++)
		{
			// If the size of the grid is zero in any dimension set the iterator
			// to the end point
			if (grid_base.size(i) == 0)
				gk.set_d(dim-1,grid_base.size(dim-1));
		}
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_STENCIL_HPP_ */
