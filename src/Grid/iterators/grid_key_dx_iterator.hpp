/*
 * grid_key_dx_iterator.hpp
 *
 *  Created on: Dec 15, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_HPP_
#define OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_HPP_



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

template<unsigned int dim>
class grid_key_dx_iterator
{
#ifdef DEBUG
	// Actual status of the iterator, when the iterator is not initialized cannot be used
	// and reinitialize must be called
	bool initialized = false;
#endif

	grid_sm<dim,void> grid_base;

	/*! \brief return the index i of the gk key
	 *
	 * \param i index to get
	 *
	 * \return index value
	 *
	 */

	size_t get_gk(size_t i) const
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
	grid_key_dx_iterator()
	{
#ifdef DEBUG
		initialized = false;
#endif
	}

	/*! \brief Constructor from a grid_key_dx_iterator<dim>
	 *
	 * \param g_it grid_key_dx_iterator<dim>
	 */
	grid_key_dx_iterator(const grid_key_dx_iterator<dim> & g_it)
	: grid_base(g_it.grid_base)
	{
		//! Initialize to 0 the index

		for (size_t i = 0 ; i < dim ; i++)
		{
			gk.set_d(i,g_it.get_gk(i));
		}

#ifdef DEBUG
		initialized = true;
#endif
	}

	/*! \brief Constructor require a grid_sm<dim,T>
	 *
	 * \param g info of the grid on which iterate
	 */
	template<typename T> grid_key_dx_iterator(const grid_sm<dim,T> & g)
	: grid_base(g)
	{
		reset();

#ifdef DEBUG
		initialized = true;
#endif
	}

	/*! \brief Constructor from another grid_key_dx_iterator
	 *
	 * \param key_it grid_key_dx_iterator
	 */
	grid_key_dx_iterator<dim> operator=(const grid_key_dx_iterator<dim> & key_it)
	{
		grid_base = key_it.grid_base;

		//! Initialize the index using key_it

		for (size_t i = 0 ; i < dim ; i++)
		{gk.set_d(i,key_it.get_gk(i));}

		return *this;
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	grid_key_dx_iterator<dim> & operator++()
	{
		//! increment the first index

		size_t id = gk.get(0);
		gk.set_d(0,id+1);

		//! check the overflow of all the index with exception of the last dimensionality

		size_t i = 0;
		for ( ; i < dim-1 ; i++)
		{
			size_t id = gk.get(i);
			if (id >= grid_base.size(i))
			{
				// ! overflow, increment the next index

				gk.set_d(i,0);
				id = gk.get(i+1);
				gk.set_d(i+1,id+1);
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
	void set(int d, size_t sz)
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

	bool isNext()
	{
		if (gk.get(dim-1) < (long int)grid_base.size(dim-1))
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

	/*! \brief Reinitialize the grid_key_dx_iterator
	 *
	 * \param key form
	 *
	 */
	void reinitialize(const grid_key_dx_iterator<dim> & key)
	{
		grid_base = key.grid_base;
		reset();
	}

	/*! \brief Reset the iterator (it restart from the beginning)
	 *
	 */
	void reset()
	{
		//! Initialize to 0 the index

		for (size_t i = 0 ; i < dim ; i++)
		{gk.set_d(i,0);}
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_HPP_ */
