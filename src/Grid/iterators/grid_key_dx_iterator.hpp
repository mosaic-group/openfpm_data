/*
 * grid_key_dx_iterator.hpp
 *
 *  Created on: Dec 15, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_HPP_
#define OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_HPP_

#include "Grid/grid_sm.hpp"
#include "stencil_type.hpp"

/**
 *
 * Grid key class iterator, iterate through the grid element
 *
 * \tparam dim dimensionality of the grid
 * \tparam stencil stencil type iterator
 *
 * \note if you have a grid you can get this object from getIterator()
 *
 * ### Grid iterator declaration and usage ###
 * \snippet grid_sm_unit_tests.hpp Grid iterator test usage
 *
 */
template<unsigned int dim, typename stencil=no_stencil>
class grid_key_dx_iterator
{
#ifdef SE_CLASS1
	// Actual status of the iterator, when the iterator is not initialized cannot be used
	// and reinitialize must be called
	bool initialized = false;
#endif

	//! information of the grid where this iterator iterate
	grid_sm<dim,void> grid_base;

	//! Additional operation and information in case we do stencil
	//! operations
	stencil stl_code;

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

	//! Actual key
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
#ifdef SE_CLASS1
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

		stl_code = g_it.stl_code;

#ifdef SE_CLASS1
		initialized = true;
#endif
	}

	/*! \brief Constructor require a grid_sm<dim,T>
	 *
	 * \param g info of the grid on which iterate
	 * \param stencil_pnt stencil points
	 *
	 */
	template<typename T>
	grid_key_dx_iterator(const grid_sm<dim,T> & g,
			             const grid_key_dx<dim> (& stencil_pnt)[stencil::nsp])
	:grid_base(g)
	{
		reset();

		grid_key_dx<dim> zero;
		for (size_t i = 0 ; i < dim ; i++)	{zero.set_d(i,0);}

		// calculate the offsets for the stencil code
		stl_code.calc_offsets(g,zero,stencil_pnt);

#ifdef SE_CLASS1
		initialized = true;
#endif
	}

	/*! \brief Constructor
	 *
	 * Using this constructor you must call reinitialize
	 *
	 * \param stencil_pnt stencil points
	 *
	 */
	grid_key_dx_iterator(const grid_key_dx<dim> (& stencil_pnt)[stencil::nsp])
	{
		// calculate the offsets for the stencil code
		stl_code.set_stencil(stencil_pnt);
	}

	/*! \brief Constructor require a grid_sm<dim,T>
	 *
	 * \param g info of the grid on which iterate
	 */
	template<typename T> grid_key_dx_iterator(const grid_sm<dim,T> & g)
	: grid_base(g)
	{
		reset();

#ifdef SE_CLASS1
		initialized = true;
#endif
	}

	/*! \brief Constructor from another grid_key_dx_iterator
	 *
	 * \param key_it grid_key_dx_iterator
	 *
	 * \return itself
	 *
	 */
	inline grid_key_dx_iterator<dim> operator=(const grid_key_dx_iterator<dim> & key_it)
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

	inline grid_key_dx_iterator<dim,stencil> & operator++()
	{
		//! increment the first index

		size_t id = gk.get(0);
		gk.set_d(0,id+1);

		stl_code.increment();

		//! check the overflow of all the index with exception of the last dimensionality

		size_t i = 0;
		for ( ; i < dim-1 ; i++)
		{
			/* coverity[dead_error_begin] */
			size_t id = gk.get(i);
			if (id >= grid_base.size(i))
			{
				// ! overflow, increment the next index

				size_t idr = gk.get(i);
				gk.set_d(i,0);
				id = gk.get(i+1);
				gk.set_d(i+1,id+1);

				stl_code.adjust_offset(i,idr,grid_base);
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

	/*! \brief Get the actual position
	 *
	 * Get the actual position
	 *
	 * \return the actual key
	 *
	 */
	template<unsigned int id> inline size_t getStencil() const
	{
		return stl_code.template getStencil<id>();
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

	/*! \brief Get the actual key
	 *
	 * Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	inline const grid_key_dx<dim> & get() const
	{
		return gk;
	}

	/*! \brief Reinitialize the grid_key_dx_iterator
	 *
	 * \param key form
	 *
	 */
	inline void reinitialize(const grid_key_dx_iterator<dim> & key)
	{
		grid_base = key.getGridInfo();
		reset();
	}

	/*! \brief Get the information about the grid
	 *
	 * \return the grid info
	 *
	 */
	inline const grid_sm<dim,void> & getGridInfo() const
	{
		return grid_base;
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

	/*! \brief Calculate the stencil offset
	 *
	 * \param start_p starting point
	 *
	 */
	void calc_stencil_offset(const grid_key_dx<dim> & start_p)
	{
		// calculate the offsets for the stencil code
		stl_code.calc_offsets(grid_base,start_p);
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_HPP_ */
