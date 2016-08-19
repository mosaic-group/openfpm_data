/*
 * CellNNIterator.hpp
 *
 *  Created on: Mar 26, 2015
 *      Author: i-bird
 */

#ifndef CELLNNITERATOR_FULL_HPP_
#define CELLNNITERATOR_FULL_HPP_

#include "util/mathutil.hpp"

#define FULL openfpm::math::pow(3,dim)
#define SYM  openfpm::math::pow(3,dim)/2 + 1
#define CRS openfpm::math::pow(2,dim)

#define NO_CHECK 1
#define SAFE 2

/*! \brief Iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 * \tparam impl implementation specific options NO_CHECK do not do check on access, SAFE do check on access
 *
 */
template<unsigned int dim, typename Cell,unsigned int NNc_size, unsigned int impl> class CellNNIterator
{
	// Cell list
	Cell & cl;

	// Actual NNc_id;
	size_t NNc_id;

	// actual cell id = NNc[NNc_id]+cell stored for performance reason
	size_t cell_id;

	// actual element id
	size_t ele_id;

	// NN cell id
	const long int (& NNc)[NNc_size];

	// Center cell, or cell for witch we are searching the NN-cell
	const long int cell;

	/*! \brief Select non-empty cell
	 *
	 */
	inline void selectValid()
	{
		while (ele_id >= cl.getNelements(cell_id))
		{
			NNc_id++;

			// No more Cell
			if (NNc_id >= NNc_size) return;

			cell_id = NNc[NNc_id] + cell;

			ele_id = 0;
		}
	}

public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 *
	 */
	inline CellNNIterator(size_t cell, const long int (&NNc)[NNc_size], Cell & cl)
	:cl(cl),NNc_id(0),cell_id(NNc[NNc_id] + cell),ele_id(0),NNc(NNc),cell(cell)
	{
		selectValid();
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	inline bool isNext()
	{
		if (NNc_id >= NNc_size)
			return false;
		return true;
	}

	/*! \brief take the next element
	 *
	 */
	inline CellNNIterator & operator++()
	{
		ele_id++;

		selectValid();

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline typename Cell::value_type & get()
	{
		return cl.get(cell_id,ele_id);
	}
};

/*! \brief it iterate through the elements of a cell
 *
 * In general you do not create this object you get it from the CellList structures
 *
 * \tparam Cell cell type
 *
 */
template<typename Cell> class CellIterator
{
	// Cell list
	Cell & cl;

	// actual element id inside the cell
	size_t ele_id;

	// selected cell
	const long int cell;

public:

	/*! \brief Cell iterator
	 *
	 * \param cell Cell id
	 * \param cl Cell on which iterate
	 *
	 */
	inline CellIterator(const size_t cell, Cell & cl)
	:cl(cl),ele_id(0),cell(cell)
	{
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	inline bool isNext()
	{
		return cl.getNelements(cell) > ele_id;
	}

	/*! \brief take the next element
	 *
	 */
	inline CellIterator & operator++()
	{
		ele_id++;

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline typename Cell::value_type & get()
	{
		return cl.get(cell,ele_id);
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline const typename Cell::value_type & get() const
	{
		return cl.get(cell,ele_id);
	}
};

#endif /* CELLNNITERATOR_FULL_HPP_ */
