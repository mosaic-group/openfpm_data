/*
 * CellNNIterator.hpp
 *
 *  Created on: Mar 26, 2015
 *      Author: i-bird
 */

#ifndef CELLNNITERATOR_FULL_HPP_
#define CELLNNITERATOR_FULL_HPP_

#include "mathutil.hpp"

#define FULL openfpm::math::pow(3,dim)
#define SYM  openfpm::math::pow(3,dim)/2
#define CRS openfpm::math::pow(2,dim)

#define NO_CHECK 1
#define SAFE 2

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
	long int (& NNc)[NNc_size];

	// Center cell, or cell for witch we are searching the NN-cell
	long int cell;

public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param Cell id
	 * \param NNc Cell NN id
	 *
	 */
	CellNNIterator(size_t cell, long int (&NNc)[NNc_size], Cell & cl)
	:cl(cl),NNc_id(0),cell_id(NNc[NNc_id] + cell),ele_id(0),NNc(NNc),cell(cell)
	{
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	bool isNext()
	{
		if (NNc_id >= NNc_size)
			return false;
		return true;
	}

	/*! \brief take the next element
	 *
	 */
	CellNNIterator & operator++()
	{
		ele_id++;

		if (ele_id >= cl.getNelements(cell_id))
		{
			NNc_id++;

			// No more Cell
			if (NNc_id >= NNc_size) return * this;

			cell_id = NNc[NNc_id] + cell;

			ele_id = 0;
		}

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	typename Cell::value_type & get()
	{
		return cl.get(cell_id,ele_id);
	}
};


#endif /* CELLNNITERATOR_FULL_HPP_ */
