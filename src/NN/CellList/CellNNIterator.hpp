/*
 * CellNNIterator.hpp
 *
 *  Created on: Mar 26, 2015
 *      Author: i-bird
 */

#ifndef CELLNNITERATOR_HPP_
#define CELLNNITERATOR_HPP_

#include "mathutil.hpp"

template<unsigned int dim, typename Cell> class CellNNIterator
{
	// Cell list
	Cell & cl;

	// Actual NNc_id;
	size_t NNc_id;

	// actual cell id = NNc[NNc_id] stored for performance reason
	size_t cell_id;

	// actual element id
	size_t ele_id;

	// NN cell id
	size_t & NNc[openfpm::math::pow(3,dim)];


	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param Cell id
	 * \param NNc Cell NN id
	 *
	 */
	CellNNIterator(size_t cell, size_t (&NNc)[openfpm::math::pow(3,dim)])
	:NNc_id(0),cell_id(NNc[NNc_id]),ele_id(0),NNc(NNc)
	{
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	bool isNext()
	{
		if (cell_id >= openfpm::math::pow(3,dim))
			return false;
		return true;
	}

	/*! \brief take the next element
	 *
	 */
	CellNNIterator & operator++()
	{
		ele_id++;

		if (ele_id > cl.getNelements(cell_id))
		{
			NNc_id++;
			ele_id = 0;
		}

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the nect element object
	 *
	 */
	Cell::value_type & get()
	{
		return cl.get(cell_id,ele_id);
	}
};


#endif /* CELLNNITERATOR_HPP_ */
