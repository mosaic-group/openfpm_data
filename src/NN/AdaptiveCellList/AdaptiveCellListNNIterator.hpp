/*
 * AdaptiveCellListNNIterator.hpp
 *
 *  Created on: May 4, 2015
 *      Author: ...
 */

#ifndef ADAPTIVECELLLISTNNITERATOR_HPP_
#define ADAPTIVECELLLISTNNITERATOR_HPP_


template<unsigned int dim, typename CellS, unsigned int NNc_size, unsigned int impl> class AdaptiveCellNNIterator
{

public:

	/*! \brief
	 *
	 * Adaptive Cell NN iterator
	 *
	 * \param Cell id
	 * \param NNc Cell NN id
	 *
	 */
	AdaptiveCellNNIterator()
	{
	}
	
	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param Cell id
	 * \param NNc Cell NN id
	 *
	 */
	AdaptiveCellNNIterator(size_t cell, long int (&NNc)[NNc_size], CellS & cl)
	{
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	bool isNext()
	{

	}

	/*! \brief take the next element
	 *
	 */
	AdaptiveCellNNIterator & operator++()
	{

		return *this;
	}

	/*! \brief Get actual element
	 *
	 * \return  the actual element
	 *
	 */
	typename CellS::value_type & get()
	{
	}
};


#endif /* ADAPTIVECELLLISTNNITERATOR_HPP_ */
