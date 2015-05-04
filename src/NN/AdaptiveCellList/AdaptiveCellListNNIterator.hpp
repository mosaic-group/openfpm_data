/*
 * AdaptiveCellListNNIterator.hpp
 *
 *  Created on: May 4, 2015
 *      Author: ...
 */

#ifndef ADAPTIVECELLLISTNNITERATOR_HPP_
#define ADAPTIVECELLLISTNNITERATOR_HPP_


template<typename Cell> class AdaptiveCellNNIterator
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
	typename Cell::value_type & get()
	{
	}
};


#endif /* ADAPTIVECELLLISTNNITERATOR_HPP_ */
