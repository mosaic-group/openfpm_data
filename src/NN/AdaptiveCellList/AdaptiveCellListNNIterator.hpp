/*
 * AdaptiveCellListNNIterator.hpp
 *
 *  Created on: May 4, 2015
 *      Author: ...
 */

#ifndef ADAPTIVECELLLISTNNITERATOR_HPP_
#define ADAPTIVECELLLISTNNITERATOR_HPP_

#include <iostream>

template<unsigned int dim, typename CellS, unsigned int NNc_size, unsigned int impl> class AdaptiveCellNNIterator
{

	CellS & cl;

	size_t ele_id;

public:
	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 */
	AdaptiveCellNNIterator(CellS & cl)
		:cl(cl), ele_id(0)
	{
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	bool isNext()
	{
		return ele_id < cl.size();
	}

	/*! \brief take the next element
	 *
	 */
	AdaptiveCellNNIterator & operator++()
	{
		if(ele_id < cl.size()) ele_id++;
		return *this;
	}

	/*! \brief Get actual element
	 *
	 * \return  the actual element
	 *
	 */
	inline size_t& get()
	{
		return cl.get(ele_id);
	}
};


#endif /* ADAPTIVECELLLISTNNITERATOR_HPP_ */
