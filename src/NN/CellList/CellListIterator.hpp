
/*
 * CellListIt.hpp
 *
 *  Created on: May 18, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_HPP_

template<typename T>
class Cell_list_iterator
{
private:
    T & NN;

	// Cells counter
	size_t cell_count;

	// Particles counter
	size_t p_count;


	/*! \brief Handles incrementing of cells and particles counters
	 *
	 */
	inline void fp()
	{
		++p_count;

		if (p_count >= NN.getNelements(NN.getKeys().get(cell_count)))
		{
			p_count = 0;
			++cell_count;

			while (cell_count < NN.getKeys().size() && NN.getNelements(NN.getKeys().get(cell_count)) == 0)
			{
				++cell_count;
			}
		}
	}

public:

	// Constructor
	Cell_list_iterator(T & NN)
	:NN(NN)
	{
		reset();
	}

	// Destructor
	~Cell_list_iterator()
	{
	}



	/*! \brief Get the next element
	 *
	 * \return cell list iterator
	 *
	 */
	inline Cell_list_iterator operator++()
	{
		fp();
		while (isNext() && this->get() >= NN.get_gm())
		{
			fp();
		}

		return *this;
	}


	/*! \brief Checks if there is a next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */
	inline bool isNext()
	{
		if (cell_count >= NN.getKeys().size())
		{
			return false;
		}

		return true;
	}


	/*! \brief Get the real particle id

	 *
	 * \return the real particle id
	 *
	 */
	inline size_t get()
	{
		auto cell_id = NN.getKeys().get(cell_count);
		auto p = NN.get(cell_id,p_count);

		return p;
	}

	/*! \brief Reset an iterator (set the counters to the first valid ones)
	 *
	 */
	void reset()
	{
		cell_count = 0;
		p_count = -1;

		this->operator++();
	}
};


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_HPP_ */
