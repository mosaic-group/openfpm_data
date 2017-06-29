
/*
 * CellListIt.hpp
 *
 *  Created on: May 18, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_HPP_

template<typename T>
class ParticleIt_CellP
{
private:

	//! Cell lisy
    T & NN;

	//! Cells counter
	size_t cell_count;

	//! Particles counter
	size_t p_count;


	/*! \brief Handles incrementing of cells and particles counters
	 *
	 */
	inline void fp()
	{
		++p_count;

		const auto & SFCkeys = NN.getCellSFC().getKeys();

		if (p_count >= NN.getNelements(SFCkeys.get(cell_count)))
		{
			p_count = 0;
			++cell_count;

			while (cell_count < SFCkeys.size() && NN.getNelements(SFCkeys.get(cell_count)) == 0)
				++cell_count;
		}
	}

public:

	// Constructor
	ParticleIt_CellP(T & NN)
	:NN(NN)
	{
		reset();
	}

	// Destructor
	~ParticleIt_CellP()
	{
	}



	/*! \brief Get the next element
	 *
	 * \return cell list iterator
	 *
	 */
	inline ParticleIt_CellP & operator++()
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
		if (cell_count >= NN.getCellSFC().getKeys().size())
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
		auto cell_id = NN.getCellSFC().getKeys().get(cell_count);
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
