/*
 * ParticleIt_Cells.hpp
 *
 *  Created on: Mar 5, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_PARTICLEIT_CELLS_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_PARTICLEIT_CELLS_HPP_

/*! \brief This iterator iterate across the particles of a Cell-list following the Cell structure
 *
 *
 * \tparam dim Dimensionality
 * \tparam CellListType type of the cell-list
 * \tparam CellIterator type of iterator over the cell
 *
 */
template<unsigned int dim,typename CellListType, typename CellIt> class ParticleIt_Cells
{
private:

	//! starting position
	const size_t * start;

	//! stop position
	const size_t * stop;

	//! Celllist type
	CellListType & cli;

	//! grid iterator for cells iteration
	CellIt it;


	/*! \brief Adjust the counters to reach a valid particle element
	 *
	 *
	 */
	void selectValid()
	{
		while (it.isNext())
		{
			auto cell = it.get();
			size_t s_cell = cli.getGrid().LinId(cell);

			// Get the starting particle
			start = &cli.getStartId(s_cell);

			// Get the stop particle
			stop = &cli.getStopId(s_cell);

			if (start == stop)
				++it;
			else
				break;
		}
	}

public:

	/*! \brief Initialize the iterator
	 *
	 * \param cli Cell-list
	 * \param size_t internal grid cells
	 * \param
	 *
	 */
	ParticleIt_Cells(CellListType & cli, CellIt & cit)
	:cli(cli),it(cit)
	{
		selectValid();
	}

	/*! \brief Increment to the next particle
	 *
	 * \return The actual particle iterator
	 *
	 */
	ParticleIt_Cells & operator++()
	{
		++start;

		if (start == stop)
		{
			++it;
			selectValid();
		}

		return *this;
	}

	/*! \brief Return true if there is the next particle
	 *
	 * \return true if there is a new point
	 *
	 */
	bool isNext()
	{
		return it.isNext();
	}

	/*! \brief Get the actual particle id
	 *
	 * \return the actual particle id
	 *
	 */
	size_t get()
	{
		return *start;
	}
};


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_PARTICLEIT_CELLS_HPP_ */
