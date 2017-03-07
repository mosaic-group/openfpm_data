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
 *
 */
template<unsigned int dim,typename CellListType> class ParticleIt_Cells
{
private:

	//! starting position
	const size_t * start;

	//! stop position
	const size_t * stop;

	//! Actual cell
	size_t cid;

	//! List of all the domain cells
	const openfpm::vector<size_t> &dom_cell;

	//! Celllist type
	CellListType & cli;

	/*! \brief Adjust the counters to reach a valid particle element
	 *
	 *
	 */
	void selectValid()
	{
		while (start == stop)
		{
			cid++;

			size_t s_cell;
			if (cid >= dom_cell.size())
				return;
			else
				s_cell = dom_cell.get(cid);

			// Get the starting particle
			start = &cli.getStartId(s_cell);

			// Get the stop particle
			stop = &cli.getStopId(s_cell);
		}
	}

public:

	/*! \brief Initialize the iterator
	 *
	 * \param cli Cell-list
	 * \param dom_cell domain cell
	 * \param anom_dom_cell anomalous domain cell
	 * \param NNc_sym symmetric neighborhood
	 *
	 *
	 */
	ParticleIt_Cells(CellListType & cli,
					 const openfpm::vector<size_t> & dom_cell)
	:cid(0),dom_cell(dom_cell),cli(cli)
	{
		size_t s_cell;
		if (dom_cell.size() != 0)
			s_cell = dom_cell.get(0);
		else
		{
			cid = 1;
			start = NULL;
			stop = NULL;

			return;
		}

		// Get the starting particle
		start = &cli.getStartId(s_cell);

		// Get the stop particle
		stop = &cli.getStopId(s_cell);

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

		selectValid();

		return *this;
	}

	/*! \brief Return true if there is the next particle
	 *
	 * \return true if there is a new point
	 *
	 */
	bool isNext()
	{
		return cid < dom_cell.size();
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
