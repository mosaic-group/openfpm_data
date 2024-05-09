
/*
 * CellListIt.hpp
 *
 *  Created on: May 18, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_HPP_

template<typename CellList_type>
class ParticleIt_CellP
{
private:

    CellList_type & cellList;

	size_t cellIndexAct;
	size_t partIndexAct;


	/*! \brief Handles incrementing of cells and particles counters
	 *
	 */
	inline void selectValid()
	{
		++partIndexAct;

		const auto & SFCkeys = cellList.getCellSFCKeys();

		if (partIndexAct >= (long int)cellList.getNelements(SFCkeys.get(cellIndexAct)))
		{
			partIndexAct = 0;
			++cellIndexAct;

			while (cellIndexAct < SFCkeys.size() && cellList.getNelements(SFCkeys.get(cellIndexAct)) == 0)
				++cellIndexAct;
		}
	}

public:

	ParticleIt_CellP(CellList_type & cellList)
	:cellList(cellList) { reset(); }


	/*! \brief Get the next element
	 *
	 * \return cell list iterator
	 *
	 */
	inline ParticleIt_CellP & operator++()
	{
		selectValid();
		while (isNext() && this->get() >= cellList.getGhostMarker())
			selectValid();

		return *this;
	}


	/*! \brief Checks if there is a next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */
	inline bool isNext()
	{
		if (cellIndexAct >= cellList.getCellSFCKeys().size())
			return false;

		return true;
	}


	/*! \brief Get the real particle id

	 *
	 * \return the real particle id
	 *
	 */
	inline size_t get()
	{
		auto cellIndex = cellList.getCellSFCKeys().get(cellIndexAct);
		auto p = cellList.get(cellIndex,partIndexAct);

		return p;
	}

	/*! \brief Reset an iterator (set the counters to the first valid ones)
	 *
	 */
	void reset()
	{
		cellIndexAct = 0;
		partIndexAct = -1;

		this->operator++();
	}
};


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_HPP_ */
