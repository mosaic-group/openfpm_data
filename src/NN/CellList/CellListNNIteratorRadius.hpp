/*
 * CellListNNIteratorRadius.hpp
 *
 *  Created on: Aug 17, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTNNITERATORRADIUS_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTNNITERATORRADIUS_HPP_


/*! \brief Iterator for the neighborhood of the cell structures with free radius
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells (inside a radius)
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 * \tparam impl implementation specific options NO_CHECK do not do check on access, SAFE do check on access
 *
 */
template<unsigned int dim, typename Cell, unsigned int impl> class CellNNIteratorRadius
{
	// Cell list
	Cell & cl;

	// Actual NNc_id;
	size_t NNc_id;

	// actual cell id = NNc[NNc_id]+cell stored for performance reason
	size_t cell_id;

	// actual element id
	size_t ele_id;

	// NN number of neighborhood cells
	const openfpm::vector<long int> & NNc;

	// Center cell, or cell for witch we are searching the NN-cell
	const long int cell;

	/*! \brief Select non-empty cell
	 *
	 */
	inline void selectValid()
	{
		while (ele_id >= cl.getNelements(cell_id))
		{
			NNc_id++;

			// No more Cell
			if (NNc_id >= NNc.size()) return;

			cell_id = NNc.get(NNc_id) + cell;

			ele_id = 0;
		}
	}

public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 *
	 */
	inline CellNNIteratorRadius(size_t cell, const openfpm::vector<long int> &NNc, Cell & cl)
	:cl(cl),NNc_id(0),cell_id(NNc.get(NNc_id) + cell),ele_id(0),NNc(NNc),cell(cell)
	{
#ifdef SE_CLASS1
		if (cell_id < 0)
			std::cerr << "Error " << __FILE__ ":" << __LINE__ << " cell_id is negative, please check the the padding is chosen correctly." <<
			                                                      "Remember, if you choose a radius that span N neighborhood cell-list, padding must be one" << std::endl;
#endif

		selectValid();
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	inline bool isNext()
	{
		if (NNc_id >= NNc.size())
			return false;
		return true;
	}

	/*! \brief take the next element
	 *
	 */
	inline CellNNIteratorRadius & operator++()
	{
		ele_id++;

		selectValid();

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline typename Cell::value_type & get()
	{
		return cl.get(cell_id,ele_id);
	}
};


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTNNITERATORRADIUS_HPP_ */
