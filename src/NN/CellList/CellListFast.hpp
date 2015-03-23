/*
 * CellListStandard.hpp
 *
 *  Created on: Mar 22, 2015
 *      Author: Pietro Incardona
 */

#ifndef CELLLISTSTANDARD_HPP_
#define CELLLISTSTANDARD_HPP_



/*! \brief Class for STANDARD cell list implementation
 *
 * This class implement the STANDARD cell list, fast but memory
 * expensive. The memory allocation is (M * N_cell_max)*sizeof(ele) + M*8
 *
 * M = number of cells
 * N_cell_max = maximum number of elements in a cell
 *
 * \note Because N_cell_max >= N/M then M * N_cell_max >= O(N)
 *
 * \warning Not not use for high asymmetric distribution
 *
 * \tparam dim Dimansionality of the space
 * \tparam T type of the space float, double, complex
 * \tparam base Base structure that store the information
 *
 */
template<unsigned int dim, typename T, typename base>
class CellList<dim,T,base,FAST>
{
	// Number of slot for each cell
	size_t slot;

	// number of particle in each cell list
	openfpm::vector<size_t> cl_n;

	// elements that each cell store (each cell can store a number
	// of elements == slot )
	base cl_base;

	// Domain of the cell list
	SpaceBox<dim,T> box;

	// Unit box of the Cell list
	SpaceBox<dim,T> box_unit;

	// Grid structure of the Cell list
	grid<dim,void> gr_cell;

public:

	/*! \brief Cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param origin of the Cell list
	 * \param div grid size on each dimension
	 *
	 */
	CellList(SpaceBox<dim,T> & box, size_t div[dim], Point<dim,T> org)
	:slot(16),box(box),gr_cell(div)
	{
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the
	 * \param ele element to store
	 *
	 */
	void addElement(T (& pos)[dim], base::element ele)
	{
		typedef SpaceBox<dim,T> sb;

		// calculate the Cell id

		size_t cell_id = getCell(pos);

		// Get the number of element the cell is storing

		size_t nl = getNelements(cell_id);

		if (nl + 1 >= slot)
		{
			// we do not have enough slots reallocate the basic structure with more
			// slots

			// Create a cell-list with double of the slots

			CellList cl_tmp(slot*2);

			// copy cl_base

			for (size_t i = 0 ; i < cl_n.size() ; i++)
			{
				for (size_t j = 0 ; j < cl_n.get(i) ; j++)
					cl_tmp.cl_base.get(i*slot + j) = cl_base.get(2*slot * i + j);
			}

			// swap the memory
			cl.swap(cl_tmp);
		}

		// we have enough slot to store another neighbor element

		cl_base.get(slot * ele_id + cl_n.get(ele_id)) = ele;
		cl_n.get(ele_id)++;
	}

	/*! \brief Get the cell-id
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \param pos Point position
	 *
	 * \return the cell-id
	 *
	 */
	size_t getCell(T (& pos)[dim])
	{
		size_t cell_id = 0;

		for (size_t s = 0 ; s < dim ; s++)
		{
			ele_id += box_unit.template get<sb::p2>()[s] * gr_cell.size(s);
		}

		return cell_id;
	}

	/*! \brief Return the number of element in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return number of elements in the cell
	 *
	 */
	size_t getNelements(size_t cell_id)
	{
		return cl_n.get(cell_id);
	}

	/*! \brief Get an element in the cell
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 */
	base::element getElement(size_t cell, size_t ele)
	{
		return cl_base.get(cell * slot + ele);
	}
};


#endif /* CELLLISTSTANDARD_HPP_ */
