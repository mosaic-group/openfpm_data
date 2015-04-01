/*
 * CelListMem.hpp
 *
 *  Created on: Mar 22, 2015
 *      Author: Pietro Incardona
 */

#ifndef CELLISTMEM_HPP_
#define CELLISTMEM_HPP_

#include <unordered_map>

/*! \brief Class for BALANCED cell list implementation
 *
 * This class implement the BALANCED cell list is fast (not best)
 * the memory allocation is small (not best).
 * The memory allocation is (in byte) Size = M*16 + N*sizeof(ele)
 *
 * Where
 *
 * N = total number of elements
 * M = number of cells
 * sizeof(ele) = the size of the element the cell list is storing, example if
 *               the cell list store the particle id (64bit) is 8 byte
 *
 * \warning Do not use for extremely fine cell list (M big)
 *
 * \tparam dim Dimensionality of the space
 * \tparam T type of the space float, double, complex
 *
 */
template<unsigned int dim, typename T, typename base>
class CellList<dim,T,MEMORY,base>
{
	// each cell has a pointer to a dynamic structure
	// that store the elements in the cell
	std::unordered_map<size_t,base *> cl_base;

	// Domain of the cell list
	SpaceBox<dim,T> box;

	// Unit box of the Cell list
	SpaceBox<dim,T> box_unit;

	// Grid structure of the Cell list
	grid_sm<dim,void> gr_cell;

public:

	/*! \brief Cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param origin of the Cell list
	 * \param div grid size on each dimension
	 *
	 */
	CellList(SpaceBox<dim,T> & box, size_t div[dim], Point<dim,T> org)
	:box(box),gr_cell(div)
	{
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the
	 * \param ele element to store
	 *
	 */
	void addElement(T (& pos)[dim], typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = getCell(pos);

		// Get the number of element the cell is storing

		cl_base[cell_id].add(ele);
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
		typedef SpaceBox<dim,T> sb;

		size_t cell_id = 0;

		for (size_t s = 0 ; s < dim ; s++)
		{
			cell_id += box_unit.template get<sb::p2>()[s] * gr_cell.size(s);
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
		return cl_base.get(cell_id)->size();
	}

	/*! \brief Get an element in the cell
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 */
	typename base::value_type getElement(size_t cell, size_t ele)
	{
		return cl_base.get(cell)->get(ele);
	}
};


#endif /* CELLISTMEM_HPP_ */
