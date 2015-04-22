/*
 * CellDecomposer.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: Pietro Incardona
 */

#ifndef CELLDECOMPOSER_HPP_
#define CELLDECOMPOSER_HPP_

#include "Space/SpaceBox.hpp"

/*! \brief Decompose a cell into space
 *
 * It is a convenient class for cell cell decomposition and index linearization
 *
 */

template<unsigned int dim,typename T>
class CellDecomposer_sm
{
protected:

	// Total number of cell
	size_t tot_n_cell;

	// Domain of the cell list
	SpaceBox<dim,T> box;

	// Unit box of the Cell list
	SpaceBox<dim,T> box_unit;

	// Grid structure of the Cell list
	grid_sm<dim,void> gr_cell;

	/*! \brief Initialize
	 *
	 */
	void Initialize()
	{
		tot_n_cell = 1;

		// Total number of cells and calculate the unt cell size

		for (size_t i = 0 ; i < dim ; i++)
		{
			tot_n_cell *= gr_cell.size(i);
			box_unit.setHigh(i,box.getHigh(i) / gr_cell.size(i));
		}
	}

public:

	/*! \brief Get the cell-id
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \param pos Point position
	 *
	 * \return the cell-id
	 *
	 */
	size_t getCell(const T (& pos)[dim])
	{
		size_t cell_id = pos[0] / box_unit.getHigh(0);

		for (size_t s = 1 ; s < dim ; s++)
		{
			cell_id += gr_cell.size(s) * (size_t)(pos[s] / box_unit.getHigh(s));
		}

		return cell_id;
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
	size_t getCell(const Point<dim,T> & pos)
	{
		size_t cell_id = pos.get(0) / box_unit.getHigh(0);

		for (size_t s = 1 ; s < dim ; s++)
		{
			cell_id += gr_cell.size_s(s-1) * (size_t)(pos.get(s) / box_unit.getHigh(s));
		}

		return cell_id;
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
	template<typename Mem> size_t getCell(const encapc<1,Point<dim,T>,Mem> & pos)
	{
		typedef Point<dim,T> p;

		size_t cell_id = pos.template get<p::x>()[0] / box_unit.getHigh(0);

		for (size_t s = 1 ; s < dim ; s++)
		{
			cell_id += gr_cell.size_s(s-1) * (size_t)(pos.template get<p::x>()[s] / box_unit.getHigh(s));
		}

		return cell_id;
	}

	/*! \brief Set the domain to decompose
	 *
	 * \param box Domain to decompose
	 * \param div array with the number of cells on each dimensions
	 *
	 */
	void setDimensions(SpaceBox<dim,T> & box, const size_t (&div)[dim])
	{
		this->box = box;
		this->gr_cell.setDimensions(div);
	}

	/*! \brief Set the domain to decompose
	 *
	 * \param box Domain to decompose
	 * \param div array with the number of cells on each dimensions
	 *
	 */
	void setDimensions(Box<dim,T> & box, const size_t (&div)[dim])
	{
		this->box = box;
		this->gr_cell.setDimensions(div);
		Initialize();
	}

	CellDecomposer_sm(SpaceBox<dim,T> & box, size_t (&div)[dim])
	:box(box),gr_cell(div)
	{
		Initialize();
	}


	//! default constructor
	CellDecomposer_sm()
	{

	}
};


#endif /* CELLDECOMPOSER_HPP_ */
