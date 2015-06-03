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
 * It is a convenient class for cell decomposition and index linearization
 * with getCell
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

	// Cell padding
	size_t padding;

	/*! \brief Initialize all the structures
	 *
	 */
	void Initialize(const size_t pad)
	{
		tot_n_cell = 1;

		// Total number of cells and calculate the unit cell size

		for (size_t i = 0 ; i < dim ; i++)
		{
			tot_n_cell *= gr_cell.size(i);

			// Cell are padded by 1
			box_unit.setHigh(i,box.getHigh(i) / (gr_cell.size(i)-2*pad));
		}

		size_t off[dim];
		for (size_t i = 0; i < dim ; i++)
			off[i] = pad;

		padding = gr_cell.LinId(off);
	}

public:

	/*! \brief Return the underlying grid information of the cell list
	 *
	 * \return the grid infos
	 *
	 */
	grid_sm<dim,void> & getGrid()
	{
#ifdef DEBUG
		if (total_ncell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif

		return gr_cell;
	}

	/*! \brief Get the cell-ids
	 *
	 * Convert the point coordinates into the cell ids
	 *
	 * \param pos Point position
	 *
	 * \return the cell-ids ad a grid_key_dx<dim>
	 *
	 */
	grid_key_dx<dim> getCellGrid(const T (& pos)[dim])
	{
#ifdef DEBUG
		if (total_ncell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif

		grid_key_dx<dim> key;
		key.set_d(0,pos[0] / box_unit.getHigh(0));

		for (size_t s = 1 ; s < dim ; s++)
			key.set_d(s,(size_t)(pos[s] / box_unit.getHigh(s)));

		return key;
	}

	/*! \brief Get the cell-ids
	 *
	 * Convert the point coordinates into the cell ids
	 *
	 * \param pos Point position
	 *
	 * \return the cell-ids ad a grid_key_dx<dim>
	 *
	 */
	grid_key_dx<dim> getCellGrid(const Point<dim,T> pos)
	{
#ifdef DEBUG
		if (total_ncell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif

		grid_key_dx<dim> key;
		key.set_d(0,pos.get(0) / box_unit.getHigh(0));

		for (size_t s = 1 ; s < dim ; s++)
			key.set_d(s,(size_t)(pos.get(s) / box_unit.getHigh(s)));

		return key;
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
	size_t getCell(const T (& pos)[dim])
	{
#ifdef DEBUG
		if (total_ncell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif

		size_t cell_id = pos[0] / box_unit.getHigh(0);

		for (size_t s = 1 ; s < dim ; s++)
		{
			cell_id += gr_cell.size(s) * ((size_t)(pos[s] / box_unit.getHigh(s)));
		}

		return cell_id + padding;
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
#ifdef DEBUG
		if (total_ncell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif

		size_t cell_id = pos.get(0) / box_unit.getHigh(0);

		for (size_t s = 1 ; s < dim ; s++)
		{
			cell_id += gr_cell.size_s(s-1) * ((size_t)(pos.get(s) / box_unit.getHigh(s)));
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

#ifdef DEBUG
		if (total_ncell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif
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
	 * \param pad padding cell
	 *
	 */
	void setDimensions(SpaceBox<dim,T> & box, const size_t (&div)[dim], const size_t pad)
	{
		this->box = box;
		this->gr_cell.setDimensions(div);
		Initialize(pad);
	}

	/*! \brief Set the domain to decompose
	 *
	 * \param box Domain to decompose
	 * \param div array with the number of cells on each dimensions
	 * \param pad padding cell
	 *
	 */
	void setDimensions(Box<dim,T> & box, const size_t (&div)[dim], const size_t pad)
	{
		this->box = box;
		this->gr_cell.setDimensions(div);
		Initialize(pad);
	}

	/*! \brief Constructor
	 *
	 * \param box Space where is defined the cell list (it is assumed p1 = {0, .... 0})
	 * \param div Reference array to the number of divisions on each dimensions
	 * \pad cell padding
	 *
	 *  Example for div = {7,7} and pad = 1
	 *
	 * +-----------------------+
     * |p |p |p |p |p |p |p |p |
     * +-----------------------+
     * |p |  |  |  |  |  |  |p |
     * +-----------------------+
     * |p |  |  |  |  |  |  |p |
     * +-----------------------+
     * |p |  |  |  |  |  |  |p |
     * +-----------------------+
     * |p |9 |  |  |  |  |  |p |
     * +-----------------------+
     * |p |p |p |p |p |p |p |p |
     * +-----------------------+
	 *
	 * Cell with p are padding cell cell that are around but external the box, the cell number 9 that
	 * is at the origin of the box is identified with 9
	 *
	 */
	CellDecomposer_sm(SpaceBox<dim,T> & box, size_t (&div)[dim], const size_t pad)
	:box(box),gr_cell(div)
	{
		Initialize(pad);
	}


	//! default constructor
	CellDecomposer_sm()
	:tot_n_cell(0)
	{

	}
};


#endif /* CELLDECOMPOSER_HPP_ */
