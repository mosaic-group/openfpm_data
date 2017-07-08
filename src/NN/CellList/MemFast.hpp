/*
 * MemFast.hpp
 *
 *  Created on: Mar 22, 2015
 *      Author: Pietro Incardona
 */

#ifndef MEMFAST_HPP_
#define MEMFAST_HPP_

#include "Space/SpaceBox.hpp"
#include "util/mathutil.hpp"
#include "Space/Shape/HyperCube.hpp"
#include "CellListIterator.hpp"
#include <unordered_map>
#include "util/common.hpp"
#include "Vector/map_vector.hpp"

/*! \brief It is a class that work like a vector of vector
 *
 * It is a class that work like a vector(1) of vector(2). To emulate
 * the vector of vector it use a 1D array of size N_ele * N_max_slot
 * where N_ele is the number of elements in vector(1) and N_max_slot
 * is the maximum number of elements across the vectors
 *
 */
class Mem_fast
{
	//! Number of slot for each cell
	size_t slot;

	//! number of particle in each cell list
	openfpm::vector<size_t> cl_n;

	//! base that store the data
	typedef typename openfpm::vector<size_t> base;

	//! elements that each cell store (each cell can store a number
	//! of elements == slot )
	base cl_base;

	/*! \brief realloc the data structures
	 *
	 *
	 */
	inline void realloc()
	{
		// we do not have enough slots reallocate the basic structure with more
		// slots
		base cl_base_(2*slot * cl_n.size());

		// copy cl_base
		for (size_t i = 0 ; i < cl_n.size() ; i++)
		{
			for (size_t j = 0 ; j < cl_n.get(i) ; j++)
				cl_base_.get(2*i*slot + j) = cl_base.get(slot * i + j);
		}

		// Double the number of slots
		slot *= 2;

		// swap the memory
		cl_base.swap(cl_base_);
	}

protected:

	/*! \brief Initialize the data to zero
	 *
	 * \param slot number of slot for each cell
	 * \param tot_n_cell total number of cells
	 *
	 */
	inline void init_to_zero(size_t slot, size_t tot_n_cell)
	{
		this->slot = slot;

		// create the array that store the number of particle on each cell and se it to 0

		cl_n.resize(tot_n_cell);
		cl_n.fill(0);

		// create the array that store the cell id

		cl_base.resize(tot_n_cell * slot);
	}

	/*! \brief copy an object Mem_fast
	 *
	 * \param mem Mem_fast to copy
	 *
	 */
	inline void operator=(const Mem_fast & mem)
	{
		slot = mem.slot;

		cl_n = mem.cl_n;
		cl_base = mem.cl_base;
	}

	/*! \brief copy an object Mem_fast
	 *
	 * \param mem Mem_fast to copy
	 *
	 */
	inline void operator=(Mem_fast && mem)
	{
		this->swap(mem);
	}

	/*! \brief Add an element to the cell
	 *
	 * \param cell_id id of the cell
	 * \param ele element to add
	 *
	 */
	inline void addCell(size_t cell_id, typename base::value_type ele)
	{
		// Get the number of element the cell is storing

		size_t nl = getNelements(cell_id);

		if (nl + 1 >= slot)
		{
			realloc();
		}

		// we have enough slot to store another neighbor element

		cl_base.get(slot * cell_id + cl_n.get(cell_id)) = ele;
		cl_n.get(cell_id)++;
	}

	/*! \brief Add an element to the cell
	 *
	 * \param cell_id id of the cell
	 * \param ele element to add
	 *
	 */
	inline void add(size_t cell_id, typename base::value_type ele)
	{
		// add the element to the cell

		this->addCell(cell_id,ele);
	}

	/*! \brief Get an element in the cell
	 *
	 * \param cell id of the cell
	 * \param ele element id in the cell
	 *
	 * \return the reference to the selected element
	 *
	 */
	inline auto get(size_t cell, size_t ele) -> decltype(cl_base.get(cell * slot + ele)) &
	{
		return cl_base.get(cell * slot + ele);
	}


	/*! \brief Get an element in the cell
	 *
	 * \param cell id of the cell
	 * \param ele element id in the cell
	 *
	 * \return the reference to the selected element
	 *
	 */
	inline auto get(size_t cell, size_t ele) const -> decltype(cl_base.get(cell * slot + ele)) &
	{
		return cl_base.get(cell * slot + ele);
	}


	/*! \brief Remove an element in the cell
	 *
	 * \param cell id of the cell
	 * \param ele element id to remove
	 *
	 */
	inline void remove(size_t cell, size_t ele)
	{
		cl_n.get(cell)--;
	}

	/*! \brief Get the number of elements in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return the number of elements in the cell
	 *
	 */
	inline size_t getNelements(const size_t cell_id) const
	{
		return cl_n.get(cell_id);
	}


	/*! \brief swap to Mem_fast object
	 *
	 * \param mem object to swap the memory with
	 *
	 */
	inline void swap(Mem_fast & mem)
	{
		cl_n.swap(mem.cl_n);
		cl_base.swap(mem.cl_base);

		size_t cl_slot_tmp = mem.slot;
		mem.slot = slot;
		slot = cl_slot_tmp;
	}

	/*! \brief swap to Mem_fast object
	 *
	 * \param mem object to swap the memory with
	 *
	 */
	inline void swap(Mem_fast && mem)
	{
		slot = mem.slot;

		cl_n.swap(mem.cl_n);
		cl_base.swap(mem.cl_base);
	}

	/*! \brief Delete all the elements in the Cell-list
	 *
	 *
	 *
	 */
	inline void clear()
	{
		for (size_t i = 0 ; i < cl_n.size() ; i++)
			cl_n.get(i) = 0;
	}

	/*! \brief Get the first element of a cell (as reference)
	 *
	 * \param cell_id cell-id
	 *
	 * \return a reference to the first element
	 *
	 */
	inline const size_t & getStartId(size_t cell_id) const
	{
		return cl_base.get(cell_id*slot);
	}

	/*! \brief Get the last element of a cell (as reference)
	 *
	 * \param cell_id cell-id
	 *
	 * \return a reference to the last element
	 *
	 */
	inline const size_t & getStopId(size_t cell_id) const
	{
		return cl_base.get(cell_id*slot+cl_n.get(cell_id));
	}

	/*! \brief Just return the value pointed by part_id
	 *
	 * \param part_id
	 *
	 * \return the value pointed by part_id
	 *
	 */
	inline const size_t & get_lin(const size_t * part_id) const
	{
		return *part_id;
	}

public:

	/*! \brief Constructor
	 *
	 * \param slot number of slot for each cell
	 *
	 */
	inline Mem_fast(size_t slot)
	:slot(slot)
	{}

	/*! \brief Set the number of slot for each cell
	 *
	 * \param number of slot
	 *
	 */
	inline void set_slot(size_t slot)
	{
		this->slot = slot;
	}

};


#endif /* CELLLISTSTANDARD_HPP_ */
