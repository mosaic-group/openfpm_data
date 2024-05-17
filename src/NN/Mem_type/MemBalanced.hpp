
/*
 * MemBalanced.hpp
 *
 *  Created on: Mar 22, 2015
 *  Last modified: June 25, 2015
 *      Authors: Pietro Incardona, Yaroslav Zaluzhnyi
 */

#ifndef CELLLISTBAL_HPP_
#define CELLLISTBAL_HPP_

#include "Space/Shape/Box.hpp"
#include "util/mathutil.hpp"
#include "Space/Shape/HyperCube.hpp"

/*! \brief Class for BALANCED cell list implementation
 *
 * \tparam local_index type of local index
 *
 * This class implement the BALANCED cell list is fast (not best)
 * The memory allocation is small (not best).
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
template<typename local_index = size_t>
class Mem_bal
{
	//! vector that store the information
	typedef openfpm::vector<local_index> base;

	//! each cell has a pointer to a dynamic structure
	// that store the elements in the cell
	openfpm::vector<base> cl_base;

	//! ghost marker for every cell (non-ghost particles < gm (ghost marker))
	openfpm::vector<size_t> ghostMarkers;

	//! Invalid element
	local_index invalid;

public:

	typedef void toKernel_type;

	//! expose the type of the local index
	typedef local_index local_index_type;

	/*! \brief Initialize all to zero
	 *
	 * \param slot number of slot (unused)
	 * \param tot_n_cell total number of cells
	 *
	 */
	inline void init_to_zero(size_t slot, size_t tot_n_cell)
	{
		//resize the vector to needed number of cells

		cl_base.resize(tot_n_cell);
		ghostMarkers.resize(tot_n_cell);
		clear();
	}

	/*! \brief Copy mem balanced
	 *
	 * \param cell memory to copy
	 *
	 * \return itself
	 *
	 */
	inline Mem_bal & operator=(const Mem_bal & cell)
	{
		cl_base = cell.cl_base;
		ghostMarkers = cell.ghostMarkers;

		return *this;
	}

	/*! \brief Add an element to the cell
	 *
	 * \param cell_id id of the cell
	 * \param ele element to add
	 *
	 */
	inline void addCell(size_t cell_id, typename base::value_type ele)
	{
		//add another neighbor element

		cl_base.get(cell_id).add(ele);
	}

	/*! \brief Remove an element from the cell
	 *
	 * \param cell id of the cell
	 * \param ele element to remove
	 *
	 */
	inline void remove(local_index cell, local_index ele)
	{
		cl_base.get(cell).remove(ele);
	}

	/*! \brief Get the number of elements in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return the number of elements in the cell
	 *
	 */
	inline local_index getNelements(const local_index cell_id) const
	{
		return cl_base.get(cell_id).size();
	}

	/*! \brief Return an element from the cell
	 *
	 * \param cell id of the cell
	 * \param ele element id
	 *
	 * \return reference to the element
	 *
	 */
	inline auto get(local_index cell, local_index ele) -> decltype(cl_base.get(0).get(0)) &
	{
		return cl_base.get(cell).get(ele);
	}

	/*! \brief Return an element from the cell
	 *
	 * \param cell id of the cell
	 * \param ele element id
	 *
	 * \return reference to the element
	 *
	 */
	inline auto get(local_index cell, local_index ele) const -> decltype(cl_base.get(0).get(0)) &
	{
		return cl_base.get(cell).get(ele);
	}

	/*! \brief Add ghost marker to the cell
	 *
	 * \param cell_id id of the cell
	 * \param g_m ghost marker to add
	 *
	 */
	inline void addCellGhostMarkers()
	{
		ghostMarkers.resize(cl_base.size());

		for (int i = 0; i < cl_base.size(); ++i)
		{
			ghostMarkers.get(i) = cl_base.get(i).size();
		}
	}

	/*! \brief Get ghost marker of the cell
	 *
	 */
	inline size_t getGhostMarker(local_index cell_id) const
	{
		return ghostMarkers.get(cell_id);
	}

	/*! \brief Swap two Mem_bal
	 *
	 * \param cl element to swap with
	 *
	 */
	inline void swap(Mem_bal & cl)
	{
		cl_base.swap(cl.cl_base);
		ghostMarkers.swap(cl.ghostMarkers);
	}

	/*! \brief Swap two Mem_bal
	 *
	 * \param cell element to swap with
	 *
	 */
	inline void swap(Mem_bal && cell)
	{
		cl_base.swap(cell.cl_base);
		ghostMarkers.swap(cell.ghostMarkers);
	}

	/*! \brief Reset the object
	 *
	 *
	 */
	inline void clear()
	{
		for (size_t i = 0 ; i < cl_base.size() ; i++) {
			cl_base.get(i).clear();
			ghostMarkers.get(i) = 0;
		}
	}

	/*! \brief Get the start index of the selected element
	 *
	 * \param cell_id element
	 *
	 */
	inline const local_index & getStartId(local_index cell_id) const
	{
		if (cl_base.get(cell_id).size() == 0)
			return invalid;

		return cl_base.get(cell_id).get(0);
	}

	/*! \brief Get the index of the first ghost element
	 *
	 * \param cell_id element
	 *
	 */
	inline const local_index & getGhostId(local_index cell_id) const
	{
		if (cl_base.get(cell_id).size() == 0)
			return invalid;

		return cl_base.get(cell_id).get(this->getGhostMarker(cell_id));
	}

	/*! \brief Get the stop index of the selected element
	 *
	 * \param cell_id element
	 *
	 */
	inline const local_index & getStopId(local_index cell_id) const
	{
		if (cl_base.get(cell_id).size() == 0)
			return invalid;

		return *(&cl_base.get(cell_id).last() + 1);
	}

	/*! \brief get_lin
	 *
	 * It just return the element pointed by cell_id
	 *
	 * \param cell_id element
	 *
	 * \return the element pointed
	 *
	 */
	inline const local_index & get_lin(const local_index * cell_id) const
	{
		return *cell_id;
	}

public:

	inline Mem_bal(size_t slot)
	:invalid(0)
	{}

	inline void set_slot(size_t slot)
	{}

};


#endif /* CELLLISTBAL_HPP_ */
