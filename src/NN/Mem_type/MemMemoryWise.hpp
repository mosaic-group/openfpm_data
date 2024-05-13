/*
 * MemMemoryWise.hpp
 *
 *  Created on: Mar 22, 2015
 *  Last modified: June 25, 2015
 *      Authors: Pietro Incardona, Yaroslav Zaluzhnyi
 */

#ifndef CELLISTMEM_HPP_
#define CELLISTMEM_HPP_

#include "NN/CellList/CellList.hpp"

/*! \brief Class for MEMORY-WISE cell list implementation
 *
 * This class implement the MEMORY-WISE cell list
 * The memory allocation is small.
 * The memory allocation is (in byte) Size = O(N*size_of(ele))
 *
 * Where
 *
 * N = total number of elements
 * M = number of cells
 * sizeof(ele) = the size of the element the cell list is storing, example if
 *               the cell list store the particle id (64bit) is 8 byte
 *
 * \note It is useful when M >> N
 *
 * \tparam dim Dimensionality of the space
 * \tparam T type of the space float, double, complex
 *
 */
template<typename local_index = size_t>
class Mem_mw
{
	//! Base type storing information
	typedef openfpm::vector<local_index> base;

	//! each cell has a dynamic structure
	//! that store the elements in the cell
	std::unordered_map<local_index,base> cl_base;

	//! In case of invalid element return this
	typename std::remove_reference<decltype(std::declval<openfpm::vector<local_index>>().get(0))>::type invalid;

	//! ghost marker for every cell (non-ghost particles < gm (ghost marker))
	std::unordered_map<local_index,size_t> ghostMarkers;

public:

	typedef void toKernel_type;

	//! expose the type of the local index
	typedef local_index local_index_type;

	/*! \brief Initialize the data structure to zeros
	 *
	 * In this case it does nothing
	 *
	 * \param slots
	 * \param tot_n_cell total number of cells
	 *
	 */
	inline void init_to_zero(local_index slot, local_index tot_n_cell)
	{
		clear();
	}

	/*! \brief Copy two data-structure
	 *
	 * \param cell data-structure to copy
	 *
	 * \return itself
	 *
	 */
	inline Mem_mw & operator=(const Mem_mw & cell)
	{
		cl_base = cell.cl_base;
		ghostMarkers = cell.ghostMarkers;
		return *this;
	}

	/*! \brief Add an element to the cell
	 *
	 * \param cell_id cell-id
	 * \param ele element to add
	 *
	 */
	inline void addCell(local_index cell_id, typename base::value_type ele)
	{
		//add another neighbor element

		cl_base[cell_id].add(ele);
	}

	/*! \brief Remove an element from the cell
	 *
	 * \param cell cell-id
	 * \param ele element to remove
	 *
	 */
	inline void remove(local_index cell, local_index ele)
	{
		cl_base[cell].remove(ele);
	}

	/*! \brief Get the number of elements in the cell
	 *
	 * \param cell_id
	 *
	 * \return the number of elements
	 *
	 */
	inline size_t getNelements(const local_index cell_id) const
	{
		auto it = cl_base.find(cell_id);
		if (it == cl_base.end())
			return 0;

		return it->second.size();
	}

	inline auto get(local_index cell, local_index ele) -> decltype(cl_base[0].get(0)) &
	{
		auto it = cl_base.find(cell);
		if (it == cl_base.end())
			return invalid;

		return it->second.get(ele);
	}

	inline auto get(local_index cell, local_index ele) const -> decltype(cl_base.find(cell)->second.get(0)) &
	{
		auto it = cl_base.find(cell);
		if (it == cl_base.end())
			return invalid;

		return it->second.get(ele);
	}

	/*! \brief Add ghost marker to the cell
	 *
	 * \param cell_id id of the cell
	 * \param g_m ghost marker to add
	 *
	 */
	inline void addCellGhostMarkers()
	{
		for (auto it = cl_base.begin(); it != cl_base.end(); it++) {
		    ghostMarkers[it->first] = it->second.size();
		}
	}

	/*! \brief Get ghost marker of the cell
	 *
	 */
	inline size_t getGhostMarker(local_index cell_id) const
	{
		auto it = ghostMarkers.find(cell_id);
		if (it == ghostMarkers.end())
			return invalid;

		return it->second;
	}

	inline void swap(Mem_mw & cl)
	{
		cl_base.swap(cl.cl_base);
		ghostMarkers.swap(cl.ghostMarkers);
	}

	inline void swap(Mem_mw && cell)
	{
		cl_base.swap(cell.cl_base);
		ghostMarkers.swap(cell.ghostMarkers);
	}

	inline void clear()
	{
		cl_base.clear();
		ghostMarkers.clear();
	}

	inline const local_index & getStartId(size_t cell_id) const
	{
		auto it = cl_base.find(cell_id);
		if (it == cl_base.end())
			return invalid;

		return it->second.get(0);
	}

	inline const local_index & getGhostId(size_t cell_id) const
	{
		auto gm_it = ghostMarkers.find(cell_id);
		if (gm_it == ghostMarkers.end())
			return invalid;

		auto cell_it = cl_base.find(cell_id);
		if (cell_it == cl_base.end())
			return invalid;

		return cell_it->second.get(gm_it->second);
	}

	inline const local_index & getStopId(size_t cell_id) const
	{
		auto it = cl_base.find(cell_id);
		if (it == cl_base.end())
			return invalid;

		return *(&it->second.last() + 1);
	}

	inline const local_index & get_lin(const local_index * cell_id) const
	{
		return *cell_id;
	}

public:

	/*! \brief constructor
	 *
	 * \param slot number of slots (unused)
	 *
	 */
	inline Mem_mw(size_t slot)
	:invalid(0)
	{
	}

	/*! \brief Set the number of slots
	 *
	 * \param slot unused
	 *
	 */
	inline void set_slot(size_t slot)
	{}

};


#endif /* CELLISTMEM_HPP_ */
