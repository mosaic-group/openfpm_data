/*
 * MemFast.hpp
 *
 *  Created on: Mar 22, 2015
 *      Author: Pietro Incardona
 */

#ifndef MEMFAST_HPP_
#define MEMFAST_HPP_

#include "config.h"
#include "Space/SpaceBox.hpp"
#include "util/mathutil.hpp"
#include "Space/Shape/HyperCube.hpp"
#include "NN/CellList/CellListIterator.hpp"
#include <unordered_map>
#include "util/common.hpp"
#include "Vector/map_vector.hpp"

template <typename Memory, typename local_index>
class Mem_fast_ker
{
	//! Number of slot for each cell
	local_index slot;

	//! number of particle in each cell list
	openfpm::vector_gpu_ker<aggregate<local_index>> cl_n;

	//! base that store the data
	typedef openfpm::vector_gpu_ker<aggregate<local_index>> base;

	//! elements that each cell store (each cell can store a number
	//! of elements == slot )
	base cl_base;

public:

	Mem_fast_ker(openfpm::vector_gpu_ker<aggregate<local_index>> cl_n, openfpm::vector_gpu_ker<aggregate<local_index>> cl_base)
	:cl_n(cl_n),cl_base(cl_base)
	{}
};

/*! \brief It is a class that work like a vector of vector
 *
 * \tparam local_index type used for the local index
 *
 * It is a class that work like a vector(1) of vector(2). To emulate
 * the vector of vector it use a 1D array of size N_ele * N_max_slot
 * where N_ele is the number of elements in vector(1) and N_max_slot
 * is the maximum number of elements across the vectors
 *
 */
template <typename Memory = HeapMemory, typename local_index = size_t>
class Mem_fast
{
	//! Number of slot for each cell
	local_index slot;

	//! number of particle in each cell list
	openfpm::vector<aggregate<local_index>,Memory> cl_n;

	//! base that store the data
	typedef typename openfpm::vector<aggregate<local_index>,Memory> base;

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
			for (local_index j = 0 ; j < cl_n.template get<0>(i) ; j++)
			{cl_base_.template get<0>(2*i*slot + j) = cl_base.template get<0>(slot * i + j);}
		}

		// Double the number of slots
		slot *= 2;

		// swap the memory
		cl_base.swap(cl_base_);
	}


public:

	typedef Mem_fast_ker<Memory,local_index> toKernel_type;

	/*! \brief return the number of elements
	 *
	 * \return the number of elements
	 *
	 */
	inline size_t size()
	{
		return cl_n.size();
	}

	/*! \brief Initialize the data to zero
	 *
	 * \param slot number of slot for each cell
	 * \param tot_n_cell total number of cells
	 *
	 */
	inline void init_to_zero(local_index slot, local_index tot_n_cell)
	{
		this->slot = slot;

		// create the array that store the number of particle on each cell and se it to 0

		cl_n.resize(tot_n_cell);
		cl_n.template fill<0>(0);

		// create the array that store the cell id

		cl_base.resize(tot_n_cell * slot);
	}

	/*! \brief copy an object Mem_fast
	 *
	 * \param mem Mem_fast to copy
	 *
	 */
	inline void operator=(const Mem_fast<Memory,local_index> & mem)
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
	inline void operator=(Mem_fast<Memory,local_index> && mem)
	{
		this->swap(mem);
	}

	/*! \brief Add an element to the cell
	 *
	 * \param cell_id id of the cell
	 * \param ele element to add
	 *
	 */
	inline void addCell(local_index cell_id, local_index ele)
	{
		// Get the number of element the cell is storing

		local_index nl = getNelements(cell_id);

		if (nl + 1 >= slot)
		{
			realloc();
		}

		// we have enough slot to store another neighbor element

		cl_base.template get<0>(slot * cell_id + cl_n.template get<0>(cell_id)) = ele;
		cl_n.template get<0>(cell_id)++;
	}

	/*! \brief Add an element to the cell
	 *
	 * \param cell_id id of the cell
	 * \param ele element to add
	 *
	 */
	inline void add(local_index cell_id, local_index ele)
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
	inline auto get(local_index cell, local_index ele) -> decltype(cl_base.template get<0>(cell * slot + ele)) &
	{
		return cl_base.template get<0>(cell * slot + ele);
	}


	/*! \brief Get an element in the cell
	 *
	 * \param cell id of the cell
	 * \param ele element id in the cell
	 *
	 * \return the reference to the selected element
	 *
	 */
	inline auto get(local_index cell, local_index ele) const -> decltype(cl_base.template get<0>(cell * slot + ele)) &
	{
		return cl_base.template get<0>(cell * slot + ele);
	}


	/*! \brief Remove an element in the cell
	 *
	 * \param cell id of the cell
	 * \param ele element id to remove
	 *
	 */
	inline void remove(local_index cell, local_index ele)
	{
		cl_n.template get<0>(cell)--;
	}

	/*! \brief Get the number of elements in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return the number of elements in the cell
	 *
	 */
	inline size_t getNelements(const local_index cell_id) const
	{
		return cl_n.template get<0>(cell_id);
	}


	/*! \brief swap to Mem_fast object
	 *
	 * \param mem object to swap the memory with
	 *
	 */
	inline void swap(Mem_fast<Memory,local_index> & mem)
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
	inline void swap(Mem_fast<Memory,local_index> && mem)
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
		{cl_n.template get<0>(i) = 0;}
	}

	/*! \brief Get the first element of a cell (as reference)
	 *
	 * \param cell_id cell-id
	 *
	 * \return a reference to the first element
	 *
	 */
	inline const local_index & getStartId(local_index cell_id) const
	{
		return cl_base.template get<0>(cell_id*slot);
	}

	/*! \brief Get the last element of a cell (as reference)
	 *
	 * \param cell_id cell-id
	 *
	 * \return a reference to the last element
	 *
	 */
	inline const local_index & getStopId(local_index cell_id) const
	{
		return cl_base.template get<0>(cell_id*slot+cl_n.template get<0>(cell_id));
	}

	/*! \brief Just return the value pointed by part_id
	 *
	 * \param part_id
	 *
	 * \return the value pointed by part_id
	 *
	 */
	inline const local_index & get_lin(const local_index * part_id) const
	{
		return *part_id;
	}

public:

	//! expose the type of the local index
	typedef local_index loc_index;

	/*! \brief Constructor
	 *
	 * \param slot number of slot for each cell
	 *
	 */
	inline Mem_fast(local_index slot)
	:slot(slot)
	{}

	/*! \brief Set the number of slot for each cell
	 *
	 * \param number of slot
	 *
	 */
	inline void set_slot(local_index slot)
	{
		this->slot = slot;
	}

	/*! \brief Convert the structure to a structure usable into a kernel
	 *
	 * \return an object usable in the kernel
	 *
	 */
	Mem_fast_ker<Memory,local_index> toKernel()
	{
		Mem_fast_ker<Memory,local_index> mfk(cl_n.toKernel(),cl_base.toKernel());

		return mfk;
	}

};


#endif /* CELLLISTSTANDARD_HPP_ */
