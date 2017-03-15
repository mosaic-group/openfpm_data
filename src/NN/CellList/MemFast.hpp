/*
 * MemFast.hpp
 *
 *  Created on: Mar 22, 2015
 *      Author: Pietro Incardona
 */

#ifndef CELLLISTSTANDARD_HPP_
#define CELLLISTSTANDARD_HPP_

#include "CellList.hpp"
#include "Space/SpaceBox.hpp"
#include "util/mathutil.hpp"
#include "Space/Shape/HyperCube.hpp"
#include "CellListIterator.hpp"
#include <unordered_map>
#include "util/common.hpp"



template<unsigned int dim, typename T>
class Mem_fast
{
	//! Number of slot for each cell
	size_t slot;

	//! number of particle in each cell list
	openfpm::vector<size_t> cl_n;

	typedef typename openfpm::vector<size_t> base;

	//! elements that each cell store (each cell can store a number
	//! of elements == slot )
	base cl_base;

	/*! \brief realloc the data structures
	 *
	 *
	 */
	void realloc()
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

	void init_to_zero(size_t slot, size_t tot_n_cell)
	{
		this->slot = slot;

		// create the array that store the number of particle on each cell and se it to 0

		cl_n.resize(tot_n_cell);
		cl_n.fill(0);

		// create the array that store the cell id

		cl_base.resize(tot_n_cell * slot);
	}

	void operator=(const Mem_fast & mem)
	{
		slot = mem.slot;

		cl_n = mem.cl_n;
		cl_base = mem.cl_base;
	}

	void addCell(size_t cell_id, typename base::value_type ele)
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

	void add(size_t cell_id, typename base::value_type ele)
	{
		// add the element to the cell

		this->addCell(cell_id,ele);
	}

	auto get(size_t cell, size_t ele) -> decltype(cl_base.get(cell * slot + ele)) &
	{
		return cl_base.get(cell * slot + ele);
	}

	void remove(size_t cell, size_t ele)
	{
		cl_n.get(cell)--;
	}

	size_t getNelements(const size_t cell_id) const
	{
		return cl_n.get(cell_id);
	}

	void swap(Mem_fast & mem)
	{
		cl_n.swap(mem.cl_n);
		cl_base.swap(mem.cl_base);

		size_t cl_slot_tmp = mem.slot;
		mem.slot = slot;
		slot = cl_slot_tmp;
	}

	void swap(Mem_fast && mem)
	{
		slot = mem.slot;

		cl_n.swap(mem.cl_n);
		cl_base.swap(mem.cl_base);
	}

	void clear()
	{
		for (size_t i = 0 ; i < cl_n.size() ; i++)
			cl_n.get(i) = 0;
	}

	inline const size_t & getStartId(size_t cell_id) const
	{
		return cl_base.get(cell_id*slot);
	}

	inline const size_t & getStopId(size_t cell_id) const
	{
		return cl_base.get(cell_id*slot+cl_n.get(cell_id));
	}

	inline const size_t & get_lin(const size_t * part_id) const
	{
		return *part_id;
	}

public:

	Mem_fast(size_t slot)
	:slot(slot)
	{}

	void set_slot(size_t slot)
	{
		this->slot = slot;
	}
};


#endif /* CELLLISTSTANDARD_HPP_ */
