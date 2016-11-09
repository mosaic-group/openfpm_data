/*
 * CellListStandard.hpp
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

//! Wrapper of the unordered map
template<typename key,typename val>
class wrap_unordered_map: public std::unordered_map<key,val>
{
};

#ifdef HAVE_LIBQUADMATH

#include <boost/multiprecision/float128.hpp>


//! Wrapper of the unordered map
template<typename val>
class wrap_unordered_map<boost::multiprecision::float128,val>
{
};

#endif

template<unsigned int dim, typename T, typename transform = no_transform<dim,T>, typename base=openfpm::vector<size_t>>
class Mem_fast
{
	//! Number of slot for each cell
	size_t slot;

	//! number of particle in each cell list
	openfpm::vector<size_t> cl_n;

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

	void swap_cl(CellList<dim,T,Mem_fast<dim,T,transform,base>,transform,base> && cell)
	{
		slot = cell.slot;

		cl_n.swap(cell.cl_n);
		cl_base.swap(cell.cl_base);
	}

	void equal_cl(const CellList<dim,T,Mem_fast<dim,T,transform,base>,transform,base> & cell)
	{
		slot = cell.slot;

		cl_n = cell.cl_n;
		cl_base = cell.cl_base;
	}

	void cell_add(size_t cell_id, typename base::value_type ele)
	{
		// Get the number of element the cell is storing

		size_t nl = getNele(cell_id);

		if (nl + 1 >= slot)
		{
			realloc();
		}

		// we have enough slot to store another neighbor element

		cl_base.get(slot * cell_id + cl_n.get(cell_id)) = ele;
		cl_n.get(cell_id)++;
	}

	void add_ele(const Point<dim,T> & pos, typename base::value_type ele)
	{

		CellList<dim,T,Mem_fast<dim,T,transform,base>,transform,base> cl;

		// calculate the Cell id

		size_t cell_id = cl.getCell(pos);

		// add the element to the cell

		cl.addCell(cell_id,ele);
	}

	void add_ele(const T (& pos)[dim], typename base::value_type ele)
	{
		CellList<dim,T,Mem_fast<dim,T,transform,base>,transform,base> cl;

		// calculate the Cell id

		size_t cell_id = cl.getCell(pos);

		// add the element to the cell

		cl.addCell(cell_id,ele);
	}

	auto get_ele(size_t cell, size_t ele) -> decltype(cl_base.get(cell * slot + ele)) &
	{
		return cl_base.get(cell * slot + ele);
	}

	void rmv(size_t cell, size_t ele)
	{
		cl_n.get(cell)--;
	}

	size_t getNele(const size_t cell_id) const
	{
		return cl_n.get(cell_id);
	}

	void swap_mem(CellList<dim,T,Mem_fast<dim,T,transform,base>,transform,base> & cl)
	{
		cl_n.swap(cl.cl_n);
		cl_base.swap(cl.cl_base);

		size_t cl_slot_tmp = cl.slot;
		cl.slot = slot;
		slot = cl_slot_tmp;
	}

	void clr()
	{
		for (size_t i = 0 ; i < cl_n.size() ; i++)
			cl_n.get(i) = 0;
	}

	inline size_t getStrtId(size_t cell_id)
	{
		return cell_id*slot;
	}

	inline size_t getStpId(size_t cell_id)
	{
		return cell_id*slot+cl_n.get(cell_id);
	}

	inline size_t & get_neighb(size_t part_id)
	{
		return cl_base.get(part_id);
	}

public:

	Mem_fast(size_t slot)
	:slot(slot)
	{}
};


#endif /* CELLLISTSTANDARD_HPP_ */
