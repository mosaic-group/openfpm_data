/*
 * MemMemoryWise.hpp
 *
 *  Created on: Mar 22, 2015
 *  Last modified: June 25, 2015
 *      Authors: Pietro Incardona, Yaroslav Zaluzhnyi
 */

#ifndef CELLISTMEM_HPP_
#define CELLISTMEM_HPP_

#include "CellList.hpp"

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
template<unsigned int dim, typename T, typename transform = no_transform<dim,T>, typename base=openfpm::vector<size_t>>
class Mem_mw
{
	// each cell has a dynamic structure
	// that store the elements in the cell
	std::unordered_map<size_t,base> cl_base;

	//Origin point
	Point<dim,T> orig;

public:

	// Object type that the structure store
	typedef T value_type;


	void init_to_zero(size_t slot, size_t tot_n_cell)
	{
		//resize the map to needed number of cells

		cl_base.rehash(tot_n_cell);

		//filling a map with "base" structures
		for (int i = 0; i < tot_n_cell; i++)
		{   base b;
			cl_base.insert(tot_n_cell, b);
		}
	}

	void swap_cl(CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base> && cell)
	{
		cl_base.swap(cell.cl_base);
	}

	void equal_cl(const CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base> & cell)
	{
		cl_base = cell.cl_base;
	}

	void cell_add(size_t cell_id, typename base::value_type ele)
	{
		//add another neighbor element

		cl_base.get(cell_id) = ele;
	}

	void add_ele(const T (& pos)[dim], typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>::getCell(pos);

		// add the element to the cell

		CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>::addCell(cell_id,ele);
	}

	void add_ele(const Point<dim,T> & pos, typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>::getCell(pos);

		// add the element to the cell

		CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>::addCell(cell_id,ele);
	}

	void rmv(size_t cell, size_t ele)
	{
		cl_base.get(cell).remove(ele);
	}

	size_t getNele(const size_t cell_id)
	{
		return cl_base.get(cell_id).size();
	}

	auto get_ele(size_t cell, size_t ele) -> decltype(cl_base.get(cell).get(ele)) &
	{
		return cl_base.get(cell).get(ele);
	}

	void swap_mem(CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base> & cl)
	{
		cl_base.swap(cl.cl_base);
	}

	CellIterator<CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>> getCellIt(size_t cell)
	{
		return CellIterator<CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>>(cell,*this);
	}

	template<unsigned int impl=NO_CHECK> inline CellNNIterator<dim,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>,FULL,impl> getNNIt(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>,FULL,impl> cln(cell,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>::NNc_full,*this);

		return cln;
	}

	template<unsigned int impl=NO_CHECK> inline CellNNIteratorRadius<dim,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>,impl> getNNItRad(size_t cell, T r_cut, openfpm::vector<long int> & NNc)
	{
		if (NNc.size() == 0)
			NNcalc(r_cut,NNc);

		CellNNIteratorRadius<dim,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>,impl> cln(cell,NNc,*this);

		return cln;
	}

	template<unsigned int impl> inline CellNNIteratorSym<dim,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>,SYM,impl> getNNItSym(size_t cell, size_t p)
	{
		CellNNIteratorSym<dim,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>,SYM,impl> cln(cell,p,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>::NNc_sym,*this);

		return cln;
	}

	template<unsigned int impl> inline CellNNIterator<dim,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>,CRS,impl> getNNItCross(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>,CRS,impl> cln(cell,CellList<dim,T,Mem_mw<dim,T,transform,base>,transform,base>::NNc_cr,*this);

		return cln;
	}

	void clr()
	{
		for (size_t i = 0 ; i < cl_base.size() ; i++)
		{
			for (size_t j = 0; j < cl_base.get(i).size(); j++)
				cl_base.get(i).get(j) = 0;
		}
	}

	inline size_t getStrtId(size_t part_id)
	{
		return part_id;
	}

	inline size_t getStpId(size_t part_id)
	{
		return part_id;
	}

	inline size_t & get_neighb(size_t part_id)
	{
		return cl_base.get(part_id);
	}

public:

	Mem_mw()
	{}

};


#endif /* CELLISTMEM_HPP_ */
