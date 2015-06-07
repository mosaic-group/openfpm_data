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
 * The memory allocation is (in byte) Size = O(N*sizeof(ele))
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
template<unsigned int dim, typename T, typename transform, typename base>
class CellList<dim,T,MEMORY,transform,base> : public CellDecomposer_sm<dim,T,transform>
{
	// The array contain the neighborhood of the cell-id in case of asymmetric interaction
	//
	//    * * *
	//    * x *
	//    * * *

	long int NNc_full[openfpm::math::pow(3,dim)];

	// The array contain the neighborhood of the cell-id in case of symmetric interaction
	//
	//   * * *
	//     x *
	//
	long int NNc_sym[openfpm::math::pow(3,dim)/2+1];

	// The array contain the neighborhood of the cell-id in case of symmetric interaction (Optimized)
	//
	//   * *
	//   x *
	//
	long int NNc_cr[openfpm::math::pow(2,dim)];

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

	// Object type that the structure store
	typedef T value_type;

	/*! \brief Cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param org of the Cell list
	 * \param div grid size on each dimension
	 *
	 */
	CellList(SpaceBox<dim,T> & box, size_t (&div)[dim], Point<dim,T> & org)
	:box(box),gr_cell(div)
	{
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the
	 * \param ele element to store
	 *
	 */
	void add(T (& pos)[dim], typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = getCell(pos);

		// Get the number of element the cell is storing

		cl_base[cell_id]->add(ele);
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	void add(const Point<dim,T> & pos, typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = this->getCell(pos);

		// add a new element

		cl_base[cell_id]->add(ele);
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
		return cl_base[cell_id]->size();
	}

	/*! \brief remove an element from the cell
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 */
	void remove(size_t cell, size_t ele)
	{
		cl_base[cell]->remove(ele);
	}

	/*! \brief Get an element in the cell
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 */
	auto get(size_t cell, size_t ele) -> decltype(cl_base[cell]->get(ele))
	{
		return cl_base[cell]->get(ele);
	}

	/*! \brief Get the Nearest Neighborhood iterator
	 *
	 * \param cell cell id
	 *
	 */
	template<unsigned int impl> CellNNIterator<dim,CellList<dim,T,MEMORY,base>,FULL,impl> getNNIterator(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,MEMORY,base>,FULL,impl> cln(cell,NNc_full,*this);

		return cln;
	}

	template<unsigned int impl> CellNNIterator<dim,CellList<dim,T,MEMORY,base>,SYM,impl> getNNIteratorSym(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,MEMORY,base>,SYM,impl> cln(cell,NNc_sym,*this);

		return cln;
	}

	template<unsigned int impl> CellNNIterator<dim,CellList<dim,T,BALANCED,base>,CRS,impl> getNNIteratorCross(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,MEMORY,base>,CRS,impl> cln(cell,NNc_cr,*this);

		return cln;
	}
};


#endif /* CELLISTMEM_HPP_ */
