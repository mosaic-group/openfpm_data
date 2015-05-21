/*
 * CellListStandard.hpp
 *
 *  Created on: Mar 22, 2015
 *      Author: Pietro Incardona
 */

#ifndef CELLLISTSTANDARD_HPP_
#define CELLLISTSTANDARD_HPP_

#include "CellDecomposer.hpp"
#include "Space/SpaceBox.hpp"
#include "mathutil.hpp"
#include "CellNNIterator.hpp"
#include "Space/Shape/HyperCube.hpp"

// Compile time array functor needed to generate array at compile-time of type
// {0,0,0,0,0,.....}
// {3,3,3,3,3,3,.....}

 template<size_t index, size_t N> struct Fill_three {
    enum { value = 3 };
 };

 template<size_t index, size_t N> struct Fill_zero {
    enum { value = 0 };
 };

 template<size_t index, size_t N> struct Fill_two {
    enum { value = 2 };
 };

 template<size_t index, size_t N> struct Fill_one {
    enum { value = 1 };
 };

/*! \brief Class for STANDARD cell list implementation
 *
 * This class implement the STANDARD cell list, fast but memory
 * expensive. The memory allocation is (M * N_cell_max)*sizeof(ele) + M*8
 *
 * M = number of cells
 * N_cell_max = maximum number of elements in a cell
 *
 * \note Because N_cell_max >= N/M then M * N_cell_max >= O(N)
 *
 * \warning Not not use for high asymmetric distribution
 *
 * \tparam dim Dimansionality of the space
 * \tparam T type of the space float, double, complex
 * \tparam base Base structure that store the information
 *
 *
 */
template<unsigned int dim, typename T, typename base>
class CellList<dim,T,FAST,base> : public CellDecomposer_sm<dim,T>
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

	// Number of slot for each cell
	size_t slot;

	// number of particle in each cell list
	openfpm::vector<size_t> cl_n;

	// elements that each cell store (each cell can store a number
	// of elements == slot )
	base cl_base;

	//Origin point
	Point<dim,T> orig;

	void realloc()
	{
		// we do not have enough slots reallocate the basic structure with more
		// slots

		// Create a cell-list with double of the slots

		CellList cl_tmp(this->box,this->gr_cell.getSize(),orig,slot*2);

		// copy cl_base

		for (size_t i = 0 ; i < cl_n.size() ; i++)
		{
			cl_tmp.cl_n.get(i) = cl_n.get(i);
			for (size_t j = 0 ; j < cl_n.get(i) ; j++)
				cl_tmp.cl_base.get(2*i*slot + j) = cl_base.get(slot * i + j);
		}

		// swap the memory
		swap(cl_tmp);
	}

public:

	// Object type that the structure store
	typedef T value_type;

	/*! \brief Return the underlying grid information of the cell list
	 *
	 * \return the grid infos
	 *
	 */
	grid_sm<dim,void> & getGrid()
	{
		CellDecomposer_sm<dim,T>::getGrid();
	}

	/*! Initialize the cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param origin of the Cell list (UNUSED)
	 * \param padding cell
	 * \param slot maximum number of slot
	 *
	 */
	void Initialize(Box<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, const size_t pad = 1, size_t slot=16)
	{
		SpaceBox<dim,T> sbox;
		Initialize(sbox,div,orig,slot);
	}

	/*! Initialize the cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param origin of the Cell list (UNUSED)
	 * \param padding cell
	 * \param slot maximum number of slot
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
	 */
	void Initialize(SpaceBox<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, const size_t pad = 1, size_t slot=16)
	{
		// Add padding
		size_t div_pad[dim];
		for (size_t i = 0 ; i < dim ; i++)
			div_pad[i] = div[i] + 2;

		CellDecomposer_sm<dim,T>::setDimensions(box,div_pad, pad);
		this->slot = slot;
		this->orig = orig;

		// create the array that store the number of particle on each cell and se it to 0

		cl_n.resize(this->tot_n_cell);
		cl_n.fill(0);

		// create the array that store the cell id

		cl_base.resize(this->tot_n_cell * slot);

		// Calculate the NNc_full array, it is a structure to get the neighborhood array

		// compile-time array {0,0,0,....} and {3,3,3,...}

		typedef typename generate_array<size_t,dim, Fill_zero>::result NNzero;
		typedef typename generate_array<size_t,dim, Fill_two>::result NNtwo;
		typedef typename generate_array<size_t,dim, Fill_one>::result NNone;

		// Generate the sub-grid iterator

		grid_key_dx_iterator_sub<dim> gr_sub3(this->gr_cell,NNzero::data,NNtwo::data);

		// Calculate the NNc array

		size_t middle = this->gr_cell.LinId(NNone::data);
		size_t i = 0;
		while (gr_sub3.isNext())
		{
			NNc_full[i] = (long int)this->gr_cell.LinId(gr_sub3.get()) - middle;

			++gr_sub3;
			i++;
		}

		// Calculate the NNc_sym array

		i = 0;
		gr_sub3.reset();
		while (gr_sub3.isNext())
		{
			auto key = gr_sub3.get();

			size_t lin = this->gr_cell.LinId(key);

			// Only the first half is considered
			if (lin < middle)
			{
				++gr_sub3;
				continue;
			}

			NNc_sym[i] = lin - middle;

			++gr_sub3;
			i++;
		}

		// Calculate the NNc_cross array

		i = 0;
		grid_key_dx_iterator_sub<dim> gr_sub2(this->gr_cell,NNzero::data,NNone::data);

		while (gr_sub2.isNext())
		{
			auto key = gr_sub2.get();

			NNc_cr[i] = (long int)this->gr_cell.LinId(key);

			++gr_sub2;
			i++;
		}
	}

	/*! \brief Default constructor
	 *
	 */
	CellList()
	{
	}


	/*! \brief Cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param origin of the Cell list (UNUSED)
	 * \param div grid size on each dimension
	 * \param pad Cell padding
	 * \param slot maximum number of slot
	 *
	 */
	CellList(Box<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, const size_t pad = 1, size_t slot=16)
	{
		SpaceBox<dim,T> sbox(box);
		Initialize(sbox,div,orig,pad,slot);
	}

	/*! \brief Cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param origin of the Cell list (UNUSED)
	 * \param div grid size on each dimension
	 * \param pad Cell padding
	 * \param slot maximum number of slot
	 */
	CellList(SpaceBox<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, const size_t pad = 1, size_t slot=16)
	{
		Initialize(box,div,orig,pad,slot);
	}


	/*! \brief Add to the cell
	 *
	 * \param cell_id Cell id where to add
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

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void add(const T (& pos)[dim], typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = this->getCell(pos,1);

		// add the element to the cell

		addCell(cell_id,ele);
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void add(const Point<dim,T> & pos, typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = this->getCell(pos);

		// add the element to the cell

		addCell(cell_id,ele);
	}

	/*! \brief remove an element from the cell
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 */
	inline void remove(size_t cell, size_t ele)
	{
		cl_n.get(cell)--;
	}

	/*! \brief Return the number of element in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return number of elements in the cell
	 *
	 */
	inline size_t getNelements(const size_t cell_id) const
	{
		return cl_n.get(cell_id);
	}

	/*! \brief Get an element in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline auto get(size_t cell, size_t ele) -> decltype(cl_base.get(cell * slot + ele))
	{
		return cl_base.get(cell * slot + ele);
	}

	/*! \brief Get an element in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	template<unsigned int i> inline auto get(size_t cell, size_t ele) -> decltype(cl_base.get(cell * slot + ele))
	{
		return cl_base.template get<i>(cell * slot + ele);
	}

	/*! \brief Swap the memory
	 *
	 * \param cl Cell list with witch you swap the memory
	 *
	 */
	inline void swap(CellList<dim,T,FAST,base> & cl)
	{
		cl_n.swap(cl.cl_n);
		cl_base.swap(cl.cl_base);
	}

	/*! \brief Get the Cell iterator
	 *
	 * \param return the iterator to the cell
	 *
	 */
	CellIterator<CellList<dim,T,FAST,base>> getIterator(size_t cell)
	{
		return CellIterator<CellList<dim,T,FAST,base>>(cell,*this);
	}

	/*! \brief Get the Nearest Neighborhood iterator
	 *
	 * \param cell cell id
	 *
	 */
	template<unsigned int impl> inline CellNNIterator<dim,CellList<dim,T,FAST,base>,FULL,impl> getNNIterator(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,FAST,base>,FULL,impl> cln(cell,NNc_full,*this);

		return cln;
	}

	template<unsigned int impl> inline CellNNIterator<dim,CellList<dim,T,FAST,base>,SYM,impl> getNNIteratorSym(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,FAST,base>,SYM,impl> cln(cell,NNc_sym,*this);

		return cln;
	}

	template<unsigned int impl> inline CellNNIterator<dim,CellList<dim,T,FAST,base>,CRS,impl> getNNIteratorCross(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,FAST,base>,CRS,impl> cln(cell,NNc_cr,*this);

		return cln;
	}
};


#endif /* CELLLISTSTANDARD_HPP_ */
