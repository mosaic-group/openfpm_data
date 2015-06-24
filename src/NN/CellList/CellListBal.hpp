/*
 * CellListBal.hpp
 *
 *  Created on: Mar 22, 2015
 *  Last modified: June 25, 2015
 *      Authors: Pietro Incardona, Yaroslav Zaluzhnyi
 */

#ifndef CELLLISTBAL_HPP_
#define CELLLISTBAL_HPP_

#include "CellDecomposer.hpp"
#include "Space/SpaceBox.hpp"
#include "mathutil.hpp"
#include "CellNNIterator.hpp"
#include "Space/Shape/HyperCube.hpp"

/*! \brief Class for BALANCED cell list implementation
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
template<unsigned int dim, typename T, typename base>
class CellList<dim,T,BALANCED,base>: public CellDecomposer_sm<dim,T>
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
	openfpm::vector<base *> cl_base;

	//Origin point
	Point<dim,T> orig;

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
	 * \param origin of the Cell list
	 * \param div grid size on each dimension
	 *
	 */

	void Initialize(Box<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, const size_t pad = 1)
	{
		SpaceBox<dim,T> sbox;
		Initialize(sbox,div,orig);
	}

	/*! Initialize the cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param origin of the Cell list
	 * \param div grid size on each dimension
	 *
	 */

	void Initialize(SpaceBox<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, const size_t pad = 1)
	{
		// Add padding
		size_t div_pad[dim];
		for (size_t i = 0 ; i < dim ; i++)
			div_pad[i] = div[i] + 2;

		CellDecomposer_sm<dim,T>::setDimensions(box,div_pad, pad);

		this->orig = orig;

		//resize the vector to needed number of cells

		cl_base.resize(this->tot_n_cell);

		//filling a vector with pointers to base-structures
		for (int i = 0; i < this->tot_n_cell; i++)
			cl_base.get(i) =  new base;


		// Calculate the NNc-arrays (for neighborhood):

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
	 * \param origin of the Cell list
	 * \param div grid size on each dimension
	 *
	 */
	CellList(Box<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, const size_t pad = 1)
	{
		SpaceBox<dim,T> sbox(box);
		Initialize(sbox,div,orig,pad);
	}

	/*! \brief Cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param origin of the Cell list
	 * \param div grid size on each dimension
	 *
	 */
	CellList(SpaceBox<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, const size_t pad = 1)
	{
		Initialize(box,div,orig,pad);
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	void add(const T (& pos)[dim], typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = this->getCell(pos);

		// add a new element

		cl_base.get(cell_id)->add(ele);
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

		cl_base.get(cell_id)->add(ele);
	}

	/*! \brief remove an element from the cell
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 */
	void remove(size_t cell, size_t ele)
	{
		cl_base.get(cell)->remove(ele);
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
		return cl_base.get(cell_id)->size();
	}

	/*! \brief Get an element in the cell
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 */
	auto get(size_t cell, size_t ele) -> decltype(cl_base.get(cell)->get(ele))
	{
		return cl_base.get(cell)->get(ele);
	}

	//Three next functions are optional

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
	template<unsigned int i> inline auto get(size_t cell, size_t ele) -> decltype(cl_base.get(cell)->get(ele))
	{
		return cl_base.template get<i>(cell)->get(ele);
	}

	/*! \brief Swap the memory
	 *
	 * \param cl Cell list with witch you swap the memory
	 *
	 */
	void swap(CellList<dim,T,BALANCED,base> & cl)
	{
		cl_base.swap(cl.cl_base);
	}


	/*! \brief Get the Cell iterator
	 *
	 * \param return the iterator to the cell
	 *
	 */
	CellIterator<CellList<dim,T,BALANCED,base>> getIterator(size_t cell)
	{
		return CellIterator<CellList<dim,T,BALANCED,base>>(cell,*this);
	}

	/*! \brief Get the Nearest Neighborhood iterator
	 *
	 * \param cell cell id
	 *
	 */
	template<unsigned int impl> CellNNIterator<dim,CellList<dim,T,BALANCED,base>,FULL,impl> getNNIterator(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,BALANCED,base>,FULL,impl> cln(cell,NNc_full,*this);

		return cln;
	}

	template<unsigned int impl> CellNNIterator<dim,CellList<dim,T,BALANCED,base>,SYM,impl> getNNIteratorSym(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,BALANCED,base>,SYM,impl> cln(cell,NNc_sym,*this);

		return cln;
	}

	template<unsigned int impl> CellNNIterator<dim,CellList<dim,T,BALANCED,base>,CRS,impl> getNNIteratorCross(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,BALANCED,base>,CRS,impl> cln(cell,NNc_cr,*this);

		return cln;
	}
};


#endif /* CELLLISTBAL_HPP_ */
