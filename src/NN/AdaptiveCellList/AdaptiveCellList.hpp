/*
 * CellListStandard.hpp
 *
 *  Created on: Mar 22, 2015
 *      Author: Pietro Incardona
 */

#ifndef ADAPTIVECELLLIST_HPP_
#define ADAPTIVECELLLIST_HPP_

#include "../CellList/CellList.hpp"
#include "../CellList/CellDecomposer.hpp"
#include "AdaptiveCellListNNIterator.hpp"
#include "Space/SpaceBox.hpp"


// Stub implementation
template<unsigned int dim, typename T,  unsigned int impl=BALANCED, typename base=openfpm::vector<size_t>>
class AdaptiveCellList
{
};

/*! \brief Class for Adaptive cell list implementation
 *
 * This class implement the Adaptive cell list, ...
 *
 * \tparam dim Dimensionality of the space
 * \tparam T type of the space float, double, complex
 * \tparam base basic object
 *
 */
template<unsigned int dim, typename T, typename base>
class AdaptiveCellList<dim,T,BALANCED,base> : public CellDecomposer_sm<dim,T>
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

	//A point
	Point<dim,T> orig;

public:

	// Object type that the structure store
	typedef T value_type;

	/*! Initialize the cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param origin of the Cell list
	 * \param div grid size on each dimension
	 *
	 */
	void Initialize(SpaceBox<dim,T> & sbox, size_t (&div)[dim], Point<dim,T> & orig, size_t slot=16)
	{
		CellDecomposer_sm<dim,T>::setDimensions(sbox,div);
		this->slot = slot;
		this->orig = orig;
	}

	/*! \brief Default constructor
	 *
	 */
	AdaptiveCellList()
	{
	}


	/*! \brief Cell list
	 *
	 * \param ... (Non default constructor if needed)
	 *
	 */
	AdaptiveCellList(Box<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, size_t slot=16)
	{
		SpaceBox<dim,T> sbox(box);
		Initialize(sbox,div,orig,slot);
	}


	/*! \brief Add to the cell an element
	 *
	 * \param cell_id Cell id where to add
	 * \param ele element to add
	 *
	 */
	inline void addCell(size_t cell_id, typename base::value_type ele)
	{
	}

	/*! \brief Add to the cell an element (from points coordinate)
	 *
	 * \param pos array that contain the coordinate +1 (the last is the radius of the particle)
	 * \param ele element to store
	 *
	 */
	inline void add(const T (& pos)[dim+1], typename base::value_type ele)
	{
	}

	/*! \brief Add to the cell an element (from points coordinate)
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void add(const Point<dim+1,T> & pos, typename base::value_type ele)
	{
	}

	/*! \brief remove an element from the cell
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 */
	inline void remove(size_t cell, size_t ele)
	{
	  // i could change params later on
	}

	/*! \brief Return the number of element in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return number of elements in the cell
	 *
	 */
	inline size_t getNelements(size_t cell_id)
	{
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
	inline value_type get(size_t cell, size_t ele)
	{
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
	inline void swap(AdaptiveCellList<dim,T,BALANCED,base> & cl)
	{
		cl_n.swap(cl.cl_n);
		cl_base.swap(cl.cl_base);
	}

	/*! \brief Get the Nearest Neighborhood iterator
	 *
	 * \param cell cell id
	 *
	 */
	template<unsigned int impl> inline AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,base>,FULL,impl> getNNIterator(size_t cell)
	{
		AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,base>,FULL,impl> cln(cell,NNc_full,*this);

		return cln;
	}

	template<unsigned int impl> inline AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,base>,SYM,impl> getNNIteratorSym(size_t cell)
	{
		AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,base>,SYM,impl> cln(cell,NNc_sym,*this);

		return cln;
	}

	template<unsigned int impl> inline AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,base>,CRS,impl> getNNIteratorCross(size_t cell)
	{
		AdaptiveCellNNIterator<dim,AdaptiveCellList<dim,T,BALANCED,base>,CRS,impl> cln(cell,NNc_cr,*this);

		return cln;
	}
};


#endif /* CELLLISTSTANDARD_HPP_ */
