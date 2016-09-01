/*
 * CellNNIterator.hpp
 *
 *  Created on: Mar 26, 2015
 *      Author: i-bird
 */

#ifndef CELLNNITERATOR_FULL_HPP_
#define CELLNNITERATOR_FULL_HPP_

#include "util/mathutil.hpp"

#define FULL openfpm::math::pow(3,dim)
#define SYM  openfpm::math::pow(3,dim)/2 + 1
#define CRS openfpm::math::pow(2,dim)

#define NO_CHECK 1
#define SAFE 2

/*! \brief Iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells
 *
 * \note to calculate quantities that involve a total reduction (like energies) use the CellIteratorSymRed
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 * \tparam impl implementation specific options NO_CHECK do not do check on access, SAFE do check on access
 *
 */
template<unsigned int dim, typename Cell,unsigned int NNc_size, unsigned int impl>
class CellNNIterator
{
protected:

	//! actual element id
	size_t start_id;

	//! stop id to read the end of the cell
	size_t stop_id;

	//! Actual NNc_id;
	size_t NNc_id;

	//! Center cell, or cell for witch we are searching the NN-cell
	const long int cell;

	//! actual cell id = NNc[NNc_id]+cell stored for performance reason
	size_t cell_id;

	//! Cell list
	Cell & cl;

	/*! \brief Select non-empty cell
	 *
	 */
	inline void selectValid()
	{
		while (start_id == stop_id)
		{
			NNc_id++;

			// No more Cell
			if (NNc_id >= NNc_size) return;

			cell_id = NNc[NNc_id] + cell;

			start_id = cl.getStartId(cell_id);
			stop_id = cl.getStopId(cell_id);
		}
	}

private:


	//! NN cell id
	const long int (& NNc)[NNc_size];

public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 *
	 */
	inline CellNNIterator(size_t cell, const long int (&NNc)[NNc_size], Cell & cl)
	:NNc_id(0),cell(cell),cell_id(NNc[NNc_id] + cell),cl(cl),NNc(NNc)
	{
		start_id = cl.getStartId(cell_id);
		stop_id = cl.getStopId(cell_id);
		selectValid();
	}

	/*! \brief Check if there is the next element
	 *
	 * \return true if there is the next element
	 *
	 */
	inline bool isNext()
	{
		if (NNc_id >= NNc_size)
			return false;
		return true;
	}

	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	inline CellNNIterator & operator++()
	{
		start_id++;

		selectValid();

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline typename Cell::value_type & get()
	{
		return cl.get_lin(start_id);
	}
};


/*! \brief Symmetric iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells.
 *
 * \note if we query the neighborhood of p and q is the neighborhood of p
 *          when we will query the neighborhood of q p is not present. This is
 *          useful to implement formula like \f$  \sum_{q = neighborhood(p) and p <= q} \f$
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 * \tparam impl implementation specific options NO_CHECK do not do check on access, SAFE do check on access
 *
 */
template<unsigned int dim, typename Cell,unsigned int NNc_size, unsigned int impl> class CellNNIteratorSym : public CellNNIterator<dim,Cell,NNc_size,impl>
{
	//! index of the particle p
	size_t p;

	/*! Select the next valid element
	 *
	 */
	inline void selectValid()
	{
		if (this->NNc_id == 0)
		{
			while (this->start_id < this->stop_id)
			{
				if (this->cl.get_lin(this->start_id) >= p)	return;
				this->start_id++;
			}

			CellNNIterator<dim,Cell,NNc_size,impl>::selectValid();
		}
		else
		{
			CellNNIterator<dim,Cell,NNc_size,impl>::selectValid();
		}
	}

public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param p index of the particle from which we are searching the neighborhood particles
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 *
	 */
	inline CellNNIteratorSym(size_t cell, size_t p, const long int (&NNc)[NNc_size], Cell & cl)
	:CellNNIterator<dim,Cell,NNc_size,impl>(cell,NNc,cl),p(p)
	{
		selectValid();
	}


	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	inline CellNNIteratorSym<dim,Cell,NNc_size,impl> & operator++()
	{
		this->start_id++;

		selectValid();

		return *this;
	}
};


/*! \brief it iterate through the elements of a cell
 *
 * In general you do not create this object you get it from the CellList structures
 *
 * \tparam Cell cell type
 *
 */
template<typename Cell> class CellIterator
{
	//! Cell list
	Cell & cl;

	//! actual element id inside the cell
	size_t ele_id;

	//! selected cell
	const long int cell;

public:

	/*! \brief Cell iterator
	 *
	 * \param cell Cell id
	 * \param cl Cell on which iterate
	 *
	 */
	inline CellIterator(const size_t cell, Cell & cl)
	:cl(cl),ele_id(0),cell(cell)
	{
	}

	/*! \brief Check if there is the next element
	 *
	 * \return true if there are still neighborhood particles
	 *
	 */
	inline bool isNext()
	{
		return cl.getNelements(cell) > ele_id;
	}

	/*! \brief take the next neoghborhood particle
	 *
	 * \return itself
	 *
	 */
	inline CellIterator & operator++()
	{
		ele_id++;

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline typename Cell::value_type & get()
	{
		return cl.get(cell,ele_id);
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline const typename Cell::value_type & get() const
	{
		return cl.get(cell,ele_id);
	}
};

#endif /* CELLNNITERATOR_FULL_HPP_ */
