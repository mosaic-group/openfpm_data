/*
 * CellNNIteratorRuntime.hpp
 *
 *  Created on: Nov 18, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLNNITERATORRUNTIME_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLNNITERATORRUNTIME_HPP_

#include "util/mathutil.hpp"

#define FULL openfpm::math::pow(3,dim)
#define SYM  openfpm::math::pow(3,dim)/2 + 1
#define CRS openfpm::math::pow(2,dim)

#define RUNTIME -1

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
 * \tparam impl implementation specific options NO_CHECK do not do check on access, SAFE do check on access
 *
 */
template<unsigned int dim, typename Cell,unsigned int impl>
class CellNNIterator<dim,Cell,RUNTIME,impl>
{
protected:

	//! actual element id
	const typename Cell::Mem_type_type::loc_index * start_id;

	//! stop id to read the end of the cell
	const typename Cell::Mem_type_type::loc_index * stop_id;

	//! Actual NNc_id;
	size_t NNc_id;

	//! Size of the neighboring cells
	size_t NNc_size;

	//! Center cell, or cell for witch we are searching the NN-cell
	const long int cell;

	//! actual cell id = NNc[NNc_id]+cell stored for performance reason
	size_t cell_id;

	//! Cell list
	Cell & cl;

	//! NN cell id
	const long int * NNc;

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

			start_id = &cl.getStartId(cell_id);
			stop_id = &cl.getStopId(cell_id);
		}
	}

private:


public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param NNc_size size of the neighborhood
	 * \param cl Cell structure
	 *
	 */
	inline CellNNIterator(size_t cell, const long int * NNc, size_t NNc_size, Cell & cl)
	:NNc_id(0),NNc_size(NNc_size),cell(cell),cell_id(NNc[NNc_id] + cell),cl(cl),NNc(NNc)
	{
		start_id = &cl.getStartId(cell_id);
		stop_id = &cl.getStopId(cell_id);
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
	inline const typename Cell::Mem_type_type::loc_index & get()
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
template<unsigned int dim, typename Cell,typename vector_pos_type, unsigned int impl>
class CellNNIteratorSym<dim,Cell,vector_pos_type,RUNTIME,impl> : public CellNNIterator<dim,Cell,RUNTIME,impl>
{
	//! index of the particle p
	size_t p;

	//! Position of the particle p
	const vector_pos_type & v;

	/*! Select the next valid element
	 *
	 */
	inline void selectValid()
	{
		if (this->NNc[this->NNc_id] == 0)
		{
			while (this->start_id < this->stop_id)
			{
				size_t q = this->cl.get_lin(this->start_id);
				for (long int i = dim-1 ; i >= 0 ; i--)
				{
					if (v.template get<0>(p)[i] < v.template get<0>(q)[i])
						return;
					else if (v.template get<0>(p)[i] > v.template get<0>(q)[i])
						goto next;
				}
				if (q >= p)	return;
next:
				this->start_id++;
			}

			CellNNIterator<dim,Cell,RUNTIME,impl>::selectValid();
		}
		else
		{
			CellNNIterator<dim,Cell,RUNTIME,impl>::selectValid();
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
	inline CellNNIteratorSym(size_t cell,
			                 size_t p,
							 const long int * NNc,
							 size_t NNc_size,
							 Cell & cl,
							 const vector_pos_type & v)
	:CellNNIterator<dim,Cell,RUNTIME,impl>(cell,NNc,NNc_size,cl),p(p),v(v)
	{
		if (this->NNc_id >= this->NNc_size)
			return;

		selectValid();
	}


	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	inline CellNNIteratorSym<dim,Cell,vector_pos_type,RUNTIME,impl> & operator++()
	{
		this->start_id++;

		selectValid();

		return *this;
	}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLNNITERATORRUNTIME_HPP_ */
