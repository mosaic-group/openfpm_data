/*
 * CellNNIteratorM.hpp
 *
 *  Created on: Jun 23, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLNNITERATORM_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLNNITERATORM_HPP_


#include "util/mathutil.hpp"
#include <boost/integer/integer_mask.hpp>


/*! \brief Iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 * \tparam impl implementation specific options NO_CHECK do not do check on access, SAFE do check on access
 *
 */
template<unsigned int dim, typename Cell, unsigned int sh_byte,unsigned int NNc_size, unsigned int impl> class CellNNIteratorM : public CellNNIterator<dim,Cell,NNc_size,impl>
{
	typedef boost::low_bits_mask_t<sizeof(size_t)*8-sh_byte>  mask_low;

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
	CellNNIteratorM(size_t cell, const long int (&NNc)[NNc_size], Cell & cl)
	:CellNNIterator<dim,Cell,NNc_size,impl>(cell,NNc,cl)
	{}


	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getP()
	{
		return CellNNIterator<dim,Cell,NNc_size,impl>::get() & mask_low::sig_bits_fast;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getV()
	{
		return (CellNNIterator<dim,Cell,NNc_size,impl>::get()) >> (sizeof(size_t)*8-sh_byte);
	}
};

/*! \brief it iterate through the elements of a cell
 *
 * In general you do not create this object you get it from the CellList structures
 *
 * \tparam Cell cell type
 *
 */
template<typename Cell, unsigned int sh_byte> class CellIteratorM : public CellIterator<Cell>
{

	typedef boost::low_bits_mask_t<sizeof(size_t)*8-sh_byte>  mask_low;

public:

	/*! \brief Cell iterator
	 *
	 * \param cell Cell id
	 * \param cl Cell on which iterate
	 *
	 */
	CellIteratorM(const size_t cell, Cell & cl)
	:CellIterator<Cell>(cell,cl)
	{}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getP()
	{
		return CellIterator<Cell>::get() & mask_low::sig_bits_fast;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getV()
	{
		return (CellIterator<Cell>::get()) >> (sizeof(size_t)*8-sh_byte);
	}
};


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLNNITERATORM_HPP_ */
