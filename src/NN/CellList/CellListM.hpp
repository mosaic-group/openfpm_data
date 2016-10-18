/*
 * CellListM.hpp
 *
 *  Created on: Jun 23, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTM_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTM_HPP_

#include "CellNNIteratorM.hpp"

struct PV_cl
{
	//! particle id
	size_t ele;
	//! phase id
	size_t v;
};

/*! \brief Class for Multi-Phase cell-list
 *
 * This class implement a Multi-Phase cell list. In practice this Cell list can contain
 * the particles from multiple vector distributed. By default this cell list is based on
 * Cell list fast with shifting
 *
 * \verbatim
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
 * \endverbatim
 *
 * \tparam dim dimensionality
 * \tparam T type of the space
 * \tparam sh_byte bit to dedicate to the phases informations
 * \tparam CellBase Base cell list used for the implementation
 *
 * ### Declaration of a Multi-Phase cell list and usage
 *
 */
template<unsigned int dim, typename T, unsigned int sh_byte, typename CellBase=CellList<dim,T,FAST,shift<dim, T>> >
class CellListM : public CellBase
{
	//! Mask to get the high bits of a number
	typedef boost::high_bit_mask_t<sh_byte>  mask_high;

	//! Mask to get the low bits of a number
	typedef boost::low_bits_mask_t<sizeof(size_t)*8-sh_byte>  mask_low;

public:

	//! Default Constructor
	CellListM()
	{};

	//! Copy constructor
	CellListM(const CellListM<dim,T,sh_byte,CellBase> & cell)
	{
		this->operator=(cell);
	}

	//! Copy constructor
	CellListM(CellListM<dim,T,sh_byte,CellBase> && cell)
	{
		this->operator=(cell);
	}


	/*! \brief Cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param mat Matrix transformation
	 * \param pad Cell padding
	 * \param slot maximum number of slot
	 *
	 */
	CellListM(Box<dim,T> & box, const size_t (&div)[dim], Matrix<dim,T> mat, const size_t pad = 1, size_t slot=STARTING_NSLOT)
	:CellBase(box,div,mat,pad,slot)
	{}

	/*! \brief Cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param pad Cell padding
	 * \param slot maximum number of slot
	 *
	 */
	CellListM(Box<dim,T> & box, const size_t (&div)[dim], const size_t pad = 1, size_t slot=STARTING_NSLOT)
	:CellBase(box,div,pad,slot)
	{}


	/*! \brief Destructor
	 *
	 *
	 */
	~CellListM()
	{}

	/*! \brief Add to the cell
	 *
	 * \param cell_id Cell id where to add
	 * \param ele element to add
	 * \param v_id phase id
	 *
	 */
	inline void addCell(size_t cell_id, size_t ele, size_t v_id)
	{
		size_t ele_k = ele | (v_id << (sizeof(size_t)*8-sh_byte));

		CellBase::addCell(cell_id,ele_k);
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 * \param v_id phase id
	 *
	 */
	inline void add(const T (& pos)[dim], size_t ele, size_t v_id)
	{
		// calculate the Cell id

		size_t cell_id = this->getCell(pos);
		size_t ele_k = ele | (v_id << (sizeof(size_t)*8-sh_byte));

		// add the element to the cell

		CellBase::addCell(cell_id,ele_k);
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 * \param v_id phase id
	 *
	 */
	inline void add(const Point<dim,T> & pos, size_t ele, size_t v_id)
	{
		// calculate the Cell id

		size_t cell_id = this->getCell(pos);
		size_t ele_k = ele | (v_id << (sizeof(size_t)*8-sh_byte));

		// add the element to the cell

		CellBase::addCell(cell_id,ele_k);
	}

	/*! \brief Get the element-id in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline size_t getP(size_t cell, size_t ele) const
	{
		return CellBase::get(cell,ele) & mask_low::sig_bits_fast;
	}

	/*! \brief Get the element vector in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline size_t getV(size_t cell, size_t ele) const
	{
		return (CellBase::get(cell,ele)) >> (sizeof(size_t)*8-sh_byte);
	}

	/*! \brief Swap the memory
	 *
	 * \param cl Cell list with witch you swap the memory
	 *
	 */
	inline void swap(CellListM<dim,T,sh_byte,CellBase> & cl)
	{
		CellBase::swap(*this);
	}

	/*! \brief Get the Cell iterator
	 *
	 * \param cell
	 *
	 * \return the iterator to the elements inside cell
	 *
	 */
	CellIterator<CellListM<dim,T,sh_byte,CellBase>> getIterator(size_t cell)
	{
		return CellBase::getIterator(cell);
	}

	/*! \brief Get the Neighborhood iterator
	 *
	 * It iterate across all the element of the selected cell and the near cells
	 *
	 *  \verbatim

	     * * *
	     * x *
	     * * *

	   \endverbatim
	 *
	 * * x is the selected cell
	 * * * are the near cell
	 *
	 * \param cell cell id
	 *
	 * \return an iterator over the particle of the selected cell
	 *
	 */
	template<unsigned int impl=NO_CHECK> inline CellNNIteratorM<dim,CellListM<dim,T,sh_byte,CellBase>,sh_byte,FULL,impl> getNNIterator(size_t cell)
	{
		CellNNIteratorM<dim,CellListM<dim,T,sh_byte,CellBase>,sh_byte,FULL,impl> cln(cell,CellListM<dim,T,sh_byte,CellBase>::NNc_full,*this);

		return cln;
	}


	/*! \brief Get the Neighborhood iterator
	 *
	 * It iterate across all the element of the selected cell and the near cells
	 *
	 *  \verbatim

	   * * *
	     x *

	   \endverbatim
	 *
	 * * x is the selected cell
	 * * * are the near cell
	 *
	 * \param cell cell id
	 *
	 * \return Cell-list structure
	 *
	 */
	template<unsigned int impl> inline CellNNIteratorM<dim,CellListM<dim,T,sh_byte,CellBase>,sh_byte,SYM,impl> getNNIteratorSym(size_t cell)
	{
		CellNNIteratorM<dim,CellListM<dim,T,sh_byte,CellBase>,sh_byte,SYM,impl> cln(cell,CellListM<dim,T,sh_byte,CellBase>::NNc_sym,*this);

		return cln;
	}


	/*! \brief Get the Neighborhood iterator
	 *
	 * It iterate across all the element of the selected cell and the near cells
	 *
	 *  \verbatim

	   * *
	   x *

	   \endverbatim
	 *
	 * * x is the selected cell
	 * * * are the near cell
	 *
	 * \param cell cell id
	 *
	 * \return Cell-list structure
	 *
	 */
	template<unsigned int impl> inline CellNNIteratorM<dim,CellListM<dim,T,sh_byte,CellBase>,sh_byte,CRS,impl> getNNIteratorCross(size_t cell)
	{
		CellNNIteratorM<dim,CellListM<dim,T,sh_byte,CellBase>,sh_byte,CRS,impl> cln(cell,CellListM<dim,T,sh_byte,CellBase>::NNc_cr,*this);

		return cln;
	}


	/*! \brief operator=
	 *
	 * \param clm Cell list to copy
	 *
	 * \return Cell-list structure
	 *
	 */
	CellListM<dim,T,sh_byte,CellBase> & operator=(CellListM<dim,T,sh_byte,CellBase> && clm)
	{
		CellBase::swap(clm);

		return this;
	}

	/*! \brief operator=
	 *
	 * \param clm Cell list to copy
	 *
	 * \return Cell-list structure
	 *
	 */
	CellListM<dim,T,sh_byte,CellBase> & operator=(const CellListM<dim,T,sh_byte,CellBase> & clm)
	{
		static_cast<CellBase *>(this)->operator=(*static_cast<const CellBase *>(&clm));

		return *this;
	}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTM_HPP_ */
