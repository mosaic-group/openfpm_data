/*
 * CellListM.hpp
 *
 *  Created on: Jun 23, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTM_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTM_HPP_

struct PV_cl
{
	size_t ele;
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
 *
 * \tparam CellBase Cell list type for the Multi-Phase cell list implementation
 *
 * ### Declaration of a Multi-Phase cell list and usage
 *
 */
template<unsigned int dim, typename T, unsigned int sh_byte, typename CellBase=CellList<dim,T,FAST,shift<dim, T>> >
class CellListM : public CellBase
{
	typedef boost::high_bit_mask_t<sh_byte>  mask_high;
	typedef boost::low_bit_mask_t<64-sh_byte>  mask_low;

	//! Default Constructor
	CellListM()
	{};

	//! Copy constructor
	CellListM(const CellListM<CellBase> & cell)
	{
		this->operator=(cell);
	}

	//! Copy constructor
	CellListM(CellListM<CellBase> && cell)
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

	/*! \brief Cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param pad Cell padding
	 * \param slot maximum number of slot
	 *
	 */
	CellListM(SpaceBox<dim,T> & box, const size_t (&div)[dim], const size_t pad = 1, size_t slot=STARTING_NSLOT)
	:CellBase(box,div,mat,pad,slot)
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
	 *
	 */
	inline void addCell(size_t cell_id, size_t ele, size_t v_id)
	{
		size_t ele_k = ele | v_id << sizeof(size_t)*8-sh_byte;

		CellBase::addCell(cell_id,ele_k);
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void add(const T (& pos)[dim], size_t ele, size_t v_id)
	{
		// calculate the Cell id

		size_t cell_id = this->getCell(pos);
		size_t ele_k = ele | v_id << sizeof(size_t)*8-sh_byte;

		// add the element to the cell

		CellBase::addCell(cell_id,ele_k);
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void add(const Point<dim,T> & pos, size_t ele, size_t v_id)
	{
		// calculate the Cell id

		size_t cell_id = this->getCell(pos);
		size_t ele_k = ele | v_id << sizeof(size_t)*8-sh_byte;

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
		return CellBase::get(cell * slot + ele) & mask_low;
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
		return (cl_base.get(cell * slot + ele) & mask_high) >> 64-sh_byte;
	}

	/*! \brief Swap the memory
	 *
	 * \param cl Cell list with witch you swap the memory
	 *
	 */
	inline void swap(CellListM<dim,T,FAST,transform,base> & cl)
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
	CellIterator<CellList<dim,T,FAST,transform,base>> getIterator(size_t cell)
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
	 */
	template<unsigned int impl=NO_CHECK> inline CellNNIterator<dim,CellList<dim,T,FAST,transform,base>,FULL,impl> getNNIterator(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,FAST,transform,base>,FULL,impl> cln(cell,NNc_full,*this);

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
	 */
	template<unsigned int impl> inline CellNNIterator<dim,CellList<dim,T,FAST,transform,base>,SYM,impl> getNNIteratorSym(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,FAST,transform,base>,SYM,impl> cln(cell,NNc_sym,*this);

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
	 */
	template<unsigned int impl> inline CellNNIterator<dim,CellList<dim,T,FAST,transform,base>,CRS,impl> getNNIteratorCross(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,FAST,transform,base>,CRS,impl> cln(cell,NNc_cr,*this);

		return cln;
	}

	/*! \brief Clear the cell list
	 *
	 */
	void clear()
	{
		slot = STARTING_NSLOT;
		for (size_t i = 0 ; i < cl_n.size() ; i++)
			cl_n.get(i) = 0;
	}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTM_HPP_ */
