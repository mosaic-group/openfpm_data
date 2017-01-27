/*
 * CellNNIteratorRuntimeM.hpp
 *
 *  Created on: Jan 26, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLNNITERATORRUNTIMEM_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLNNITERATORRUNTIMEM_HPP_


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
template<unsigned int dim, typename Cell, unsigned int sh_byte, unsigned int impl>
class CellNNIteratorM<dim,Cell,sh_byte,RUNTIME,impl> : public CellNNIterator<dim,Cell,RUNTIME,impl>
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
	CellNNIteratorM(size_t cell, const long int * NNc, size_t NNc_size, Cell & cl)
	:CellNNIterator<dim,Cell,RUNTIME,impl>(cell,NNc,NNc_size,cl)
	{}


	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getP()
	{
		return CellNNIterator<dim,Cell,RUNTIME,impl>::get() & mask_low::sig_bits_fast;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getV()
	{
		return (CellNNIterator<dim,Cell,RUNTIME,impl>::get()) >> (sizeof(size_t)*8-sh_byte);
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
template<unsigned int dim, typename Cell, unsigned int sh_byte, unsigned int impl>
class CellNNIteratorSymM<dim,Cell,sh_byte,RUNTIME,impl> : public CellNNIterator<dim,Cell,RUNTIME,impl>
{
	typedef boost::low_bits_mask_t<sizeof(size_t)*8-sh_byte>  mask_low;

	//! phase of the particle p
	size_t pp;

	//! index of the particle p
	size_t p;

	const openfpm::vector<Point<dim,typename Cell::stype>> & pos;

	//! Position of the particles in the phases
	const openfpm::vector<pos_v<dim,typename Cell::stype>> & ps;

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
					if (pos.template get<0>(p)[i] < ps.get(q >> (sizeof(size_t)*8-sh_byte)).pos.template get<0>(q & mask_low::sig_bits_fast)[i])
						return;
					else if (pos.template get<0>(p)[i] > ps.get(q >> (sizeof(size_t)*8-sh_byte)).pos.template get<0>(q & mask_low::sig_bits_fast)[i])
						goto next;
				}
				if (q >> (sizeof(size_t)*8-sh_byte) != pp)	return;
				if ((q & mask_low::sig_bits_fast) >= p)	return;
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
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 *
	 */
	CellNNIteratorSymM(size_t cell,
					   size_t pp,
			           size_t p,
					   const long int * NNc,
					   size_t NNc_size,
					   Cell & cl,
					   const openfpm::vector<Point<dim,typename Cell::stype>> & pos,
					   const openfpm::vector<pos_v<dim,typename Cell::stype>> & ps)
	:CellNNIterator<dim,Cell,RUNTIME,impl>(cell,NNc,NNc_size,cl),pp(pp),p(p),pos(pos),ps(ps)
	{}


	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getP()
	{
		return CellNNIterator<dim,Cell,RUNTIME,impl>::get() & mask_low::sig_bits_fast;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getV()
	{
		return (CellNNIterator<dim,Cell,RUNTIME,impl>::get()) >> (sizeof(size_t)*8-sh_byte);
	}

	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	inline CellNNIteratorSymM<dim,Cell,sh_byte,RUNTIME,impl> & operator++()
	{
		this->start_id++;

		selectValid();

		return *this;
	}
};


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLNNITERATORRUNTIMEM_HPP_ */
