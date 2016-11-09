/*
 * VerletListM.hpp
 *
 *  Created on: Oct 15, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTM_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTM_HPP_

#include "NN/VerletList/VerletNNIteratorM.hpp"

/*! \brief Class for Verlet list implementation with Multiphase
 *
 * \tparam dim Dimensionality of the space
 * \tparam T type of the space float, double ...
 * \tparam CellListImpl Base structure that store the information
 *
 */
template<unsigned int dim, typename T, unsigned int sh_byte , typename transform = shift<dim,T>, typename VerletBase=VerletList<dim,T,FAST,transform> >
class VerletListM : public VerletBase
{

	//! Mask to get the high bits of a number
	typedef boost::high_bit_mask_t<sh_byte>  mask_high;

	//! Mask to get the low bits of a number
	typedef boost::low_bits_mask_t<sizeof(size_t)*8-sh_byte>  mask_low;

public:

	/*! \brief Get the element-id in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param part part id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline size_t getP(size_t part, size_t ele) const
	{
		return VerletBase::get(part,ele) & mask_low::sig_bits_fast;
	}

	/*! \brief Get the element vector in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param part particle id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline size_t getV(size_t part, size_t ele) const
	{
		return (VerletBase::get(part,ele)) >> (sizeof(size_t)*8-sh_byte);
	}

	/*! \brief Get the Neighborhood iterator
	 *
	 * It iterate across all the neighborhood particles of a selected particle
	 *
	 * \param part_id particle id
	 *
	 * \return an interator across the neighborhood particles
	 *
	 */
	template<unsigned int impl=NO_CHECK> inline VerletNNIteratorM<dim,VerletListM<dim,T,sh_byte,transform,VerletBase>,sh_byte> getNNIterator(size_t part_id)
	{
		VerletNNIteratorM<dim,VerletListM<dim,T,sh_byte,transform,VerletBase>,sh_byte> vln(part_id,*this);

		return vln;
	}
};


#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTM_HPP_ */
