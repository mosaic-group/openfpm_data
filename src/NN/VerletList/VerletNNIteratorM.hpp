/*
 * VerletNNIteratorM.hpp
 *
 *  Created on: Oct 15, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETNNITERATORM_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETNNITERATORM_HPP_


/*! \brief Iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Ver Base verlet structure
 * \tparam sh_byte number of bits for the phase number
 *
 */
template<unsigned int dim, typename Ver, unsigned int sh_byte> class VerletNNIteratorM : public VerletNNIterator<dim,Ver>
{
	//! Mask to get the high bits of a number
	typedef boost::high_bit_mask_t<sh_byte>  mask_high;

	//! Mask to get the low bits of a number
	typedef boost::low_bits_mask_t<sizeof(size_t)*8-sh_byte>  mask_low;


public:

	/*! \brief Constructor for the Verlet iterator Multi-phase
	 *
	 *
	 */
	VerletNNIteratorM(size_t part_id, Ver & ver)
	:VerletNNIterator<dim,Ver>(part_id,ver)
	{
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getP()
	{
		return VerletNNIterator<dim,Ver>::get() & mask_low::sig_bits_fast;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline size_t getV()
	{
		return (VerletNNIterator<dim,Ver>::get()) >> (sizeof(size_t)*8-sh_byte);
	}


};


#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETNNITERATORM_HPP_ */
