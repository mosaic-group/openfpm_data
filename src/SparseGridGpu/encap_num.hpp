/*
 * encap_num.hpp
 *
 *  Created on: Apr 25, 2020
 *      Author: i-bird
 */

#ifndef ENCAP_NUM_HPP_
#define ENCAP_NUM_HPP_

template<typename encap>
struct enc_num
{
	encap ec;
	unsigned int offset;

	__device__ enc_num(const encap & ec, unsigned int offset)
	:ec(ec),offset(offset)
	{}

	template<unsigned int p>
	__device__ auto get() -> decltype(ec.template get<p>()[offset])
	{
		return ec.template get<p>()[offset];
	}

	__device__ unsigned int getOffset()
	{
		return offset;
	}
};


#endif /* ENCAP_NUM_HPP_ */
