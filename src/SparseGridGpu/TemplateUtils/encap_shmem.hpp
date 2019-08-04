/*
 * memcpy_shmem.hpp
 *
 *  Created on: Aug 3, 2019
 *      Author: i-bird
 */

#ifndef ENCAP_SHMEM_HPP_
#define ENCAP_SHMEM_HPP_

/*! \brief memcpy it split the copy across threads
 *
 * it assume that the memory to copy is smaller than blockDim.x * 4 byte
 *
 */
template<unsigned int sz>
struct encap_shmem
{
	static const int size = (sz / 4)*4 + (sz % 4 != 0)*4;
	static const int nthr = (sz / 4) + (sz % 4 != 0);

	static void copy(int * src, int * dst)
	{
		if (threadIdx.x < nthr)
		{dst[threadIdx.x] = src[threadIdx.x];}
	}
};


#endif /* MEMCPY_SHMEM_HPP_ */
