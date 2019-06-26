/*
 * SparseGridGpu_iterator_sub.hpp
 *
 *  Created on: Jun 24, 2019
 *      Author: i-bird
 */

#ifndef SPARSEGRIDGPU_ITERATOR_SUB_HPP_
#define SPARSEGRIDGPU_ITERATOR_SUB_HPP_

template<unsigned int dim/*, typename blockMap_type*/>
class SparseGridGpu_iterator_sub
{
	unsigned int chunk;
	unsigned int pnt;

	Box<dim,size_t> sub;

//	blockMap_type & blockMap;

public:

	/*! \brief Reinitialize the iterator
	 *
	 * it re-initialize the iterator with the passed grid_key_dx_iterator_sub
	 * the actual position of the grid_key_dx_iterator_sub is ignored
	 *
	 * \param g_s_it grid_key_dx_iterator_sub
	 *
	 */
	inline void reinitialize(const SparseGridGpu_iterator_sub & g_s_it)
	{
		// Reinitialize the iterator
		chunk = g_s_it.chunk;
		pnt = g_s_it.pnt;
	}


};

class SparseGridGpu_iterator
{

};


#endif /* SPARSEGRIDGPU_ITERATOR_SUB_HPP_ */
