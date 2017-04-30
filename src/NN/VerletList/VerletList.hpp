/*
 * VerletList.hpp
 *
 *  Created on: Aug 16, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_

#include "Vector/map_vector.hpp"
#include "NN/CellList/CellList.hpp"

#ifdef LOCAL_INDEX64
typedef size_t local_index_;
#else
typedef unsigned int local_index_;
#endif

#define FAST 1


/*! \brief Cell list structure
 *
 * Stub object see spacialization
 *
 * \see CellList<dim,T,FAST,transform,base>
 *
 */
template<unsigned int dim,
         typename T,
		 unsigned int impl=FAST,
		 typename transform = no_transform<dim,T>,
		 typename local_index=local_index_,
		 typename CellListImpl = CellList<dim,T,Mem_fast,transform> >
class VerletList
{
};

#include "VerletListFast.hpp"


#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_ */
