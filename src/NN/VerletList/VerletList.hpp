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

#define VERLETLIST_FAST() VerletList<dim,St,Mem_fast<>,shift<dim,St> >
#define VERLETLIST_BAL VerletList<dim,St,Mem_bal,shift<dim,St> >
#define VERLETLIST_MEM VerletList<dim,St,Mem_mem<>,shift<dim,St> >

/*! \brief Cell list structure
 *
 * Stub object see spacialization
 *
 * \see CellList<dim,T,FAST,transform,base>
 *
 */
/*template<unsigned int dim,
         typename T,
		 typename impl=Mem_fast,
		 typename transform = no_transform<dim,T>,
		 typename local_index=local_index_,
		 typename CellListImpl = CellList<dim,T,Mem_fast,transform> >
class VerletList
{
};*/

#include "VerletListFast.hpp"


#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_ */
