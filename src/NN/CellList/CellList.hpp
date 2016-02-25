/*
 * CellList.hpp
 *
 *  Created on: Mar 21, 2015
 *      Author: Pietro Incardona
 */

#ifndef CELLLIST_HPP_
#define CELLLIST_HPP_

#include "Vector/map_vector.hpp"
#include "CellDecomposer.hpp"

//! Point at witch the cell do a reallocation (it should the the maximum for all the implementations)
#define CELL_REALLOC 16ul

/* NOTE all the implementations
 *
 * has complexity O(1) in getting the cell id and the elements in a cell
 * but with different constants
 *
 */

// Fastest implementation of the Cell list
#define FAST 1
// Balanced memory and speed wise implementation of the cell list
#define BALANCED 2
// Memory wise of the Cell list
#define MEMORY 3


/*! \brief Cell list structure
 *
 * Stub object see spacialization
 *
 * \see CellList<dim,T,FAST,transform,base>
 *
 */
template<unsigned int dim, typename T,  unsigned int impl=FAST, typename transform = no_transform<dim,T>, typename base=openfpm::vector<size_t>>
class CellList
{
};

#include "CellListFast.hpp"
#include "CellListBal.hpp"
#include "CellListMem.hpp"

#endif /* CELLLIST_HPP_ */
