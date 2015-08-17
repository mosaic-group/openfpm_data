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
#define CELL_REALLOC 16

/* NOTE all the implementations
 *
 * has complexity O(1) in getting the cell id and the elements in a cell
 * but with different constants
 * 
 * So how does it work?
 * 
 * I have a particle and want to know all its interactions with particles in its subdomain.
 * Since I don't want to check against all other particles, I only check the proximity condition for those that are "nearby enough".
 * We could do so by dividing the subdomain into many cells of equal size (maximal cutoff/radius of the particles as the length of each side).
 * Once we accomplished that we only have to check for possible interactions with all the particles that are in neighbor cells of the cell our initially chosen particle was in.
 * This is what the NNIterator does: given some start cell, list all elements that elements in this start cell could maybe interact with.
 * There are multiple ways to define neighbor cells, I can check all neighbor cells (full) or only those that are "further ahaed" (sym).
 * If we know that interactions are symmetric, this makes them unique (and saves us nearly half the computations!)
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
