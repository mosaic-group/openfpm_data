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

/*! \brief Calculate parameters for the cell list
 *
 * \param pbox Processor box
 * \param div Division array
 * \param r_cut interation radius or size of each cell
 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost. This parameter says how much must be enlarged
 *
 * \return the processor bounding box
 */
template<unsigned int dim, typename St> static inline void cl_param_calculate(Box<dim, St> & pbox, size_t (&div)[dim], St r_cut, const Ghost<dim, St> & enlarge)
{
	// calculate the parameters of the cell list

	// extend by the ghost
	pbox.enlarge(enlarge);

	// Calculate the division array and the cell box
	for (size_t i = 0; i < dim; i++)
	{
		div[i] = static_cast<size_t>((pbox.getP2().get(i) - pbox.getP1().get(i)) / r_cut);
		div[i]++;
		pbox.setHigh(i,pbox.getLow(i) + div[i]*r_cut);
	}
}

/*! \brief Calculate parameters for the symmetric cell list
 *
 * \param[in] dom Simulation domain
 * \param[output] cd_sm This cell-decomposer is set according to the needed division
 * \param[in] g Ghost dimensions
 * \param[output] pad required padding for the cell-list
 *
 * \return the processor bounding box
 */
template<unsigned int dim, typename St> static inline void cl_param_calculateSym(const Box<dim,St> & dom, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, Ghost<dim,St> g, St r_cut, size_t & pad)
{
	size_t div[dim];

	for (size_t i = 0 ; i < dim ; i++)
		div[i] = (dom.getHigh(i) - dom.getLow(i)) / r_cut;

	g.magnify(1.013);

	// Calculate the maximum padding
	for (size_t i = 0 ; i < dim ; i++)
	{
		size_t tmp = std::ceil(fabs(g.getLow(i)) / r_cut);
		pad = (pad > tmp)?pad:tmp;

		tmp = std::ceil(fabs(g.getHigh(i)) / r_cut);
		pad = (pad > tmp)?pad:tmp;
	}

	cd_sm.setDimensions(dom,div,pad);
}

/*! \brief Calculate parameters for the symmetric cell list
 *
 * \param[in] dom Simulation domain
 * \param[out] cd_sm This cell-decomposer is set according to the needed division
 * \param[in] g Ghost part extension
 * \param[out] pad required padding for the cell-list
 *
 * \return the processor bounding box
 */
template<unsigned int dim, typename St> static inline void cl_param_calculateSym(const Box<dim,St> & dom, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, Ghost<dim,St> g, size_t & pad)
{
	size_t div[dim];

	for (size_t i = 0 ; i < dim ; i++)
		div[i] = (dom.getHigh(i) - dom.getLow(i)) / g.getHigh(i);

	g.magnify(1.013);

	pad = 1;

	cd_sm.setDimensions(dom,div,pad);
}

#include "CellListFast.hpp"
#include "CellListBal.hpp"
#include "CellListMem.hpp"

#endif /* CELLLIST_HPP_ */
