/*
 * CellList_util.hpp
 *
 *  Created on: Oct 20, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_

#define CL_SYMMETRIC 1
#define CL_NON_SYMMETRIC 2

#include "Vector/map_vector.hpp"

/*! \brief populate the Cell-list with particles non symmetric case
 *
 * \tparam dim dimensionality of the space
 * \tparam T type of the space
 * \tparam CellList type of cell-list
 *
 * \param pos vector of positions
 * \param cli Cell-list
 * \param g_m marker (particle below this marker must be inside the domain, particles outside this marker must be outside the domain)
 *
 */
template<unsigned int dim, typename T, typename CellList> void populate_cell_list_no_sym(openfpm::vector<Point<dim,T>> & pos, CellList & cli, size_t g_m)
{
	cli.clear();

	for (size_t i = 0; i < pos.size() ; i++)
	{
		cli.add(pos.get(i), i);
	}
}

/*! \brief populate the Cell-list with particles symmetric case
 *
 * \tparam dim dimensionality of the space
 * \tparam T type of the space
 * \tparam CellList type of cell-list
 *
 * \param pos vector of positions
 * \param cli Cell-list
 * \param g_m marker (particle below this marker must be inside the domain, particles outside this marker must be outside the domain)
 *
 */
template<unsigned int dim, typename T, typename CellList> void populate_cell_list_sym(openfpm::vector<Point<dim,T>> & pos, CellList & cli, size_t g_m)
{
	cli.clear();

	for (size_t i = 0; i < g_m ; i++)
	{
		cli.addDom(pos.get(i), i);
	}

	for (size_t i = g_m; i < pos.size() ; i++)
	{
		cli.addPad(pos.get(i), i);
	}
}

/*! \brief populate the Cell-list with particles generic case
 *
 * \tparam dim dimensionality of the space
 * \tparam T type of the space
 * \tparam CellList type of cell-list
 *
 * \param pos vector of positions
 * \param cli Cell-list
 * \param opt option like CL_SYMMETRIC or CL_NON_SYMMETRIC
 * \param g_m marker (particle below this marker must be inside the domain, particles outside this marker must be outside the domain)
 *
 */
template<unsigned int dim, typename T, typename CellList> void populate_cell_list(openfpm::vector<Point<dim,T>> & pos, CellList & cli, size_t g_m, size_t opt)
{
	if (opt == CL_NON_SYMMETRIC)
		populate_cell_list_no_sym(pos,cli,g_m);
	else
		populate_cell_list_sym(pos,cli,g_m);
}

/*! \brief Structure that contain a reference to a vector of particles
 *
 *
 */
template<unsigned int dim, typename T>
struct pos_v
{
	openfpm::vector<Point<dim,T>> & pos;

	pos_v(openfpm::vector<Point<dim,T>> & pos)
	:pos(pos)
	{}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_ */
