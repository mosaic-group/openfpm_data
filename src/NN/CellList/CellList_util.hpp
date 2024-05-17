/*
 * CellList_util.hpp
 *
 *  Created on: Oct 20, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_

#include "util/ofp_context.hpp"

#include "Vector/map_vector.hpp"

/*! \brief Structure that contain a reference to a vector of particles
 *
 *
 */
template<typename vector_pos_type>
struct pos_v
{
	vector_pos_type & pos;

	pos_v(vector_pos_type & pos)
	:pos(pos)
	{}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_ */
