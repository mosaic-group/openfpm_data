/*
 * CellList_util.hpp
 *
 *  Created on: Oct 20, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_

#include "util/ofp_context.hpp"

// Cell list config options
constexpr int CL_SYMMETRIC = 1;
constexpr int CL_NON_SYMMETRIC = 2;
constexpr int CL_LOCAL_SYMMETRIC = 8;
constexpr int CL_LINEAR_CELL_KEYS = 16;
constexpr int CL_HILBERT_CELL_KEYS = 32;
constexpr int CL_GPU_REORDER_POSITION = 64;
constexpr int CL_GPU_REORDER_PROPERTY = 128;
constexpr int CL_GPU_REORDER = CL_GPU_REORDER_POSITION | CL_GPU_REORDER_PROPERTY;

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
