/*
 * HDF5_loader.hpp
 *
 *  Created on: May 1, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_HDF5_WR_HDF5_READER_HPP_
#define OPENFPM_IO_SRC_HDF5_WR_HDF5_READER_HPP_

#include "VCluster/VCluster.hpp"

template <unsigned int type>
class HDF5_reader
{
	void load()
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " Error: we do not know how to write this type of data" << std::endl;
	}
};

#include "HDF5_reader_vd.hpp"
#include "HDF5_reader_gd.hpp"

#endif /* OPENFPM_IO_SRC_HDF5_WR_HDF5_READER_HPP_ */
