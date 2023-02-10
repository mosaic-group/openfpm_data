/*
 * HDF5_writer.hpp
 *
 *  Created on: May 1, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_HDF5_WR_HDF5_WRITER_HPP_
#define OPENFPM_IO_SRC_HDF5_WR_HDF5_WRITER_HPP_


#include "VCluster/VCluster.hpp"
#include "hdf5.h"

template <unsigned int type>
class HDF5_writer
{

	void save()
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " Error: we do not know how to write this type of data" << std::endl;
	}

};

#include "HDF5_writer_vd.hpp"
#include "HDF5_writer_gd.hpp"

#endif /* OPENFPM_IO_SRC_HDF5_WR_HDF5_WRITER_HPP_ */
