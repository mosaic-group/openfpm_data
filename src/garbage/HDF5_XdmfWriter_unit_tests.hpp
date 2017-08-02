/*
 * H5PartWriter_unit_tests.hpp
 *
 *  Created on: Feb 22, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_HDF5_XDMFWRITER_HDF5_XDMFWRITER_UNIT_TESTS_HPP_
#define OPENFPM_IO_SRC_HDF5_XDMFWRITER_HDF5_XDMFWRITER_UNIT_TESTS_HPP_

#include "VCluster.hpp"
#include "util/SimpleRNG.hpp"
#include "HDF5_XdmfWriter.hpp"

BOOST_AUTO_TEST_SUITE( HDF5_writer_test )


BOOST_AUTO_TEST_CASE( HDF5_writer_use)
{
	openfpm::vector<Point<3,double>> pv;
	openfpm::vector<Point_test<double>> pvp;

	SimpleRNG rng;

	Vcluster & v_cl = *global_v_cluster;

	if (v_cl.getProcessingUnits() != 3)
		return;

	double z_base = v_cl.getProcessUnitID();

	// fill 1000 particles for each processors

	for (size_t i = 0 ; i < 1000 ; i++)
	{
		Point<3,double> p;
		p[0] = rng.GetUniform();
		p[1] = rng.GetUniform();
		p[2] = z_base+rng.GetUniform();

		pv.add(p);

		p[0] += 2.0;

		Point_test<double> pt;
		pt.fill();

		pvp.add(pt);
	}

	HDF5_XdmfWriter<H5_POINTSET> h5p;
	h5p.template write<Point<3,double>,Point_test<double>,0,1,4,5>("h5part.h5",pv,pvp,1000);

	// check that match

	bool test = compare("test_h5part.h5part","test_h5part_test.h5part");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_SUITE_END()



#endif /* OPENFPM_IO_SRC_HDF5_XDMFWRITER_HDF5_XDMFWRITER_UNIT_TESTS_HPP_ */
