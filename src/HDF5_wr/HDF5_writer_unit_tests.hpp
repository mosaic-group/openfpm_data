/*
 * HDF5_writer_unit_test.hpp
 *
 *  Created on: May 1, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_HDF5_WR_HDF5_WRITER_UNIT_TESTS_HPP_
#define OPENFPM_IO_SRC_HDF5_WR_HDF5_WRITER_UNIT_TESTS_HPP_

#include "HDF5_wr.hpp"

#include "hdf5.h"

BOOST_AUTO_TEST_SUITE( vd_hdf5_chckpnt_rstrt_test_io )

// Dimensionality
const size_t dim = 3;

BOOST_AUTO_TEST_CASE( vector_dist_hdf5_save_test )
{
	openfpm::vector<Point<3,float>> vpos;
	openfpm::vector<aggregate<float[dim]>> vprp;

	// Put forces

	for (size_t i = 0 ; i < 1024 ; i++)
	{
		Point<3,float> p;

		p.get(0) = i;
		p.get(1) = i+13;
		p.get(2) = i+17;

		vpos.add(p);

		vprp.add();
		vprp.template get<0>(vprp.size()-1)[0] = p.get(0) + 100.0;
		vprp.template get<0>(vprp.size()-1)[1] = p.get(1) + 200.0;
		vprp.template get<0>(vprp.size()-1)[2] = p.get(2) + 300.0;
	}

	HDF5_writer<VECTOR_DIST> h5;

	// Save the vector
    h5.save("vector_dist.h5",vpos,vprp);

    HDF5_reader<VECTOR_DIST> h5r;

	openfpm::vector<Point<3,float>> vpos2;
	openfpm::vector<aggregate<float[dim]>> vprp2;

    size_t g_m = 0;
    h5r.load("vector_dist.h5",vpos2,vprp2,g_m);

    BOOST_REQUIRE_EQUAL(1024ul,vpos2.size());
    BOOST_REQUIRE_EQUAL(1024ul,vprp2.size());

    BOOST_REQUIRE_EQUAL(1024ul,g_m);

    // Check that vpos == vpos2 and vprp2 == vprp2

    bool check = true;
    for (size_t i = 0 ; i < vpos.size() ; i++)
    {
    	check &= (vpos.get(i) == vpos2.get(i));
    	check &= (vprp.get_o(i) == vprp2.get_o(i));
    }

    BOOST_REQUIRE_EQUAL(check,true);
}



BOOST_AUTO_TEST_CASE( vector_dist_hdf5_load_test )
{
	Vcluster<> & v_cl = create_vcluster();

#ifdef OPENFPM_PDATA

	std::string c2 = std::string("openfpm_io/test_data/vector_dist_24.h5");

#else

	std::string c2 = std::string("test_data/vector_dist_24.h5");

#endif

	openfpm::vector<Point<3,float>> vpos;
	openfpm::vector<aggregate<float[dim]>> vprp;

	HDF5_reader<VECTOR_DIST> h5;

	size_t g_m = 0;

	// Load the vector
    h5.load(c2,vpos,vprp,g_m);

	/////////////////// Checking data ///////////////////////

	// Check total number of particles
	size_t n_part = vpos.size();
	v_cl.sum(n_part);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(n_part,1024ul*24ul);

	BOOST_REQUIRE_EQUAL(vpos.size(),vprp.size());

	bool check = true;

	for (size_t i = 0 ; i < vpos.size() ; i++)
	{
		check &= (vprp.template get<0>(i)[0] == vpos.template get<0>(i)[0] + 100.0);
		check &= (vprp.template get<0>(i)[1] == vpos.template get<0>(i)[1] + 200.0);
		check &= (vprp.template get<0>(i)[2] == vpos.template get<0>(i)[2] + 300.0);
	}


}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_IO_SRC_HDF5_WR_HDF5_WRITER_UNIT_TESTS_HPP_ */
