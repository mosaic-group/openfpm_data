#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "HDF5_wr.hpp"

#include "hdf5.h"

BOOST_AUTO_TEST_SUITE( vd_hdf5_chckpnt_rstrt_test_gpu )

// Dimensionality
const size_t dim = 3;

BOOST_AUTO_TEST_CASE( vector_dist_hdf5_save_test_gpu )
{
	openfpm::vector_gpu<Point<3,float>> vpos;
	openfpm::vector_gpu<aggregate<double,float[dim],size_t,float[dim][dim]>> vprp;

	// Put forces

	for (size_t i = 0 ; i < 1024 ; i++)
	{
		Point<3,float> p;

		p.get(0) = i;
		p.get(1) = i+13;
		p.get(2) = i+17;

		vpos.add(p);

		vprp.add();

		vprp.template get<0>(vprp.size()-1) = p.get(0) + 113.0;

		vprp.template get<1>(vprp.size()-1)[0] = p.get(0) + 100.0;
		vprp.template get<1>(vprp.size()-1)[1] = p.get(1) + 200.0;
		vprp.template get<1>(vprp.size()-1)[2] = p.get(2) + 300.0;

		vprp.template get<2>(vprp.size()-1) = p.get(0) + 10013.0;

		vprp.template get<3>(vprp.size()-1)[0][0] = p.get(0) + 600.0;
		vprp.template get<3>(vprp.size()-1)[0][1] = p.get(0) + 900.0;
		vprp.template get<3>(vprp.size()-1)[0][2] = p.get(0) + 1600.0;
		vprp.template get<3>(vprp.size()-1)[1][0] = p.get(0) + 5600.0;
		vprp.template get<3>(vprp.size()-1)[1][1] = p.get(0) + 6900.0;
		vprp.template get<3>(vprp.size()-1)[1][2] = p.get(0) + 7600.0;
		vprp.template get<3>(vprp.size()-1)[2][0] = p.get(0) + 9600.0;
		vprp.template get<3>(vprp.size()-1)[2][1] = p.get(0) + 1900.0;
		vprp.template get<3>(vprp.size()-1)[2][2] = p.get(0) + 101600.0;
	}

	HDF5_writer<VECTOR_DIST> h5;

	// Save the vector
    h5.save("vector_dist.h5",vpos,vprp);

    HDF5_reader<VECTOR_DIST> h5r;

	openfpm::vector_gpu<Point<3,float>> vpos2;
	openfpm::vector_gpu<aggregate<double,float[dim],size_t,float[dim][dim]>> vprp2;

    size_t g_m = 0;
    h5r.load("vector_dist.h5",vpos2,vprp2,g_m);

    BOOST_REQUIRE_EQUAL(1024ul,vpos2.size());
    BOOST_REQUIRE_EQUAL(1024ul,vprp2.size());

    BOOST_REQUIRE_EQUAL(1024ul,g_m);

    // Check that vpos == vpos2 and vprp2 == vprp2

    bool check = true;
    for (size_t i = 0 ; i < vpos.size() ; i++)
    {
    	Point<3,float> p1 = vpos.get(i);
    	Point<3,float> p2 = vpos2.get(i);

    	check &= (p1 == p2);

    	check &= (vprp.template get<1>(i)[0] == vprp2.template get<1>(i)[0]);
    	check &= (vprp.template get<1>(i)[1] == vprp2.template get<1>(i)[1]);
    	check &= (vprp.template get<1>(i)[2] == vprp2.template get<1>(i)[2]);

    	check &= (vprp.template get<0>(i) == vprp2.template get<0>(i));
    	check &= (vprp.template get<2>(i) == vprp2.template get<2>(i));

    	check &= (vprp.template get<3>(i)[0][0] == vprp2.template get<3>(i)[0][0]);
    	check &= (vprp.template get<3>(i)[0][1] == vprp2.template get<3>(i)[0][1]);
    	check &= (vprp.template get<3>(i)[0][2] == vprp2.template get<3>(i)[0][2]);
    	check &= (vprp.template get<3>(i)[1][0] == vprp2.template get<3>(i)[1][0]);
    	check &= (vprp.template get<3>(i)[1][1] == vprp2.template get<3>(i)[1][1]);
    	check &= (vprp.template get<3>(i)[1][2] == vprp2.template get<3>(i)[1][2]);
    	check &= (vprp.template get<3>(i)[2][0] == vprp2.template get<3>(i)[2][0]);
    	check &= (vprp.template get<3>(i)[2][1] == vprp2.template get<3>(i)[2][1]);
    	check &= (vprp.template get<3>(i)[2][2] == vprp2.template get<3>(i)[2][2]);
    }

    BOOST_REQUIRE_EQUAL(check,true);
}

BOOST_AUTO_TEST_SUITE_END()
