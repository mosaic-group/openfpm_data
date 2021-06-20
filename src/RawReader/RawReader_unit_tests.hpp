/*
 * RawReader_unit_tests.hpp
 *
 *  Created on: April 16, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_RAW_READER_UNIT_TESTS_HPP_
#define OPENFPM_IO_RAW_READER_UNIT_TESTS_HPP_

#include "RawReader.hpp"

BOOST_AUTO_TEST_SUITE( raw_reader_unit_test )


BOOST_AUTO_TEST_CASE( raw_reader_read_test )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

#ifdef OPENFPM_PDATA

	if (v_cl.rank() != 0) {return;}
	std::string c2 = std::string("openfpm_io/test_data/raw_read_sv_test.bin");

#else

	std::string c2 = std::string("test_data/raw_read_sv_test.bin");

#endif

	grid_cpu<3,aggregate<float,float[3]>> read_bin_test;

	GridRawReader<3,aggregate<float,float[3]>,int> rr;

#ifndef SE_CLASS3

	rr.read(c2,read_bin_test,FORTRAN_STYLE | STRUCT_OF_ARRAY,12);

	auto it = read_bin_test.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		BOOST_REQUIRE_EQUAL(read_bin_test.template get<0>(key),1.5f);

		BOOST_REQUIRE_EQUAL(read_bin_test.template get<1>(key)[0],1.5f);
		BOOST_REQUIRE_EQUAL(read_bin_test.template get<1>(key)[1],2.5f);
		BOOST_REQUIRE_EQUAL(read_bin_test.template get<1>(key)[2],3.5f);

		++it;
	}

#endif
}



BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_PLOT_PLOT_UNIT_TESTS_HPP_ */
