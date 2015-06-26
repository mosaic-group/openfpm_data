#ifndef CSVWRITER_UNIT_TESTS_HPP_
#define CSVWRITER_UNIT_TESTS_HPP_

#include "CSVWriter.hpp"
#include "Vector/vector_test_util.hpp"

BOOST_AUTO_TEST_SUITE( csv_writer_test )


BOOST_AUTO_TEST_CASE( csv_writer_particles )
{
	// Allocate a property vector
	auto v_prp = allocate_openfpm(16);
	// Vector of position
	openfpm::vector<Point<3,float>> v_pos;

	// create a positional vector
	for (size_t i = 0 ; i < v_prp.size() ; i++)
	{
		Point<3,float> p({1.0,2.0,3.0});

		v_pos.add(p);
	}

	// CSVWriter test
	CSVWriter<openfpm::vector<Point<3,float>>, openfpm::vector<Point_test<float>> > csv_writer;

	// Write the CSV
	csv_writer.write("csv_out.csv",v_pos,v_prp);

	bool test = compare("csv_out.csv","csv_out_test.csv");
	BOOST_REQUIRE_EQUAL(true,test);


}

BOOST_AUTO_TEST_SUITE_END()

#endif
