#ifndef CSVWRITER_UNIT_TESTS_HPP_
#define CSVWRITER_UNIT_TESTS_HPP_

#include "CSVWriter.hpp"
#include "Vector/vector_test_util.hpp"

BOOST_AUTO_TEST_SUITE( csv_writer_test )


BOOST_AUTO_TEST_CASE( csv_writer_particles )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	{
	// Allocate a property vector
	auto v_prp = allocate_openfpm_prp(16);
	// Vector of position
	openfpm::vector<Point<3,float>> v_pos;

	// create a positional vector
	for (size_t i = 0 ; i < v_prp.size() ; i++)
	{
		Point<3,float> p({1.0,2.0,3.0});

		v_pos.add(p);
	}

	// CSVWriter test
	CSVWriter<openfpm::vector<Point<3,float>>, openfpm::vector<Point_test_prp<float>> > csv_writer;

	// Write the CSV
	csv_writer.write("csv_out.csv",v_pos,v_prp);

	bool test = compare("csv_out.csv","test_data/csv_out_test.csv");
	BOOST_REQUIRE_EQUAL(true,test);
	}

	{
	// Allocate a property vector
	auto v_prp = allocate_openfpm_aggregate_with_complex(16);
	// Vector of position
	openfpm::vector<Point<3,float>> v_pos;

	// create a positional vector
	for (size_t i = 0 ; i < v_prp.size() ; i++)
	{
		Point<3,float> p({1.0,2.0,3.0});

		v_pos.add(p);
	}

	// CSVWriter test
	CSVWriter<openfpm::vector<Point<3,float>>, openfpm::vector< aggregate<float,float,float,float,float[3],float[3][3],openfpm::vector<int>> > > csv_writer;

	// Write the CSV
	csv_writer.write("csv_out_unk.csv",v_pos,v_prp);

	// In case of SE_CLASS3 enabled the number of properties change

#ifndef SE_CLASS3
	bool test = compare("csv_out_unk.csv","test_data/csv_out_unk_test.csv");
	BOOST_REQUIRE_EQUAL(true,test);
#endif
	}

}

BOOST_AUTO_TEST_SUITE_END()

#endif
