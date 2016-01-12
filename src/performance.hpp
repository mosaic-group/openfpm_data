/*
 * performance.hpp
 *
 *  Created on: Jan 11, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_PERFORMANCE_HPP_
#define OPENFPM_DATA_SRC_PERFORMANCE_HPP_

#include "Plot/GoogleChart.hpp"
#include "timer.hpp"

#define N_STAT 256
#define N_STAT_SMALL 32
#define N_TRY 8

#ifdef PERFORMANCE_TEST

BOOST_AUTO_TEST_SUITE( performance )

//// Include tests ////////

#include "Grid/grid_performance_tests.hpp"
#include "Vector/vector_performance_test.hpp"

#define MEASURE_SET 5

/*! \brief Load the previous test result combine with the actual result and save
 *
 * \param file file that contain the previous result
 * \param yn vector with the name of the dataset
 * \param y vector with the data to load
 * \param nc number of colums
 *
 */
void load_and_combine(std::string file, openfpm::vector<std::string> & yn, openfpm::vector<openfpm::vector<float>> & y, openfpm::vector<float> & per_times, size_t nc)
{
	// Load the previous measure and combine the previous measure with the actual measure
	y.clear();
	yn.clear();

	y.load(file);
	y.resize(nc);

	for(size_t i = 0 ; i < y.size() ; i++)
	{
		if (y.get(i).size() >= MEASURE_SET)
			y.get(i).remove(0);
		y.get(i).add(per_times.get(i));
	}

	for (size_t j = 0; j < y.get(0).size(); j++)
	{
		if (j < y.get(0).size() - 1)
			yn.add("previous " + std::to_string(j));
		else
			yn.add("actual");
	}

	y.save(file);
}

BOOST_AUTO_TEST_CASE(performance_report_out)
{
	GoogleChart cg;

	openfpm::vector<std::string> yn;
	openfpm::vector<openfpm::vector<float>> y;

	load_and_combine("/home/i-bird/Desktop/MOSAIC/OpenFPM_project/openfpm_pdata/openfpm_data/previous_measureg",yn,y,per_timesg,testsg.size());

	// Google charts options
	GCoptions options;

	options.title = std::string("Grid Performances");
	options.yAxis = std::string("Time (seconds)");
	options.xAxis = std::string("Benchmark");
	options.stype = std::string("bars");


	cg.addHTML("<h2>Grid performance test</h2>");
	cg.AddColumsGraph(testsg,y,yn,options);


	load_and_combine("/home/i-bird/Desktop/MOSAIC/OpenFPM_project/openfpm_pdata/openfpm_data/previous_measurev",yn,y,per_timesv,testsv.size());

	options.title = std::string("Vector Performances");
	options.yAxis = std::string("Time (seconds)");
	options.xAxis = std::string("Benchmark");
	options.stype = std::string("bars");

	cg.addHTML("<h2>Vector performance test</h2>");
	cg.AddColumsGraph(testsv,y,yn,options);
	cg.write("gc_out.html");
}

BOOST_AUTO_TEST_SUITE_END()

#endif

#endif /* OPENFPM_DATA_SRC_PERFORMANCE_HPP_ */
