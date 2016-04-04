/*
 * grid_performance_tests.hpp
 *
 *  Created on: Nov 1, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_PERFORMANCE_TESTS_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_PERFORMANCE_TESTS_HPP_

#include "grid_util_test.hpp"

openfpm::vector<std::string> testsg;
openfpm::vector<float> per_timesg;

BOOST_AUTO_TEST_CASE(grid_performance_set_obj)
{
	size_t sz[] = {128,128,128};

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);

	Point_test<float> f;
	f.fill();

	std::vector<double> times(N_STAT + 1);
	times[0] = 1000;

	for (size_t j = 0 ; j < 8 ; j++)
	{
		for (size_t i = 1 ; i < N_STAT+1 ; i++)
		{
			timer t;
			t.start();

			auto it = c3.getIterator();

			while (it.isNext())
			{
				c3.set(it.get(),f);

				++it;
			}

			t.stop();

			times[i] = t.getwct();
		}
		std::sort(times.begin(),times.end());
		sleep(5);
	}

	testsg.add("Grid so");
	per_timesg.add(times[0]);

}

BOOST_AUTO_TEST_CASE(grid_performance_set_other_grid)
{
	size_t sz[] = {128,128,128};

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);
	grid_cpu<3, Point_test<float> > c1(sz);
	c1.setMemory();

	std::vector<double> times(N_STAT + 1);
	times[0] = 1000;

	for (size_t j = 0 ; j < 8 ; j++)
	{
		for (size_t i = 1 ; i < N_STAT+1 ; i++)
		{
			timer t;
			t.start();

			auto it = c3.getIterator();

			while (it.isNext())
			{
				c3.set(it.get(),c1,it.get());

				++it;
			}

			t.stop();

			times[i] = t.getwct();
		}
		std::sort(times.begin(),times.end());
		sleep(5);
	}

	testsg.add("Grid sog");
	per_timesg.add(times[0]);

}

BOOST_AUTO_TEST_CASE(grid_performance_set_other_grid_encap)
{
	size_t sz[] = {128,128,128};

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);
	grid_cpu<3, Point_test<float> > c1(sz);
	c1.setMemory();

	std::vector<double> times(N_STAT + 1);
	times[0] = 1000;

	for (size_t j = 0 ; j < 8 ; j++)
	{
		for (size_t i = 1 ; i < N_STAT+1 ; i++)
		{
			timer t;
			t.start();

			auto it = c3.getIterator();

			while (it.isNext())
			{
				c3.set(it.get(),c1.get_o(it.get()));

				++it;
			}

			t.stop();

			times[i] = t.getwct();
		}
		std::sort(times.begin(),times.end());
		sleep(5);
	}

	testsg.add("Grid soge");
	per_timesg.add(times[0]);
}

BOOST_AUTO_TEST_CASE(grid_performance_duplicate)
{
	size_t sz[] = {128,128,128};

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);
	grid_cpu<3, Point_test<float> > c1;

	std::vector<double> times(N_STAT_SMALL + 1);
	times[0] = 1000;

	for (size_t j = 0 ; j < 8 ; j++)
	{
		for (size_t i = 1 ; i < N_STAT_SMALL+1 ; i++)
		{
			timer t;
			t.start();

			c1 = c3.duplicate();

			t.stop();

			times[i] = t.getwct();
		}
		std::sort(times.begin(),times.end());
		sleep(5);
	}

	testsg.add("Grid dup");
	per_timesg.add(times[0]);
}

/////// THIS IS NOT A TEST IT WRITE THE PERFORMANCE RESULT ///////

BOOST_AUTO_TEST_CASE(grid_performance_write_report)
{
	openfpm::vector<std::string> yn;
	openfpm::vector<openfpm::vector<float>> y;

	// Get the directory of the performance test files
	std::string per_dir(test_dir);

	// Reference time
	openfpm::vector<openfpm::vector<float>> y_ref;
	y_ref.load(per_dir + std::string("/ref_timesg"));

	load_and_combine(per_dir + std::string("/previous_measureg"),y,per_timesg);

	// Adding the dataset names
	if (y.size() != 0)
	{
		for (size_t j = 0; j < y.get(0).size(); j++)
			yn.add("config " + std::to_string(j));
	}

	// Google charts options
	GCoptions options;

	options.title = std::string("Grid Performances");
	options.yAxis = std::string("Time (seconds)");
	options.xAxis = std::string("Benchmark");
	options.stype = std::string("bars");

	std::stringstream g_test_desc;
	g_test_desc << "<h2>Grid performance test</h2>\n";
	g_test_desc << "<strong>128x128x128 Grid containing a Point_test<float></strong><br>";
	g_test_desc << "<strong>Grid so:</strong> Initialize each element of the grid<br>";
	g_test_desc << "<strong>Grid sog:</strong> Manual copy of two grids<br>";
	g_test_desc << "<strong>Grid soge:</strong> Manual copy of two grids in a different way<br>";
	g_test_desc << "<strong>Grid dup:</strong> Duplication of the grid (Duplication include grid creation time)<br>";


	cg.addHTML(g_test_desc.str());
	cg.AddHistGraph(testsg,y,yn,options);

	// compare the reference times with the actual times

	// calculate speed-up
	openfpm::vector<openfpm::vector<float>> y_ref_sup;

	speedup_calculate(y_ref_sup,y,y_ref,yn);

	std::stringstream g_test_spdesc;
	g_test_spdesc << "<h2>Grid speedup</h2>\n";
	g_test_spdesc << "The previous tests are compared with the best performances ever registered, ";
	g_test_spdesc << "the banded area indicate the upper and lower bounds of the best registrered performances.<br>";
	g_test_spdesc << "The lines are the latest 5 tests<br>";
	g_test_spdesc << "<strong>Line inside the area</strong>: The tested configuration has no improvement or degradation in performance<br>";
	g_test_spdesc << "<strong>Line break the upper bound</strong>: The tested configuration has improvement in performance<br>";
	g_test_spdesc << "<strong>Line break the lower bound</strong>: The tested configuration has degradation in performance<br>";
	g_test_spdesc << "<strong>Y axis:</strong> Performance change in percentage from the average of the best registered performances<br>";


	cg.addHTML(g_test_spdesc.str());
	cg.AddLinesGraph(testsg,y_ref_sup,yn,options);
}

#endif /* OPENFPM_DATA_SRC_GRID_GRID_PERFORMANCE_TESTS_HPP_ */
