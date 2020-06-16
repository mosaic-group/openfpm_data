/*
 * grid_performance_tests.hpp
 *
 *  Created on: Nov 1, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_PERFORMANCE_TESTS_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_PERFORMANCE_TESTS_HPP_

#include "Grid/grid_util_test.hpp"
#include "util/stat/common_statistics.hpp"

// Property tree
struct report_grid_copy_func_tests
{
	boost::property_tree::ptree graphs;
};

report_grid_copy_func_tests report_grid_funcs;

BOOST_AUTO_TEST_SUITE( grid_performance )

BOOST_AUTO_TEST_CASE(grid_performance_set_obj)
{
	size_t sz[] = {128,128,128};

	report_grid_funcs.graphs.put("performance.grid.set(0).grid.x",sz[0]);
	report_grid_funcs.graphs.put("performance.grid.set(0).grid.y",sz[1]);
	report_grid_funcs.graphs.put("performance.grid.set(0).grid.z",sz[2]);

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);

	Point_test<float> f __attribute__((aligned(16)));
	f.fill();

	std::vector<double> times(N_STAT + 1);

	for (size_t i = 0 ; i < N_STAT+1 ; i++)
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

	double mean;
	double dev;
	standard_deviation(times,mean,dev);

	report_grid_funcs.graphs.put("performance.grid.set(0).x.data.name","Grid_so");
	report_grid_funcs.graphs.put("performance.grid.set(0).y.data.mean",mean);
	report_grid_funcs.graphs.put("performance.grid.set(0).y.data.dev",dev);
}

BOOST_AUTO_TEST_CASE(grid_performance_set_other_grid)
{
	size_t sz[] = {128,128,128};

	report_grid_funcs.graphs.put("performance.grid.set(1).grid.x",sz[0]);
	report_grid_funcs.graphs.put("performance.grid.set(1).grid.y",sz[1]);
	report_grid_funcs.graphs.put("performance.grid.set(1).grid.z",sz[2]);

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);
	grid_cpu<3, Point_test<float> > c1(sz);
	c1.setMemory();

	std::vector<double> times(N_STAT + 1);

	for (size_t i = 0 ; i < N_STAT+1 ; i++)
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

	double mean;
	double dev;
	standard_deviation(times,mean,dev);

	report_grid_funcs.graphs.put("performance.grid.set(1).x.data.name","Grid_sog");
	report_grid_funcs.graphs.put("performance.grid.set(1).y.data.mean",mean);
	report_grid_funcs.graphs.put("performance.grid.set(1).y.data.dev",dev);
}

BOOST_AUTO_TEST_CASE(grid_performance_set_other_grid_encap)
{
	size_t sz[] = {128,128,128};

	report_grid_funcs.graphs.put("performance.grid.set(2).grid.x",sz[0]);
	report_grid_funcs.graphs.put("performance.grid.set(2).grid.y",sz[1]);
	report_grid_funcs.graphs.put("performance.grid.set(2).grid.z",sz[2]);

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);
	grid_cpu<3, Point_test<float> > c1(sz);
	c1.setMemory();

	std::vector<double> times(N_STAT + 1);

	for (size_t i = 0 ; i < N_STAT+1 ; i++)
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

	double mean;
	double dev;
	standard_deviation(times,mean,dev);

	report_grid_funcs.graphs.put("performance.grid.set(2).x.data.name","Grid_soge");
	report_grid_funcs.graphs.put("performance.grid.set(2).y.data.mean",mean);
	report_grid_funcs.graphs.put("performance.grid.set(2).y.data.dev",dev);
}

BOOST_AUTO_TEST_CASE(grid_performance_duplicate)
{
	size_t sz[] = {128,128,128};

	report_grid_funcs.graphs.put("performance.grid.set(3).grid.x",sz[0]);
	report_grid_funcs.graphs.put("performance.grid.set(3).grid.y",sz[1]);
	report_grid_funcs.graphs.put("performance.grid.set(3).grid.z",sz[2]);

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);
	grid_cpu<3, Point_test<float> > c1;

	std::vector<double> times(N_STAT_SMALL + 1);

	for (size_t i = 0 ; i < N_STAT_SMALL+1 ; i++)
	{
		timer t;
		t.start();

		c1 = c3.duplicate();

		t.stop();

		times[i] = t.getwct();
	}

	double mean;
	double dev;
	standard_deviation(times,mean,dev);

	report_grid_funcs.graphs.put("performance.grid.set(3).x.data.name","Grid_dup");
	report_grid_funcs.graphs.put("performance.grid.set(3).y.data.mean",mean);
	report_grid_funcs.graphs.put("performance.grid.set(3).y.data.dev",dev);
}

/////// THIS IS NOT A TEST IT WRITE THE PERFORMANCE RESULT ///////

BOOST_AUTO_TEST_CASE(grid_performance_write_report)
{
	// Create a graphs

	report_grid_funcs.graphs.put("graphs.graph(0).type","line");
	report_grid_funcs.graphs.add("graphs.graph(0).title","Grid set functions (so/sog/soge) and duplicate (dup) performance");
	report_grid_funcs.graphs.add("graphs.graph(0).x.title","Tests");
	report_grid_funcs.graphs.add("graphs.graph(0).y.title","Time seconds");
	report_grid_funcs.graphs.add("graphs.graph(0).y.data(0).source","performance.grid.set(#).y.data.mean");
	report_grid_funcs.graphs.add("graphs.graph(0).x.data(0).source","performance.grid.set(#).x.data.name");
	report_grid_funcs.graphs.add("graphs.graph(0).y.data(0).title","Actual");
	report_grid_funcs.graphs.add("graphs.graph(0).interpolation","lines");

	boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
	boost::property_tree::write_xml("grid_performance_funcs.xml", report_grid_funcs.graphs,std::locale(),settings);

	GoogleChart cg;

	std::string file_xml_ref(test_dir);
	file_xml_ref += std::string("/openfpm_data/grid_performance_funcs_ref.xml");

	StandardXMLPerformanceGraph("grid_performance_funcs.xml",file_xml_ref,cg);

	addUpdtateTime(cg,1);

	cg.write("grid_performance_funcs.html");
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_GRID_GRID_PERFORMANCE_TESTS_HPP_ */
