/*
 * Plot_unit_tests.hpp
 *
 *  Created on: Jan 9, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_PLOT_PLOT_UNIT_TESTS_HPP_
#define OPENFPM_DATA_SRC_PLOT_PLOT_UNIT_TESTS_HPP_

#include "GoogleChart.hpp"
#include "Plot/util.hpp"

BOOST_AUTO_TEST_SUITE( plot_unit_test )


BOOST_AUTO_TEST_CASE( google_chart_bar_string )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	//! [Producing an Histogram graph]

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<size_t>> y;
	openfpm::vector<std::string> yn;

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	// Each colum can have multiple data set (in this case 4 dataset)
	// Each dataset can have a name
	yn.add("dataset1");
	yn.add("dataset2");
	yn.add("dataset3");
	yn.add("dataset4");

	// Each colums can have multiple data-set
	y.add({2,3,5,6});
	y.add({5,6,1,6});
	y.add({2,1,6,9});
	y.add({1,6,3,2});
	y.add({3,3,0,6});
	y.add({2,1,4,6});

	// Google charts options
	GCoptions options;

	options.title = std::string("Example");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.stype = std::string("bars");
	options.barWD = true;

	// it say that the colum4 must me represented with a line
	options.stypeext = std::string("{3: {type: 'line'}}");

	GoogleChart cg;
	cg.AddHistGraph(x,y,yn,options);
	cg.write("gc_out_sc.html");

	//! [Producing an Histogram graph]

	bool test = compare("gc_out_sc.html","test_data/gc_out_sc_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}


BOOST_AUTO_TEST_CASE( google_chart )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	//! [Producing an Histogram graph]

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<size_t>> y;
	openfpm::vector<std::string> yn;

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	// Each colum can have multiple data set (in this case 4 dataset)
	// Each dataset can have a name
	yn.add("dataset1");
	yn.add("dataset2");
	yn.add("dataset3");
	yn.add("dataset4");

	// Each colums can have multiple data-set
	y.add({2,3,5,6});
	y.add({5,6,1,6});
	y.add({2,1,6,9});
	y.add({1,6,3,2});
	y.add({3,3,0,6});
	y.add({2,1,4,6});

	// Google charts options
	GCoptions options;

	options.title = std::string("Example");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.stype = std::string("bars");

	// it say that the colum4 must me represented with a line
	options.stypeext = std::string("{3: {type: 'line'}}");

	GoogleChart cg;
	cg.AddHistGraph(x,y,yn,options);
	cg.write("gc_out.html");

	//! [Producing an Histogram graph]

	bool test = compare("gc_out.html","gc_out_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart2 )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<float>> y;
	openfpm::vector<std::string> yn;

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	// Each colum can have multiple data set (in this case 4 dataset)
	// Each dataset can have a name
	yn.add("dataset1");
	yn.add("dataset2");
	yn.add("dataset3");
	yn.add("dataset4");

	// Each colums can have multiple data-set
	y.add({2.2,1.3,4.5,0.6});
	y.add({5.0,6.1,1.3,2.6});
	y.add({2.1,1.0,6.1,9.3});
	y.add({1.1,6.1,3.0,2.0});
	y.add({3.3,0.3,0.0,6.2});
	y.add({2.0,1.1,4.0,6.1});

	// Google charts options
	GCoptions options;

	options.title = std::string("Example");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.stype = std::string("bars");

	GoogleChart cg;
	cg.AddHistGraph(x,y,yn,options);
	cg.write("gc_out2.html");

	bool test = compare("gc_out2.html","gc_out2_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart3 )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<float>> y;
	openfpm::vector<std::string> yn;

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	// Each colum can have multiple data set (in this case 4 dataset)
	// Each dataset can have a name
	yn.add("dataset1");
	yn.add("dataset2");
	yn.add("dataset3");
	yn.add("dataset4");

	// Each colums can have multiple data-set
	y.add({2.2,1.3,4.5,0.6});
	y.add({5.0,6.1,1.3,2.6});
	y.add({2.1,1.0,6.1,9.3});
	y.add({1.1,6.1,3.0,2.0});
	y.add({3.3,0.3,0.0,6.2});
	y.add({2.0,1.1,4.0,6.1});

	// Google charts options
	GCoptions options;

	options.title = std::string("Example");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");

	GoogleChart cg;
	cg.AddHistGraph(x,y,yn,options);
	cg.write("gc_out3.html");

	bool test = compare("gc_out3.html","gc_out3_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart4 )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<float>> y;
	openfpm::vector<std::string> yn;

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	// Each colum can have multiple data set (in this case 4 dataset)
	// Each dataset can have a name
	yn.add("dataset1");
	yn.add("dataset2");
	yn.add("dataset3");
	yn.add("dataset4");

	// Each colums can have multiple data-set
	y.add({2.2,1.3,4.5,0.6});
	y.add({5.0,6.1,1.3,2.6});
	y.add({2.1,1.0,6.1,9.3});
	y.add({1.1,6.1,3.0,2.0});
	y.add({3.3,0.3,0.0,6.2});
	y.add({2.0,1.1,4.0,6.1});

	GoogleChart cg;
	cg.AddHistGraph(x,y,yn);
	cg.write("gc_out4.html");

	bool test = compare("gc_out4.html","gc_out4_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart5 )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<float>> y;

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	// Each colums can have multiple data-set
	y.add({2.2,1.3,4.5,0.6});
	y.add({5.0,6.1,1.3,2.6});
	y.add({2.1,1.0,6.1,9.3});
	y.add({1.1,6.1,3.0,2.0});
	y.add({3.3,0.3,0.0,6.2});
	y.add({2.0,1.1,4.0,6.1});

	GoogleChart cg;
	cg.AddHistGraph(x,y);
	cg.write("gc_out5.html");

	bool test = compare("gc_out5.html","gc_out5_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart6 )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	openfpm::vector<openfpm::vector<float>> y;

	// Each colums can have multiple data-set
	y.add({2.2,1.3,4.5,0.6});
	y.add({5.0,6.1,1.3,2.6});
	y.add({2.1,1.0,6.1,9.3});
	y.add({1.1,6.1,3.0,2.0});
	y.add({3.3,0.3,0.0,6.2});
	y.add({2.0,1.1,4.0,6.1});

	GoogleChart cg;
	cg.AddHistGraph(y);
	cg.write("gc_out6.html");

	bool test = compare("gc_out6.html","gc_out6_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart_with_inject_HTML )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	//! [Producing a set of histograms graphs]

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<size_t>> y;
	openfpm::vector<std::string> yn;

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	// Each colum can have multiple data set (in this case 4 dataset)
	// Each dataset can have a name
	yn.add("dataset1");
	yn.add("dataset2");
	yn.add("dataset3");
	yn.add("dataset4");

	// Each colums can have multiple data-set
	y.add({2,3,5,6});
	y.add({5,6,1,6});
	y.add({2,1,6,9});
	y.add({1,6,3,2});
	y.add({3,3,0,6});
	y.add({2,1,4,6});

	// Google charts options
	GCoptions options;

	options.title = std::string("Example");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.stype = std::string("bars");

	// it say that the colum4 must me represented with a line
	options.stypeext = std::string("{3: {type: 'line'}}");

	GoogleChart cg;
	//
	cg.addHTML("<h2>Before first graph</h2>");
	cg.AddHistGraph(x,y,yn,options);
	cg.addHTML("<h2>Before second graph</h2>");
	cg.AddHistGraph(x,y,yn,options);
	cg.addHTML("<h2>Before third graph</h2>");
	cg.AddHistGraph(x,y,yn,options);
	cg.addHTML("<h2>At the end</h2>");
	cg.write("gc_out7.html");

	//! [Producing a set of histograms graphs]

	bool test = compare("gc_out7.html","gc_out7_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart_linear_plot )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	//! [Producing lines graph with style]

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<double>> y;
	openfpm::vector<std::string> yn;

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	// Here we specify how many lines we have
	// first Line
	yn.add("line1");

	// second line + 2 intervals (Error bands)
	yn.add("line2");
	yn.add("interval");
	yn.add("interval");
	yn.add("interval");
	yn.add("interval");

	// third line + 1 interval (Error bands)
	yn.add("line3");
	yn.add("interval");
	yn.add("interval");

	// Each line can have multiple intervals or error bars
	// The first number specify the bottom line
	// The last three numbers specify the top line + error band (min, max)
	// The middle 5 line specify the middle lines + one external error band + one internal error band

	y.add({0.10,0.20,0.19,0.22,0.195,0.215,0.35,0.34,0.36});
	y.add({0.11,0.21,0.18,0.22,0.19,0.215,0.36,0.35,0.37});
	y.add({0.12,0.22,0.21,0.23,0.215,0.225,0.35,0.34,0.36});
	y.add({0.15,0.25,0.20,0.26,0.22,0.255,0.36,0.35,0.37});
	y.add({0.09,0.29,0.25,0.30,0.26,0.295,0.35,0.34,0.36});
	y.add({0.08,0.28,0.27,0.29,0.275,0.285,0.36,0.35,0.37});

	// Google charts options
	GCoptions options;

	options.title = std::string("Example");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.lineWidth = 1.0;
	options.intervalext = std::string("{'i2': { 'color': '#4374E0', 'style':'bars', 'lineWidth':4, 'fillOpacity':1 } }");

	GoogleChart cg;
	cg.AddLinesGraph(x,y,yn,options);
	cg.write("gc_plot_out.html");

	//! [Producing lines graph with style]

	bool test = compare("gc_plot_out.html","gc_plot_out_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart_linear_plot2 )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	//! [Producing lines]

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<double>> y;

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	// Each line can have multiple intervals or error bars
	// The first number specify the bottom line
	// The last three numbers specify the top line + error band (min, max)
	// The middle 5 line specify the middle lines + one external error band + one internal error band

	y.add({0.10,0.20,0.19,0.22,0.195,0.215,0.35,0.34,0.36});
	y.add({0.11,0.21,0.18,0.22,0.19,0.215,0.36,0.35,0.37});
	y.add({0.12,0.22,0.21,0.23,0.215,0.225,0.35,0.34,0.36});
	y.add({0.15,0.25,0.20,0.26,0.22,0.255,0.36,0.35,0.37});
	y.add({0.09,0.29,0.25,0.30,0.26,0.295,0.35,0.34,0.36});
	y.add({0.08,0.28,0.27,0.29,0.275,0.285,0.36,0.35,0.37});

	// Google charts options
	GCoptions options;

	options.title = std::string("Example");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.lineWidth = 1.0;

	GoogleChart cg;
	cg.AddLinesGraph(x,y,options);
	cg.write("gc_plot2_out.html");

	//! [Producing lines]

	bool test = compare("gc_plot2_out.html","gc_plot2_out_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

//! [Definition of a function]

double f(double x)
{
	return x*x;
}

//! [Definition of a function]

BOOST_AUTO_TEST_CASE( plot_util )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() != 0)
		return;

	//! [fill a vector]

	openfpm::vector<double> x;

	Fill1D(0.0,2.0,5,x);

	BOOST_REQUIRE_EQUAL(x.get(0),0.0);
	BOOST_REQUIRE_EQUAL(x.get(1),0.5);
	BOOST_REQUIRE_EQUAL(x.get(2),1.0);
	BOOST_REQUIRE_EQUAL(x.get(3),1.5);
	BOOST_REQUIRE_EQUAL(x.get(4),2.0);

	//! [fill a vector]

	x.clear();

	//! [fill a vector with a function]

	Fill1D(0.0,2.0,5,x,f);

	BOOST_REQUIRE_EQUAL(x.get(0),0.0);
	BOOST_REQUIRE_EQUAL(x.get(1),0.25);
	BOOST_REQUIRE_EQUAL(x.get(2),1.0);
	BOOST_REQUIRE_EQUAL(x.get(3),2.25);
	BOOST_REQUIRE_EQUAL(x.get(4),4.0);

	//! [fill a vector function]
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_PLOT_PLOT_UNIT_TESTS_HPP_ */
