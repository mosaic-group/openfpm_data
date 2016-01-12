/*
 * Plot_unit_tests.hpp
 *
 *  Created on: Jan 9, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_PLOT_PLOT_UNIT_TESTS_HPP_
#define OPENFPM_DATA_SRC_PLOT_PLOT_UNIT_TESTS_HPP_

#include "GoogleChart.hpp"

BOOST_AUTO_TEST_SUITE( plot_unit_test )

BOOST_AUTO_TEST_CASE( google_chart )
{
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
	cg.AddColumsGraph(x,y,yn,options);
	cg.write("gc_out.html");

	bool test = compare("gc_out.html","gc_out_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart2 )
{
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
	cg.AddColumsGraph(x,y,yn,options);
	cg.write("gc_out2.html");

	bool test = compare("gc_out2.html","gc_out2_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart3 )
{
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
	cg.AddColumsGraph(x,y,yn,options);
	cg.write("gc_out3.html");

	bool test = compare("gc_out3.html","gc_out3_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart4 )
{
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
	cg.AddColumsGraph(x,y,yn);
	cg.write("gc_out4.html");

	bool test = compare("gc_out4.html","gc_out4_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart5 )
{
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
	cg.AddColumsGraph(x,y);
	cg.write("gc_out5.html");

	bool test = compare("gc_out5.html","gc_out5_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart6 )
{
	openfpm::vector<openfpm::vector<float>> y;

	// Each colums can have multiple data-set
	y.add({2.2,1.3,4.5,0.6});
	y.add({5.0,6.1,1.3,2.6});
	y.add({2.1,1.0,6.1,9.3});
	y.add({1.1,6.1,3.0,2.0});
	y.add({3.3,0.3,0.0,6.2});
	y.add({2.0,1.1,4.0,6.1});

	GoogleChart cg;
	cg.AddColumsGraph(y);
	cg.write("gc_out6.html");

	bool test = compare("gc_out6.html","gc_out6_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( google_chart_with_inject_HTML )
{
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
	cg.AddColumsGraph(x,y,yn,options);
	cg.addHTML("<h2>Before second graph</h2>");
	cg.AddColumsGraph(x,y,yn,options);
	cg.addHTML("<h2>Before third graph</h2>");
	cg.AddColumsGraph(x,y,yn,options);
	cg.addHTML("<h2>At the end</h2>");
	cg.write("gc_out7.html");

	bool test = compare("gc_out7.html","gc_out7_test.html");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_PLOT_PLOT_UNIT_TESTS_HPP_ */
