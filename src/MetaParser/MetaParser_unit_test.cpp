/*
 * MetaParser_unit_test.cpp
 *
 *  Created on: Mar 3, 2019
 *      Author: i-bird
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "MetaParser.hpp"

BOOST_AUTO_TEST_SUITE( vtk_writer_tests )

BOOST_AUTO_TEST_CASE( vtk_writer_meta_parser_use )
{
	std::string test_options("time = 5.0");
	std::string test_options2("time=6.0");
	std::string test_options3("time =7.0");
	std::string test_options4("time= 8.0");
	std::string test_options5("time=     9.0");

	double time = 1.0;

	MetaParser_options opts;
	opts.add_options()
    	    ("time", MetaParser_def::value<double>(&time));

	MetaParser mp(opts);
	mp.parse(test_options);

	BOOST_REQUIRE_EQUAL(time,5.0);

	time = 0.0;
	bool exist = mp.getOption("time",time);

	BOOST_REQUIRE_EQUAL(exist,true);
	BOOST_REQUIRE_EQUAL(time,5.0);

	exist = mp.getOption("invalid",time);
	BOOST_REQUIRE_EQUAL(exist,false);

	// parse another

	mp.clear();

	time = 0.0;
	mp.parse(test_options2);

	BOOST_REQUIRE_EQUAL(time,6.0);

	time = 0.0;
	exist = mp.getOption("time",time);

	BOOST_REQUIRE_EQUAL(exist,true);
	BOOST_REQUIRE_EQUAL(time,6.0);

	exist = mp.getOption("invalid",time);
	BOOST_REQUIRE_EQUAL(exist,false);

	// parse another

	mp.clear();

	time = 0.0;
	mp.parse(test_options3);

	BOOST_REQUIRE_EQUAL(time,7.0);

	time = 0.0;
	exist = mp.getOption("time",time);

	BOOST_REQUIRE_EQUAL(exist,true);
	BOOST_REQUIRE_EQUAL(time,7.0);

	exist = mp.getOption("invalid",time);
	BOOST_REQUIRE_EQUAL(exist,false);

	// parse another

	mp.clear();

	time = 0.0;
	mp.parse(test_options4);

	BOOST_REQUIRE_EQUAL(time,8.0);

	time = 0.0;
	exist = mp.getOption("time",time);

	BOOST_REQUIRE_EQUAL(exist,true);
	BOOST_REQUIRE_EQUAL(time,8.0);

	exist = mp.getOption("invalid",time);
	BOOST_REQUIRE_EQUAL(exist,false);

	// parse another

	mp.clear();

	time = 0.0;
	mp.parse(test_options5);

	BOOST_REQUIRE_EQUAL(time,9.0);

	time = 0.0;
	exist = mp.getOption("time",time);

	BOOST_REQUIRE_EQUAL(exist,true);
	BOOST_REQUIRE_EQUAL(time,9.0);

	exist = mp.getOption("invalid",time);
	BOOST_REQUIRE_EQUAL(exist,false);

}

BOOST_AUTO_TEST_SUITE_END()
