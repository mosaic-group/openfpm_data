/*
 * GraphMLWriter_unit_tests.hpp
 *
 *  Created on: Dec 9, 2014
 *      Author: i-bird
 */

#ifndef GRAPHMLWRITER_UNIT_TESTS_HPP_
#define GRAPHMLWRITER_UNIT_TESTS_HPP_

#define GS_SIZE 128

#include "GraphMLWriter.hpp"
#include "Graph/CartesianGraphFactory.hpp"

BOOST_AUTO_TEST_SUITE( graphml_writer_test )

/*!
 *
 * Test node and edge
 *
 */

struct ne_cp
{
	//! The node contain several properties
	typedef boost::fusion::vector<float,float,float,double,long int,int,std::string> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! define attributes names
	struct attributes
	{
		static const std::string name[];
	};

	//! The data
	type data;

	//! x property id in boost::fusion::vector
	static const unsigned int x = 0;
	//! y property id in boost::fusion::vector
	static const unsigned int y = 1;
	//! z property id in boost::fusion::vector
	static const unsigned int z = 2;
	//! float_num property id in boost::fusion::vector
	static const unsigned int float_num = 3;
	//! double_num property id in boost::fusion::vector
	static const unsigned int double_num = 4;
	//! long_num property id in boost::fusion::vector
	static const unsigned int long_num = 5;
	//! integer property id in boost::fusion::vector
	static const unsigned int integer = 6;
	//! string property id in boost::fusion::vector
	static const unsigned int string = 7;
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 8;
};

// Initialize the attributes strings array
const std::string ne_cp::attributes::name[] = {"x","y","z","float_num","double_num","long_num","integer","string"};

BOOST_AUTO_TEST_CASE( graphml_writer_use)
{
	//! Create a graph

	CartesianGraphFactory<3,Graph_CSR<ne_cp,ne_cp>> g_factory;

	// Cartesian grid
	std::vector<size_t> sz;
	sz.push_back(GS_SIZE);
	sz.push_back(GS_SIZE);
	sz.push_back(GS_SIZE);

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	Graph_CSR<ne_cp,ne_cp> g_csr = g_factory.construct<5,float,2>(sz,box);

	GraphMLWriter<Graph_CSR<ne_cp,ne_cp>> gw(g_csr);

	gw.write("test_graph.gml");
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* GRAPHMLWRITER_UNIT_TESTS_HPP_ */
