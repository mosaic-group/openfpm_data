/*
 * graph_unit_test.hpp
 *
 *  Created on: Nov 22, 2014
 *      Author: i-bird
 */

#ifndef GRAPH_UNIT_TEST_HPP_
#define GRAPH_UNIT_TEST_HPP_

#include "map_graph.hpp"
#include "Point_test.hpp"

#define GS_SIZE 128

BOOST_AUTO_TEST_SUITE( graph_test )

BOOST_AUTO_TEST_CASE( graph_use)
{
	std::cout << "Graph unit test start" << "\n";

	 typedef Point_test<float> V;
	 typedef Point_test<float> E;

	 // Point
	 Point_test<float> p;
	 p.setx(1.0);
	 p.sety(2.0);
	 p.setz(3.0);
	 p.sets(4.0);

	 p.get<V::v>()[0] = 1.0;
	 p.get<V::v>()[1] = 2.0;
	 p.get<V::v>()[2] = 7.0;

	 p.get<V::t>()[0][0] = 10.0;
	 p.get<V::t>()[0][1] = 13.0;
	 p.get<V::t>()[0][2] = 8.0;
	 p.get<V::t>()[1][0] = 19.0;
	 p.get<V::t>()[1][1] = 23.0;
	 p.get<V::t>()[1][2] = 5.0;
	 p.get<V::t>()[2][0] = 4.0;
	 p.get<V::t>()[2][1] = 3.0;
	 p.get<V::t>()[2][2] = 11.0;

	//! define a test graph where the vertex and edge store several scalar vectorial and tensorial quantities

	Graph_CSR<V,E> g;

	//! Create a point

	//! try to construct a 2D dimensional cartesian graph


	// first create the vertex

	for (size_t i = 0 ; i < GS_SIZE ; i++)
	{
		for (size_t j = 0 ; j < GS_SIZE ; j++)
		{
			// Add vertex

			g.addVertex(p);
		}
	}

	// than create the edge

	// Create a grid, // NOTE this does not produce any memory grid, it retain only the information
	// and give a set of usefull function


	std::vector<size_t> gs;
	gs.push_back(GS_SIZE);
	gs.push_back(GS_SIZE);

	grid<2,void> g2(gs);

	// Create the edge 4 for each vertex

	for (size_t i = 0 ; i < GS_SIZE ; i++)
	{
		for (size_t j = 0 ; j < GS_SIZE ; j++)
		{
			// Add 4 edge

			g.addEdge<CheckExistence>(g2.LinId(i,j),g2.LinId(i+1,j),p);
			g.addEdge<CheckExistence>(g2.LinId(i,j),g2.LinId(i,j+1),p);
			g.addEdge<CheckExistence>(g2.LinId(i,j),g2.LinId(i-1,j),p);
			g.addEdge<CheckExistence>(g2.LinId(i,j),g2.LinId(i,j-1),p);
		}
	}

	std::cout << "Graph unit test end" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* GRAPH_UNIT_TEST_HPP_ */
