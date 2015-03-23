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

	// Test if duplicate work

	Graph_CSR<V,E> g_dup = g.duplicate();

	// check if the two graph matches

	for (size_t i = 0 ; i < g.getNVertex() ; i++)
	{
		BOOST_REQUIRE_EQUAL(g.vertex(i).template get<V::x>(),g_dup.vertex(i).template get<V::x>());
		BOOST_REQUIRE_EQUAL(g.vertex(i).template get<V::y>(),g_dup.vertex(i).template get<V::y>());
		BOOST_REQUIRE_EQUAL(g.vertex(i).template get<V::z>(),g_dup.vertex(i).template get<V::z>());
		BOOST_REQUIRE_EQUAL(g.vertex(i).template get<V::s>(),g_dup.vertex(i).template get<V::s>());

		for (int j = 0 ; j < 3 ; j++)
		{
			BOOST_REQUIRE_EQUAL(g.vertex(i).template get<V::v>()[j],g_dup.vertex(i).template get<V::v>()[j]);
		}

		for (int j = 0 ; j < 3 ; j++)
		{
			for (int k = 0 ; k < 3 ; k++)
			{
				BOOST_REQUIRE_EQUAL(g.vertex(i).template get<V::t>()[j][k],g_dup.vertex(i).template get<V::t>()[j][k]);
			}
		}
	}

	// Test the no edge case

	Graph_CSR<V> g_no_edge;

	// Create a tree

	// root
	g_no_edge.addVertex(p);

	// 2 leaf

	g_no_edge.addVertex(p);
	g_no_edge.addEdge(0,1);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(0,2);

	// 4 leaf

	g_no_edge.addVertex(p);
	g_no_edge.addEdge(1,3);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(1,4);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(2,5);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(2,6);

	// 8 leaf

	g_no_edge.addVertex(p);
	g_no_edge.addEdge(3,7);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(3,8);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(4,9);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(4,10);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(5,11);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(5,12);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(6,13);
	g_no_edge.addVertex(p);
	g_no_edge.addEdge(6,14);

	// Check that each vertex has 2 child

	size_t start = 0;

	for (size_t i = 0 ; i < g_no_edge.getNVertex() ; i++)
	{
		if (g_no_edge.getNChilds(i) == 0)
			continue;

		size_t s1 = g_no_edge.getChild(i,0);
		size_t s2 = g_no_edge.getChild(i,1);

		BOOST_REQUIRE_EQUAL(s1+1,s2);

		size_t a = log2(i + 1);
		size_t start = (a == 0)?1:(2 << (a-1));
		start -= 1;
		size_t start_p = 2 << a;
		start_p -= 1;

		BOOST_REQUIRE_EQUAL(s1,start_p + (i - start)*2);
	}

	std::cout << "Graph unit test end" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* GRAPH_UNIT_TEST_HPP_ */
