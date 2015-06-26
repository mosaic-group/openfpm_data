/*
 * VTKWriter_unit_tests.hpp
 *
 *  Created on: May 6, 2015
 *      Author: Pietro Incardona
 */

#ifndef VTKWRITER_UNIT_TESTS_HPP_
#define VTKWRITER_UNIT_TESTS_HPP_

BOOST_AUTO_TEST_SUITE( vtk_writer_test )

/* \brief Sub-domain vertex graph node
 *
 */

struct vertex
{
	//! The node contain 3 unsigned long integer for communication computation memory and id
	typedef boost::fusion::vector<float,float,float,float,size_t,double,unsigned char,long int> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! type of the positional field
	typedef float s_type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data
	type data;

	//! computation property id in boost::fusion::vector
	static const unsigned int x = 0;
	//! computation property id in boost::fusion::vector
	static const unsigned int y = 1;
	//! memory property id in boost::fusion::vector
	static const unsigned int z = 2;
	//! computation property id in boost::fusion::vector
	static const unsigned int prp1 = 3;
	//! computation property id in boost::fusion::vector
	static const unsigned int prp2 = 4;
	//! memory property id in boost::fusion::vector
	static const unsigned int prp3 = 5;
	//! memory property id in boost::fusion::vector
	static const unsigned int prp4 = 6;
	//! memory property sub_id in boost::fusion::vector
	static const unsigned int prp5 = 7;

	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 8;

	/*!
	 * Default constructor
	 *
	 */
	vertex()
	{

	}

	/*! \brief Initialize the VTKVertex
	 *
	 * \param
	 *
	 */
	vertex(float x, float y, float z)
	{
		boost::fusion::at_c<vertex::x>(data) = x;
		boost::fusion::at_c<vertex::y>(data) = y;
		boost::fusion::at_c<vertex::z>(data) = z;
	}
};

// use the vertex like the edge
typedef vertex edge;

const std::string vertex::attributes::name[] = {"x","y","z","prp1","prp2","prp3","prp4","prp5"};

BOOST_AUTO_TEST_CASE( vtk_writer_use_graph)
{
	// Create some graphs and output them

	std::cout << "Graph unit test start" << "\n";

	// Graph

	Graph_CSR<vertex,edge> gr;

	// Create a cube graph

	gr.addVertex(vertex(0.0,0.0,0.0));
	gr.addVertex(vertex(0.0,0.0,1.0));
	gr.addVertex(vertex(0.0,1.0,0.0));
	gr.addVertex(vertex(0.0,1.0,1.0));
	gr.addVertex(vertex(1.0,0.0,0.0));
	gr.addVertex(vertex(1.0,0.0,1.0));
	gr.addVertex(vertex(1.0,1.0,0.0));
	gr.addVertex(vertex(1.0,1.0,1.0));

	gr.addEdge(0,6);
	gr.addEdge(6,4);
	gr.addEdge(4,0);

	gr.addEdge(0,2);
	gr.addEdge(2,6);
	gr.addEdge(6,0);

	gr.addEdge(0,3);
	gr.addEdge(3,2);
	gr.addEdge(2,0);

	gr.addEdge(0,1);
	gr.addEdge(1,3);
	gr.addEdge(3,0);

	gr.addEdge(2,7);
	gr.addEdge(7,6);
	gr.addEdge(6,2);

	gr.addEdge(2,3);
	gr.addEdge(3,7);
	gr.addEdge(7,2);

	gr.addEdge(4,6);
	gr.addEdge(6,7);
	gr.addEdge(7,4);

	gr.addEdge(4,7);
	gr.addEdge(7,5);
	gr.addEdge(5,4);

	gr.addEdge(0,4);
	gr.addEdge(4,5);
	gr.addEdge(5,0);

	gr.addEdge(0,5);
	gr.addEdge(5,1);
	gr.addEdge(1,0);

	gr.addEdge(1,5);
	gr.addEdge(5,7);
	gr.addEdge(7,1);

	gr.addEdge(1,7);
	gr.addEdge(7,3);
	gr.addEdge(3,1);

	// Write the VTK file

	VTKWriter<Graph_CSR<vertex,edge>,GRAPH> vtk(gr);
	vtk.write("vtk_graph.vtk");

	// check that match

	bool test = compare("vtk_graph.vtk","vtk_graph_test.vtk");
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_CASE( vtk_writer_use_vector_box)
{
	// Create a vector of boxes
	openfpm::vector<Box<2,float>> vb;

	vb.add(Box<2,float>({0.2,0.2},{1.0,0.5}));
	vb.add(Box<2,float>({0.0,0.0},{0.2,0.2}));
	vb.add(Box<2,float>({0.2,0.0},{0.5,0.2}));
	vb.add(Box<2,float>({0.5,0.0},{1.0,0.2}));
	vb.add(Box<2,float>({0.0,0.2},{0.2,0.5}));
	vb.add(Box<2,float>({0.0,0.5},{1.0,1.0}));

	// Create a writer and write
	VTKWriter<openfpm::vector<Box<2,float>>,VECTOR_BOX> vtk_box;
	vtk_box.add(vb);
	vtk_box.write("vtk_box.vtk");

	// Check that match
	bool test = compare("vtk_box.vtk","vtk_box_test.vtk");
	BOOST_REQUIRE_EQUAL(test,true);

	// Create a vector of boxes
	openfpm::vector<Box<3,float>> vb2;

	vb2.add(Box<3,float>({0.2,0.2,0.0},{1.0,0.5,0.5}));
	vb2.add(Box<3,float>({0.0,0.0,0.0},{0.2,0.2,0.5}));
	vb2.add(Box<3,float>({0.2,0.0,0.0},{0.5,0.2,0.5}));
	vb2.add(Box<3,float>({0.5,0.0,0.0},{1.0,0.2,0.5}));
	vb2.add(Box<3,float>({0.0,0.2,0.0},{0.2,0.5,0.5}));
	vb2.add(Box<3,float>({0.0,0.5,0.0},{1.0,1.0,0.5}));

	// Create a writer and write
	VTKWriter<openfpm::vector<Box<3,float>>,VECTOR_BOX> vtk_box2;
	vtk_box2.add(vb2);
	vtk_box2.write("vtk_box_3D.vtk");

	// Check that match
	test = compare("vtk_box_3D.vtk","vtk_box_3D_test.vtk");
	BOOST_REQUIRE_EQUAL(test,true);

	// Create a vector of boxes
	openfpm::vector<Box<3,float>> vb3;
	vb3.add(Box<3,float>({0.2,0.2,0.5},{1.0,0.5,1.0}));
	vb3.add(Box<3,float>({0.0,0.0,0.5},{0.2,0.2,1.0}));
	vb3.add(Box<3,float>({0.2,0.0,0.5},{0.5,0.2,1.0}));
	vb3.add(Box<3,float>({0.5,0.0,0.5},{1.0,0.2,1.0}));
	vb3.add(Box<3,float>({0.0,0.2,0.5},{0.2,0.5,1.0}));
	vb3.add(Box<3,float>({0.0,0.5,0.5},{1.0,1.0,1.0}));

	// Create a writer and write
	VTKWriter<openfpm::vector<Box<3,float>>,VECTOR_BOX> vtk_box3;
	vtk_box3.add(vb2);
	vtk_box3.add(vb3);
	vtk_box3.write("vtk_box_3D_2.vtk");

	// Check that match
	test = compare("vtk_box_3D_2.vtk","vtk_box_3D_2_test.vtk");
	BOOST_REQUIRE_EQUAL(test,true);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* VTKWRITER_UNIT_TESTS_HPP_ */
