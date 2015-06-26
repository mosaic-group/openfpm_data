/*
 * VTKWriter_graph.hpp
 *
 *  Created on: May 5, 2015
 *      Author: i-bird
 */

#ifndef VTKWRITER_GRAPH_HPP_
#define VTKWRITER_GRAPH_HPP_

/*! \brief This class specialize functions in the case the type T
 * has or not defined attributes
 *
 * In C++ partial specialization of a function is not allowed so we have to
 * encapsulate this function in a class
 *
 * \tparam has_attributes parameter that specialize the function in case the vertex
 *         define or not attributes name
 *
 * \tparam Graph type of graph we are processing
 * \tparam p the property we are going to write
 *
 */

template<bool has_attributes, typename Graph, unsigned int i>
class prop_output
{
public:

	/*! \brief For each vertex set the value
	 *
	 * \tparam i vertex property to print
	 *
	 */

	static std::string get_point_data(Graph & g)
	{
		//! vertex node output string
		std::string v_out;

		//! Get a vertex iterator
		auto it = g.getVertexIterator();

		// if there is the next element
		while (it.isNext())
		{
			// Print the property
			v_out += std::to_string(g.vertex(it.get()).template get<i>()) + "\n";

			// increment the iterator and counter
			++it;
		}

		return v_out;
	}

	/*! \brief Given a Graph return the point data header for a typename T
	 *
	 * \tparam T type to write
	 * \param n_node number of the node
	 *
	 */

	static std::string get_point_property_header(size_t prop)
	{
		//! vertex node output string
		std::string v_out;

		// Check if T is a supported format
		// for now we support only scalar of native type

		std::string type = getType<typename boost::fusion::result_of::at<typename Graph::V_type::type,boost::mpl::int_<i>>::type>();

		// if the type is not supported return
		if (type.size() == 0)
		{return v_out;}

		// Create point data properties
		v_out += "SCALARS " + get_attributes() + " " + type + "\n";

		// Default lookup table
		v_out += "LOOKUP_TABLE default\n";

		// return the vertex list
		return v_out;
	}

	/*! \brief Get the attributes name
	 *
	 */

	static std::string get_attributes()
	{
		return Graph::V_type::attributes::name[i];
	}
};

/*! \brief This class specialize functions in the case the type T
 * has not defined attributes
 *
 * In C++ partial specialization of a function is not allowed so we have to
 * encapsulate this function in a class
 *
 * \tparam has_attributes parameter that specialize the function in case the vertex
 *         define or not attributes name
 *
 * \tparam i id of the property we are going to write
 *
 */

template<typename Graph, unsigned int i>
class prop_output<false,Graph,i>
{
	/*! \brief For each vertex set the value
	 *
	 * \tparam i vertex property to print
	 *
	 */

	static std::string get_point_data(Graph & g)
	{
		//! vertex node output string
		std::string v_out;

		//! Get a vertex iterator
		auto it = g.getVertexIterator();

		// if there is the next element
		while (it.isNext())
		{
			// Print the property
			v_out += std::to_string(g.vertex(it.get()).template get<i>()) + "\n";

			// increment the iterator and counter
			++it;
		}

		return v_out;
	}

	/*! \brief Given a Graph return the point data header for a typename T
	 *
	 * \tparam T type to write
	 *
	 * \param n_node number of the node
	 * \param prop id of the property
	 *
	 */

	static std::string get_point_property_header(size_t prop)
	{
		//! vertex node output string
		std::string v_out;

		// Check if T is a supported format
		// for now we support only scalar of native type

		std::string type = getType<boost::fusion::result_of::at<typename Graph::V_type::type,boost::mpl::int_<i>>>("attr" + std::to_string(prop));

		// if the type is not supported return
		if (type.size() == 0)
		{return v_out;}

		// Create point data properties
		v_out += "SCALARS " + get_attributes() + " " + type + "\n";

		// Default lookup table
		v_out += "LOOKUP_TABLE default\n";

		// return the vertex list
		return v_out;
	}

	/*! \brief Get the attributes name
	 *
	 * \tparam has_prop true if T has properties name defined
	 * \tparam T type to process
	 *
	 * \param i attribute to get
	 *
	 */

	static std::string get_attributes()
	{
		return "attr" + std::to_string(i);
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to produce at output for each property
 *
 * \tparam Graph graph we are processing
 *
 * \param dim Dimensionality
 * \param S type of grid
 *
 */

template<typename Graph>
struct prop_out
{
	// property output string
	std::string & v_out;

	// Graph that we are processing
	Graph & g;

	/*! \brief constructor
	 *
	 * \param v_out string to fill with the vertex properties
	 *
	 */
	prop_out(std::string & v_out, Graph & g)
	:v_out(v_out),g(g)
	{};

	//! It produce an output for each property
    template<typename T>
    void operator()(T& t) const
    {
    	// actual string size
    	size_t sz = v_out.size();

		// Produce the point properties header
		v_out += prop_output<has_attributes<typename Graph::V_type>::value ,Graph,T::value>::get_point_property_header(t);

		// If the output has changed, we have to write the properties
		if (v_out.size() != sz)
		{
			std::string attr = prop_output<has_attributes<typename Graph::V_type>::value,Graph,T::value>::get_attributes();

			// Produce point data
			v_out += prop_output<has_attributes<typename Graph::V_type>::value ,Graph ,T::value>::get_point_data(g);
		}
    }
};

/*!
 *
 * It write a VTK format file in case for a graph
 *
 * \tparam Type of graph
 *
 */

template <typename Graph>
class VTKWriter<Graph,GRAPH>
{
	Graph & g;

	/*! \brief It get the vertex properties list
	 *
	 * It get the vertex properties list of the vertex defined as VTK header
	 *
	 * \return a string that define the vertex properties in graphML format
	 *
	 */

	std::string get_vertex_properties_list()
	{
		//! vertex property output string
		std::string v_out;

		// write the number of vertex
		v_out += "VERTICES " + std::to_string(g.getNVertex()) + " " + std::to_string(g.getNVertex() * 2) + "\n";

		// return the vertex properties string
		return v_out;
	}

	/*! \brief It get the vertex properties list
	 *
	 * It get the vertex properties list of the vertex defined as a VTK header
	 *
	 * \return a string that define the vertex properties in graphML format
	 *
	 */

	std::string get_point_properties_list()
	{
		//! vertex property output string
		std::string v_out;

		// write the number of vertex
		v_out += "POINTS " + std::to_string(g.getNVertex()) + " float" + "\n";

		// return the vertex properties string
		return v_out;
	}

	/*! \brief It get the edge properties list
	 *
	 * It get the edge properties list of the edge defined as a GraphML header
	 *
	 * \return a string that define the edge properties in graphML format
	 *
	 */

	std::string get_edge_properties_list()
	{
		//! vertex property output string
		std::string e_out;

		// write the number of lines
		e_out += "LINES " + std::to_string(g.getNEdge()) + " " + std::to_string(3*g.getNEdge()) + "\n";

		// return the vertex properties string
		return e_out;
	}

	/*! \brief Create the VTK point definition
	 *
	 * \tparam s_type spatial type of the data
	 * \tparam attr false x,y,z are set to 0 for each vertex
	 *
	 */

	template <bool attr> std::string get_point_list()
	{
		//! VTK spatial information
		typename Graph::V_type::s_type x[3] = {0,0,0};

		//! vertex node output string
		std::string v_out;

		//! Get a vertex iterator
		auto it = g.getVertexIterator();

		// if there is the next element
		while (it.isNext())
		{
			// Get vtk vertex node
			auto obj = g.vertex(it.get());

			// create a vertex list functor
			vtk_vertex_node<Graph,attr> vn(v_out,obj,x);

			// Iterate through all the vertex and create the vertex list
			boost::mpl::for_each< boost::mpl::range_c<int,0,Graph::V_type::max_prop-1> >(vn);

			// write the node string
			vn.write();

			// increment the iterator and counter
			++it;
		}

		// return the vertex list
		return v_out;
	}

	/*! \brief Create the VTK vertex definition
	 *
	 * \tparam s_type spatial type of the data
	 * \tparam attr false x,y,z are set to 0 for each vertex
	 *
	 */

	std::string get_vertex_list()
	{
		//! vertex node output string
		std::string v_out;

		//! For each point create a vertex
		for (size_t i = 0 ; i < g.getNVertex() ; i++)
		{
			v_out += "1 " + std::to_string(i) + "\n";
		}

		// return the vertex list
		return v_out;
	}

	/*! \brief Get the point data header
	 *
	 * \return a string with the point data header for VTK format
	 *
	 */

	std::string get_point_data_header()
	{
		std::string v_out;

		v_out += "POINT_DATA " + std::to_string(g.getNVertex()) + "\n";

		return v_out;
	}

	/*! \brief Return the edge list
	 *
	 * \return the edge list
	 *
	 */

	std::string get_edge_list()
	{
		//! edge node output string
		std::string e_out;

		//! Get an edge iterator
		auto it = g.getEdgeIterator();

		// if there is the next element
		while (it.isNext())
		{
			// create an edge list functor
//			edge_node<Graph> en(e_out,g.edge(it.get()));

			e_out += "2 " + std::to_string(it.source()) + " " + std::to_string(it.target()) + "\n";

			// increment the operator
			++it;
		}

		// return the edge list
		return e_out;
	}

public:

	/*!
	 *
	 * VTKWriter constructor, it take a graph and write a GraphML format
	 *
	 * \param g Graph to write
	 *
	 */
	VTKWriter(Graph & g)
	:g(g)
	{}

	/*! \brief It write a VTK file from a graph
	 *
	 * \tparam prp_out which properties to output [default = -1 (all)]
	 *
	 * \param file path where to write
	 * \param name of the graph
	 * \param file_type specify if it is a VTK BINARY or ASCII file [default = ASCII]
	 *
	 */

	template<int prp = -1> bool write(std::string file, std::string graph_name="Graph", file_type ft = file_type::ASCII)
	{
		// Check that the Vertex type define x y and z attributes

		if (has_attributes<typename Graph::V_type>::value == false)
		{
			std::cerr << "Error writing a graph: Vertex must has defines x,y,z properties" << "\n";
			return false;
		}

		// Header for the vtk
		std::string vtk_header;
		// Point list of the VTK
		std::string point_list;
		// Vertex list of the VTK
		std::string vertex_list;
		// Graph header
		std::string vtk_binary_or_ascii;
		// Edge list of the GraphML
		std::string edge_list;
		// vertex properties header
		std::string point_prop_header;
		// edge properties header
		std::string vertex_prop_header;
		// edge properties header
		std::string edge_prop_header;
		// Data point header
		std::string point_data_header;
		// Data point
		std::string point_data;

		// VTK header
		vtk_header = "# vtk DataFile Version 3.0\n"
				     + graph_name + "\n";

		// Choose if binary or ASCII
		if (ft == file_type::ASCII)
		{vtk_header += "ASCII\n";}
		else
		{vtk_header += "BINARY\n";}

		// Data type for graph is DATASET POLYDATA
		vtk_header += "DATASET POLYDATA\n";

		// point properties header
		point_prop_header = get_point_properties_list();

		// Get point list
		point_list = get_point_list<has_attributes<typename Graph::V_type>::value>();

		// vertex properties header
		vertex_prop_header = get_vertex_properties_list();

		// Get vertex list
		vertex_list = get_vertex_list();

		// Edge properties header
		edge_prop_header = get_edge_properties_list();

		// Get the edge graph list
		edge_list = get_edge_list();

		// Get the point data header
		point_data_header = get_point_data_header();

		// For each property in the vertex type produce a point data

		prop_out<Graph> pp(point_data, g);

		if (prp == -1)
			boost::mpl::for_each< boost::mpl::range_c<int,0, Graph::V_type::max_prop> >(pp);
		else
			boost::mpl::for_each< boost::mpl::range_c<int,prp, prp> >(pp);

		// write the file
		std::ofstream ofs(file);

		// Check if the file is open
		if (ofs.is_open() == false)
		{std::cerr << "Error cannot create the VTK file: " + file;}

		ofs << vtk_header << point_prop_header << point_list <<
				vertex_prop_header << vertex_list << edge_prop_header << edge_list << point_data_header << point_data;

		// Close the file

		ofs.close();

		// Completed succefully
		return true;
	}
};


#endif /* VTKWRITER_GRAPH_HPP_ */
