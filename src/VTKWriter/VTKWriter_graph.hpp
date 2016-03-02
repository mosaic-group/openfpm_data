/*
 * VTKWriter_graph.hpp
 *
 *  Created on: May 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VTKWRITER_GRAPH_HPP_
#define VTKWRITER_GRAPH_HPP_

/*! Property data store for scalar and vector
 *
 */
template<bool is_array>
struct vtk_vertex_node_array_scalar_selector
{
	/*! /brief Print the geometric informations in case it is not an array
	 *
	 * \tparam T Type of the vertex
	 * \tparam ele_v Attribute element to check
	 * \tparam G Graph of reference
	 * \tparam s_type Vertex spatial type information
	 *
	 * \param vo Vertex object container
	 * \param x Array to store geometric informations
	 * \param z_set Value set to true id z axis is in use
	 */
	template<typename T, typename ele_v, typename G, typename s_type>
	static inline void move(typename G::V_container &vo, s_type (&x)[3], bool &z_set)
	{
		if (G::V_type::attributes::name[T::value] == "x")
		{
			x[0] = convert<typename boost::remove_reference<decltype(vo.template get<T::value>())>::type>::template to<s_type>(vo.template get<T::value>());
		}
		else if (G::V_type::attributes::name[T::value] == "y")
		{
			x[1] = convert<typename boost::remove_reference<decltype(vo.template get<T::value>())>::type>::template to<s_type>(vo.template get<T::value>());
		}
		else if (G::V_type::attributes::name[T::value] == "z")
		{
			x[2] = convert<typename boost::remove_reference<decltype(vo.template get<T::value>())>::type>::template to<s_type>(vo.template get<T::value>());
			z_set = true;
		}
	}
};

/*! Template specialization in the case of array type attribute
 *
 */
template<>
struct vtk_vertex_node_array_scalar_selector<true>
{
	/*! \brief Store the geometric informations in case it is an array
	 *
	 * \tparam T Type of the vertex
	 * \tparam ele_v Attribute element to check
	 * \tparam G Graph of reference
	 * \tparam s_type Vertex spatial type information
	 *
	 * \param vo Vertex object container
	 * \param x Array to store geometric informations
	 * \param z_set Value set to true id z axis is in use
	 */
	template<typename T, typename ele_v, typename G, typename s_type>
	static inline void move(typename G::V_container &vo, s_type (&x)[3], bool &z_set)
	{
		if (std::extent<ele_v>::value == 3)
			z_set = true;

		for (size_t i = 0; i < std::extent<ele_v>::value; i++)
			x[i] = convert<typename boost::remove_reference<decltype(vo.template get<T::value>()[i])>::type>::template to<s_type>(vo.template get<T::value>()[i]);

	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to create a string containing all the vertex
 * properties
 *
 * \tparam G graph type
 * \tparam attr has the vertex attributes
 *
 */

template<typename G, bool attr>
struct vtk_vertex_node
{
	// Vertex spatial type information
	typedef typename G::V_type::s_type s_type;

	bool z_set;

	s_type (&x)[3];

	// Vertex object container
	typename G::V_container & vo;

	// vertex node string
	std::string & v_node;

	/*! \brief Constructor
	 *
	 * Create a vertex properties list
	 *
	 * \param v_node std::string that is filled with the graph properties in the GraphML format
	 * \param n_obj object container to access its properties for example encapc<...>
	 *
	 */
	vtk_vertex_node(std::string & v_node, typename G::V_container & n_obj, s_type (&x)[3])
	:z_set(false),x(x),vo(n_obj),v_node(v_node)
	{
	}

	//! \brief Write collected information
	void write()
	{
		v_node += std::to_string(x[0]) + " " + std::to_string(x[1]) + " " + std::to_string(x[2]) + "\n";
	}

	//! It call the functor for each member
	template<typename T>
	void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename G::V_type::type, boost::mpl::int_<T::value>>::type ele_v;

		// if the attribute name is x y or z, create a string with the value of the properties and store it

		vtk_vertex_node_array_scalar_selector<std::is_array<ele_v>::value>::template move<T, ele_v, G, s_type>(vo, x, z_set);

	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to create a string containing all the vertex
 * properties
 *
 * Specialization when we do not have vertex attributes
 *
 * \tparam G graph type
 *
 */

template<typename G>
struct vtk_vertex_node<G, false>
{
	// Vertex object container
	typename G::V_container & vo;

	// vertex node string
	std::string & v_node;

	/*! \brief Constructor
	 *
	 * Create a vertex properties list
	 *
	 * \param v_node std::string that is filled with the graph properties in the GraphML format
	 * \param n_obj object container to access its properties for example encapc<...>
	 *
	 */
	vtk_vertex_node(std::string & v_node, typename G::V_container & n_obj) :
			vo(n_obj), v_node(v_node)
	{
	}
	;

	//! It call the functor for each member
	template<typename T>
	void operator()(T& t)
	{
		v_node += "0 0 0\n";
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to create a string containing all the edge
 * properties
 *
 */

template<typename G>
struct vtk_edge_node
{
	// Vertex object container
	typename G::E_container & vo;

	// edge node string
	std::string & e_node;

	/*! \brief Constructor
	 *
	 * Create an edge node
	 *
	 * \param e_node std::string that is filled with the graph properties in the GraphML format
	 * \param n_obj object container to access the object properties for example encapc<...>
	 * \param n_prop number of properties
	 *
	 */
	vtk_edge_node(std::string & e_node, typename G::E_container & n_obj) :
			vo(n_obj), e_node(e_node)
	{
	}
	;

	/*! \brief Create a new node
	 *
	 * \param vc node number
	 *
	 */
	void new_node(size_t v_c, size_t s, size_t d)
	{
		// start a new node
		e_node += "2 " + std::to_string(s) + " " + std::to_string(d) + "\n";
	}
};

/*! \brief Property writer for scalar and vector
 *
 */
template<bool is_array>
struct prop_output_array_scalar_selector_vertex
{
	/*! \brief Writer in case the property is not an array
	 *
	 * \tparam ele_v Property element
	 * \tparam Graph Graph of reference
	 * \tparam i Property id
	 *
	 * \param v_out Buffer to write into
	 * \param g Graph
	 * \param p Property id
	 */
	template<typename ele_v, typename Graph, unsigned int i>
	static inline void write(std::string &v_out, const Graph &g, size_t p)
	{
		v_out += std::to_string(g.vertex(p).template get<i>()) + "\n";
	}
};

/*! \brief Property writer for vector
 *
 */
template<>
struct prop_output_array_scalar_selector_vertex<true>
{
	/*! \brief Writer in case the property is an array
	 *
	 * \tparam ele_v Property element
	 * \tparam Graph Graph of reference
	 * \tparam i Property id
	 *
	 * \param v_out Buffer to write into
	 * \param g Graph
	 * \param p Property id
	 */
	template<typename ele_v, typename Graph, unsigned int i>
	static inline void write(std::string &v_out, const Graph &g, size_t p)
	{
		for (size_t j = 0; j < 2; j++)
		{
			v_out += std::to_string(g.vertex(p).template get<i>()[j]) + " ";
		}

		if (std::extent<ele_v>::value == 2)
			v_out += "0";
		else
			v_out += std::to_string(g.vertex(p).template get<i>()[2]);

		v_out += "\n";
	}
};

/*! \brief Property writer for scalar and vector
 *
 */
template<bool is_array>
struct prop_output_array_scalar_selector_edge
{
	/*! \brief Writer in case the property is not an array
	 *
	 * \tparam ele_v Property element
	 * \tparam Graph Graph of reference
	 * \tparam i Property id
	 *
	 * \param v_out Buffer to write into
	 * \param g Graph
	 * \param p Property id
	 */
	template<typename ele_v, typename Graph, unsigned int i>
	static inline void write(std::string &v_out, const Graph &g, const typename Graph::E_container &edge)
	{
		v_out += std::to_string(edge.template get<i>()) + "\n";
	}
};

/*! \brief Property writer for vector
 *
 */
template<>
struct prop_output_array_scalar_selector_edge<true>
{
	/*! \brief Writer in case the property is an array
	 *
	 * \tparam ele_v Property element
	 * \tparam Graph Graph of reference
	 * \tparam i Property id
	 *
	 * \param v_out Buffer to write into
	 * \param g Graph
	 * \param p Property id
	 */
	template<typename ele_v, typename Graph, unsigned int i>
	static inline void write(std::string &v_out, const Graph &g, const typename Graph::E_container &edge)
	{
		for (size_t j = 0; j < 2; j++)
		{
			v_out += std::to_string(edge.template get<i>()[j]) + " ";
		}

		if (std::extent<ele_v>::value == 2)
			v_out += "0";
		else
			v_out += std::to_string(edge.template get<i>()[2]);

		v_out += "\n";
	}
};

/*! \brief Property writer for scalar and vector, it fill the vertex data (needed for edge representation in vtk)
 *
 */
template<bool is_array>
struct prop_output_array_scalar_selector_edge_fill_vertex
{
	/*! \brief Writer in case the property is not an array
	 *
	 * \param v_out Buffer to write into
	 */
	static inline void write(std::string &v_out)
	{
		v_out += "0\n";
	}
};

/*! \brief Property writer for vector
 *
 */
template<>
struct prop_output_array_scalar_selector_edge_fill_vertex<true>
{
	/*! \brief Writer in case the property is an array
	 *
	 * \param v_out Buffer to write into
	 */
	static inline void write(std::string &v_out)
	{
		v_out += "0 0 0\n";
	}
};

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
 * \tparam i the property we are going to write
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

	static std::string get_point_data(const Graph & g)
	{
		//! vertex node output string
		std::string v_out;

		//! Get a vertex iterator
		auto it = g.getVertexIterator();

		// if there is the next element
		while (it.isNext())
		{
			typedef typename boost::mpl::at<typename Graph::V_type::type, boost::mpl::int_<i>>::type ele_v;
			prop_output_array_scalar_selector_vertex<std::is_array<ele_v>::value>::template write<ele_v, Graph, i>(v_out, g, it.get());

			// increment the iterator and counter
			++it;
		}

		return v_out;
	}

	/*! \brief For each edge set the value, set 1 on vertices, needed by vtk file format
	 *
	 * \tparam i edge property to print
	 *
	 */

	static std::string get_cell_data(const Graph & g)
	{
		//! vertex node output string
		std::string e_out;

		//! Get a vertex iterator
		auto it_v = g.getVertexIterator();

		// if there is the next element
		while (it_v.isNext())
		{
			// Print the property
			typedef typename boost::mpl::at<typename Graph::E_type::type, boost::mpl::int_<i>>::type ele_v;
			prop_output_array_scalar_selector_edge_fill_vertex<std::is_array<ele_v>::value>::write(e_out);

			// increment the iterator and counter
			++it_v;
		}

		//! Get an edge iterator
		auto it_e = g.getEdgeIterator();

		// if there is the next element
		while (it_e.isNext())
		{
			typedef typename boost::mpl::at<typename Graph::E_type::type, boost::mpl::int_<i>>::type ele_v;
			prop_output_array_scalar_selector_edge<std::is_array<ele_v>::value>::template write<ele_v, Graph, i>(e_out, g, g.edge(it_e.get()));

			// increment the iterator and counter
			++it_e;
		}

		return e_out;
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

		// Type of the property
		std::string type;

		typedef typename boost::mpl::at<typename Graph::V_type::type, boost::mpl::int_<i>>::type T;

		// Check if T is a supported format
		// for now we support only scalar of native type
		if (std::rank<T>::value == 1)
		{
			//Get type of the property
			type = getType<typename std::remove_all_extents<T>::type>();

			// if the type is not supported return
			if (type.size() == 0)
				return v_out;

			// Create point data properties
			v_out += "VECTORS " + get_attributes_vertex() + " " + type + "\n";
		}
		else
		{
			type = getType<T>();

			// if the type is not supported return
			if (type.size() == 0)
				return v_out;

			// Create point data properties
			v_out += "SCALARS " + get_attributes_vertex() + " " + type + "\n";

			// Default lookup table
			v_out += "LOOKUP_TABLE default\n";

		}

		// return the vertex list
		return v_out;
	}

	/*! \brief Given a Graph return the cell data header for a typename T
	 *
	 * \tparam T type to write
	 * \param n_node number of the node
	 *
	 */

	static std::string get_cell_property_header(size_t prop)
	{
		//! edge node output string
		std::string e_out;

		// Type of the property
		std::string type;

		typedef typename boost::mpl::at<typename Graph::E_type::type, boost::mpl::int_<i>>::type T;

		// Check if T is a supported format
		// for now we support only scalar of native type
		if (std::is_array<T>::value == true && std::is_array<typename std::remove_extent<T>::type>::value == false)
		{
			//Get type of the property
			type = getType<typename std::remove_all_extents<T>::type>();

			// if the type is not supported return
			if (type.size() == 0)
				return e_out;

			// Create point data properties
			e_out += "VECTORS " + get_attributes_edge() + " " + type + "\n";
		}
		else
		{
			type = getType<T>();

			// if the type is not supported return
			if (type.size() == 0)
				return e_out;

			// Create point data properties
			e_out += "SCALARS " + get_attributes_edge() + " " + type + "\n";

			// Default lookup table
			e_out += "LOOKUP_TABLE default\n";

		}

		// return the vertex list
		return e_out;
	}

	/*! \brief Get the attributes name for vertex
	 *
	 */

	static std::string get_attributes_vertex()
	{
		return Graph::V_type::attributes::name[i];
	}

	/*! \brief Get the attributes name for edge
	 *
	 */

	static std::string get_attributes_edge()
	{
		return Graph::E_type::attributes::name[i];
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
class prop_output<false, Graph, i>
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

	/*! \brief For each edge set the value
	 *
	 * \tparam i edge property to print
	 *
	 */

	static std::string get_cell_data(const Graph & g)
	{
		//! vertex node output string
		std::string e_out;

		//! Get a vertex iterator
		auto it_v = g.getVertexIterator();

		// if there is the next element
		while (it_v.isNext())
		{
			// Print the property
			e_out += std::to_string(0) + "\n";

			// increment the iterator and counter
			++it_v;
		}

		//! Get an edge iterator
		auto it_e = g.getEdgeIterator();

		// if there is the next element
		while (it_e.isNext())
		{
			// Print the property
			e_out += std::to_string(g.edge(it_e.get()).template get<i>()) + "\n";

			// increment the iterator and counter
			++it_e;
		}

		return e_out;
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

		std::string type = getType<boost::fusion::result_of::at<typename Graph::V_type::type, boost::mpl::int_<i>>>("attr" + std::to_string(prop));

		// if the type is not supported return
		if (type.size() == 0)
		{
			return v_out;
		}

		// Create point data properties
		v_out += "SCALARS " + get_attributes_vertex() + " " + type + "\n";

		// Default lookup table
		v_out += "LOOKUP_TABLE default\n";

		// return the vertex list
		return v_out;
	}

	/*! \brief Given a Graph return the cell data header for a typename T
	 *
	 * \param n_node number of the node
	 *
	 */

	static std::string get_cell_property_header(size_t prop)
	{
		//! edge node output string
		std::string e_out;

		// Type of the property
		std::string type;

		typedef typename boost::mpl::at<typename Graph::E_type::type, boost::mpl::int_<i>>::type T;

		// Check if T is a supported format
		// for now we support only scalar of native type
		if (std::is_array<T>::value == true && std::is_array<typename std::remove_extent<T>::type>::value == false)
		{
			//Get type of the property
			type = getType<typename std::remove_all_extents<T>::type>();

			// if the type is not supported return
			if (type.size() == 0)
				return e_out;

			// Create point data properties
			e_out += "VECTORS " + get_attributes_edge() + " " + type + "\n";
		}
		else
		{
			type = getType<T>();

			// if the type is not supported return
			if (type.size() == 0)
				return e_out;

			// Create point data properties
			e_out += "SCALARS " + get_attributes_edge() + " " + type + "\n";

			// Default lookup table
			e_out += "LOOKUP_TABLE default\n";

		}

		// return the vertex list
		return e_out;
	}

	/*! \brief Get the attributes name for vertex
	 *
	 */

	static std::string get_attributes_vertex()
	{
		return Graph::V_type::attributes::name[i];
	}

	/*! \brief Get the attributes name for edge
	 *
	 */

	static std::string get_attributes_edge()
	{
		return Graph::E_type::attributes::name[i];
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
struct prop_out_vertex
{
	// property output string
	std::string & v_out;

	// Graph that we are processing
	const Graph & g;

	/*! \brief constructor
	 *
	 * \param v_out string to fill with the vertex properties
	 *
	 */
	prop_out_vertex(std::string & v_out, const Graph & g) :
			v_out(v_out), g(g)
	{
	}
	;

	//! It produce an output for each property
	template<typename T>
	void operator()(T& t) const
	{
		// actual string size
		size_t sz = v_out.size();

		// Produce the point properties header
		v_out += prop_output<has_attributes<typename Graph::V_type>::value, Graph, T::value>::get_point_property_header(t);

		// If the output has changed, we have to write the properties
		if (v_out.size() != sz)
		{
			std::string attr = prop_output<has_attributes<typename Graph::V_type>::value, Graph, T::value>::get_attributes_vertex();

			// Produce point data
			v_out += prop_output<has_attributes<typename Graph::V_type>::value, Graph, T::value>::get_point_data(g);
		}
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
struct prop_out_edge
{
	// property output string
	std::string & e_out;

	// Graph that we are processing
	const Graph & g;

	/*! \brief constructor
	 *
	 * \param v_out string to fill with the vertex properties
	 *
	 */
	prop_out_edge(std::string & e_out, const Graph & g) :
			e_out(e_out), g(g)
	{
	}
	;

	//! It produce an output for each property
	template<typename T>
	void operator()(T& t) const
	{
		// actual string size
		size_t sz = e_out.size();

		// Produce the point properties header
		e_out += prop_output<has_attributes<typename Graph::E_type>::value, Graph, T::value>::get_cell_property_header(t);

		// If the output has changed, we have to write the properties
		if (e_out.size() != sz)
		{
			std::string attr = prop_output<has_attributes<typename Graph::E_type>::value, Graph, T::value>::get_attributes_edge();

			// Produce cell data
			e_out += prop_output<has_attributes<typename Graph::E_type>::value, Graph, T::value>::get_cell_data(g);
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

template<typename Graph>
class VTKWriter<Graph, VTK_GRAPH>
{
	const Graph & g;

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
		e_out += "LINES " + std::to_string(g.getNEdge()) + " " + std::to_string(3 * g.getNEdge()) + "\n";

		// return the vertex properties string
		return e_out;
	}

	/*! \brief Create the VTK point definition
	 *
	 * \tparam s_type spatial type of the data
	 * \tparam attr false x,y,z are set to 0 for each vertex
	 *
	 */

	template<bool attr> std::string get_point_list()
	{
		//! VTK spatial information
		typename Graph::V_type::s_type x[3] = { 0, 0, 0 };

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
			vtk_vertex_node<Graph, attr> vn(v_out, obj, x);

			// Iterate through all the vertex and create the vertex list
			boost::mpl::for_each<boost::mpl::range_c<int, 0, Graph::V_type::max_prop > >(vn);

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
		for (size_t i = 0; i < g.getNVertex(); i++)
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

	/*! \brief Get the point data header
	 *
	 * \return a string with the point data header for VTK format
	 *
	 */

	std::string get_cell_data_header()
	{
		std::string v_out;

		v_out += "CELL_DATA " + std::to_string(g.getNVertex() + g.getNEdge()) + "\n";

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
	VTKWriter(const Graph & g) :
			g(g)
	{
	}

	/*! \brief It write a VTK file from a graph
	 *
	 * \tparam prp_out which properties to output [default = -1 (all)]
	 *
	 * \param file path where to write
	 * \param name of the graph
	 * \param file_type specify if it is a VTK BINARY or ASCII file [default = ASCII]
	 *
	 */

	template<int prp = -1> bool write(std::string file, std::string graph_name = "Graph", file_type ft = file_type::ASCII)
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
		// Cell data header
		std::string cell_data_header;
		// Cell data
		std::string cell_data;

		// VTK header
		vtk_header = "# vtk DataFile Version 3.0\n" + graph_name + "\n";

		// Choose if binary or ASCII
		if (ft == file_type::ASCII)
		{
			vtk_header += "ASCII\n";
		}
		else
		{
			vtk_header += "BINARY\n";
		}

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

		// Get the cell data header
		cell_data_header = get_cell_data_header();

		// For each property in the vertex type produce a point data

		prop_out_vertex<Graph> pp(point_data, g);

		if (prp == -1)
			boost::mpl::for_each<boost::mpl::range_c<int, 0, Graph::V_type::max_prop> >(pp);
		else
			boost::mpl::for_each<boost::mpl::range_c<int, prp, prp> >(pp);

		// For each property in the edge type produce a point data

		prop_out_edge<Graph> ep(cell_data, g);

		if (prp == -1)
			boost::mpl::for_each<boost::mpl::range_c<int, 0, Graph::E_type::max_prop> >(ep);
		else
			boost::mpl::for_each<boost::mpl::range_c<int, prp, prp> >(ep);

		// write the file
		std::ofstream ofs(file);

		// Check if the file is open
		if (ofs.is_open() == false)
		{
			std::cerr << "Error cannot create the VTK file: " + file;
		}

		ofs << vtk_header << point_prop_header << point_list << vertex_prop_header << vertex_list << edge_prop_header << edge_list << point_data_header << point_data << cell_data_header << cell_data;

		// Close the file

		ofs.close();

		// Completed succefully
		return true;
	}
};

#endif /* VTKWRITER_GRAPH_HPP_ */
