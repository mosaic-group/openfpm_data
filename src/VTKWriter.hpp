/*
 * VTKWriter.hpp
 *
 *  Created on: Dec 15, 2014
 *      Author: Pietro Incardona
 */

#ifndef VTKWRITER_HPP_
#define VTKWRITER_HPP_

#include "Graph/map_graph.hpp"
#include <iostream>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <fstream>
#include "util/common.hpp"

/*! \brief Get the type
 *
 * It convert T to a string identify the corrispondent type in VTK format
 *
 */

template <typename T> std::string getType()
{
	// Create a property string based on the type of the property
	if (typeid(T) == typeid(float))
		return "float";
	else if (typeid(T) == typeid(double))
		return "double";
	else if (typeid(T) == typeid(char))
		return "char";
	else if (typeid(T) == typeid(unsigned char))
		return "unsigned_char";
	else if (typeid(T) == typeid(short))
		return "short";
	else if (typeid(T) == typeid(unsigned short))
		return "unsigned_short";
	else if (typeid(T) == typeid(int))
		return "int";
	else if (typeid(T) == typeid(unsigned int))
		return "unsigned_int";
	else if (typeid(T) == typeid(long int))
		return "long";
	else if (typeid(T) == typeid(unsigned long int))
		return "unsigned_long";
	else if (typeid(T) == typeid(bool))
		return "bit";
	return "";
}

/*! \brief Set a conversion map between A and B
 *
 * Convert A to B
 *
 * \tparam B destination type
 * \tparam A source type
 *
 */

template<typename A>
class convert
{
public:
	template<typename B> static B to(const A & data)
	{
		return static_cast<B>(data);
	}
};

/*! \brief Partial specialization when A is a string
 *
 *
 */

template<>
class convert<std::string>
{
public:
	template<typename B> static B to(const std::string & data)
	{
		return atof(data.c_str());
	}
};

/*! \brief It specify the VTK output file type
 *
 */

enum file_type
{
	BINARY,
	ASCII
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

	s_type (& x)[3];

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
	:x(x),vo(n_obj),v_node(v_node)
	{
	};

	//! \brief Write collected information
	void write()
	{
		v_node += std::to_string(x[0]) + " " + std::to_string(x[1]) + " " + std::to_string(x[2]) + "\n";
	}

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
    	// if the attribute name is x y or z, create a string with the value of the properties and store it
		if (G::V_type::attributes::name[T::value] == "x"){x[0] = convert<typename boost::remove_reference<decltype(vo.template get<T::value>())>::type>::template to<s_type>(vo.template get<T::value>());}
		else if (G::V_type::attributes::name[T::value] == "y"){x[1] = convert<typename boost::remove_reference<decltype(vo.template get<T::value>())>::type>::template to<s_type>(vo.template get<T::value>());}
		else if (G::V_type::attributes::name[T::value] == "z"){x[2] = convert<typename boost::remove_reference<decltype(vo.template get<T::value>())>::type>::template to<s_type>(vo.template get<T::value>());}
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
struct vtk_vertex_node<G,false>
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
	vtk_vertex_node(std::string & v_node, typename G::V_container & n_obj)
	:vo(n_obj),v_node(v_node)
	{
	};

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
	vtk_edge_node(std::string & e_node, typename G::E_container & n_obj)
	:vo(n_obj),e_node(e_node)
	{
	};

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

#define GRAPH 1
#define VECTOR_BOX 2

template <typename Graph, unsigned int imp>
class VTKWriter
{

};

#include "VTKWriter_graph.hpp"
#include "VTKWriter_vector_box.hpp"

#endif /* VTKWRITER_HPP_ */
