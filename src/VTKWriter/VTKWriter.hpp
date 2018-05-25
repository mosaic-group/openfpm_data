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

#define VTK_GRAPH 1
#define VECTOR_BOX 2
#define VECTOR_GRIDS 3
#define VECTOR_ST_GRIDS 4
#define DIST_GRAPH 5
#define VECTOR_POINTS 6
#define VTK_WRITER 0x10000
#define FORMAT_ASCII 0x0
#define FORMAT_BINARY 0x10000000
#define PRINT_GHOST 1

template <typename Object, unsigned int imp>
class VTKWriter
{

};

#include "VTKWriter_graph.hpp"
#include "VTKWriter_vector_box.hpp"
#include "VTKWriter_grids.hpp"
#include "VTKWriter_grids_st.hpp"
#ifndef NO_PARALLEL
#include "VTKWriter_dist_graph.hpp"
#endif
#include "VTKWriter_point_set.hpp"

#endif /* VTKWRITER_HPP_ */
