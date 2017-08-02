/*
 * VTKWriter_vector_box.hpp
 *
 *  Created on: May 5, 2015
 *      Author: i-bird
 */

#ifndef VTKWRITER_VECTOR_BOX_HPP_
#define VTKWRITER_VECTOR_BOX_HPP_

#include <boost/math/special_functions/pow.hpp>
#include "Space/Shape/HyperCube.hpp"
#include <random>
#include "util/util.hpp"

template <typename vector>
class v_box
{
public:

	v_box(const vector & v)
	:v(v)
	{}

	std::string dataset;
	const vector & v;
};

/*!
 *
 * From a basic structure it write a VTK format file in case of a vector of Boxes
 *
 * \tparam vector type
 *
 */

template <typename vector>
class VTKWriter<vector,VECTOR_BOX>
{
	openfpm::vector<v_box<vector>> v;

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

		// number of points
		size_t np = 0;

		// count the number of points
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			np += v.get(i).v.size() * boost::math::pow<vector::value_type::dims>(2);
		}

		// write the number of vertex
		v_out += "POINTS " + std::to_string(np) + " float" + "\n";

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

	std::string get_cell_properties_list()
	{
		//! vertex property output string
		std::string e_out;

		//! number of cells and box
		size_t nc = 0;
		size_t nb = 0;

		// count the number of cells
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			nb += v.get(i).v.size();
		}

		nc = nb * (boost::math::pow<vector::value_type::dims>(2) + 1);

		// write the number of lines
		e_out += "CELLS " + std::to_string(nb) + " " + std::to_string(nc) + "\n";

		// return the vertex properties string
		return e_out;
	}

	/*! \brief Create the VTK point definition
	 *
	 * \tparam s_type spatial type of the data
	 * \tparam attr false x,y,z are set to 0 for each vertex
	 *
	 */

	std::string get_point_list()
	{
		//! vertex node output string
		std::string v_out;

		//! for each vertex dataset

		for (size_t i = 0 ; i < v.size() ; i++)
		{
			auto it = v.get(i).v.getIterator();

			// if there is the next element
			while (it.isNext())
			{
				// Get the box
				auto box = v.get(i).v.get(it.get());

				// Create an hyper-cube and get the vertex combinations

				HyperCube<vector::value_type::dims> hyp;
				std::vector<comb<vector::value_type::dims>> comb =  hyp.getCombinations_R(0);

				// Create the box vertex points

				for (size_t j = 0; j < comb.size() ; j++)
				{
					Point<vector::value_type::dims,float> p;

					for (size_t k = 0 ; k < 3 ; k++)
					{
						if (k < vector::value_type::dims)
						{
							if (comb[j].value(k) < 0)
								v_out += std::to_string(box.template get<vector::value_type::p1>()[k]) + " ";
							else
								v_out += std::to_string(box.template get<vector::value_type::p2>()[k]) + " ";
						}
						else
						{
							v_out += "0.0";
						}
					}
					v_out += "\n";
				}

				// increment the iterator and counter
				++it;
			}
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

	std::string get_cell_list()
	{
		// base
		size_t base = 0;

		//! vertex node output string
		std::string v_out;

		//! for each vector in the dataset

		for (size_t i = 0 ; i < v.size() ; i++)
		{
			auto it = v.get(i).v.getIterator();

			// for each box
			while (it.isNext())
			{
				// Output the box vertex id
				v_out += std::to_string((size_t)boost::math::pow<vector::value_type::dims>(2)) + " ";
				for (size_t k = 0 ; k < boost::math::pow<vector::value_type::dims>(2) ; k++)
				{
					v_out += " " + std::to_string(base+k);
				}
				base += boost::math::pow<vector::value_type::dims>(2);
				v_out += "\n";

				++it;
			}
			v_out += "\n";
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

		// number of points
		size_t np = 0;

		// count the number of points
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			np += v.get(i).v.size() * boost::math::pow<vector::value_type::dims>(2);
		}


		v_out += "POINT_DATA " + std::to_string(np) + "\n";

		return v_out;
	}

	std::string get_cell_types_header()
	{
		//! vertex property output string
		std::string e_out;

		//! number of cells and box
		size_t nb = 0;

		// count the number of cells
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			nb += v.get(i).v.size();
		}

		// write the number of lines
		e_out += "CELL_TYPES " + std::to_string(nb) + "\n";

		// return the vertex properties string
		return e_out;
	}

	std::string get_cell_types_list()
	{
		// Cell id
		size_t cell_id;
		if (vector::value_type::dims == 2)
			cell_id = 8;
		else
			cell_id = 11;

		//! vertex node output string
		std::string v_out;

		//! for each vector in the dataset

		for (size_t i = 0 ; i < v.size() ; i++)
		{
			auto it = v.get(i).v.getIterator();

			// for each box
			while (it.isNext())
			{
				v_out += std::to_string(cell_id) + "\n";

				++it;
			}
		}

		// return the vertex list
		return v_out;
	}


	std::string get_cell_data_header()
	{
		//! vertex property output string
		std::string e_out;

		//! number of cells and box
		size_t nb = 0;

		// count the number of cells
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			nb += v.get(i).v.size();
		}

		// write the number of lines
		e_out += "CELL_DATA " + std::to_string(nb) + "\n";
		e_out += "COLOR_SCALARS data 4\n";

		// return the vertex properties string
		return e_out;
	}

	std::string get_cell_data_list()
	{
		// random engine
		SimpleRNG rng;

		//! vertex node output string
		std::string v_out;

		size_t col_group = 0;

		//! for each vector in the dataset
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			auto it = v.get(i).v.getIterator();

			// for each box
			while (it.isNext())
			{
				// write a color
				v_out += getColor(col_group,rng).toString() + " 1.0" + "\n";

				++it;
			}
			v_out += "\n";
			col_group++;
		}

		// return the vertex list
		return v_out;
	}

public:

	/*!
	 *
	 * VTKWriter constructor
	 *
	 */
	VTKWriter()
	{}

	/*! \brief Add box vector dataset
	 *
	 * \param v vector to add
	 *
	 */
	void add(const vector & vc)
	{
		v_box<vector> t(vc);

		v.add(t);
	}

	/*! \brief It write a VTK file from a graph
	 *
	 * \tparam prp_out which properties to output [default = -1 (all)]
	 *
	 * \param file path where to write
	 * \param name of the graph
	 * \param ft specify if it is a VTK BINARY or ASCII file [default = ASCII]
	 *
	 */

	template<int prp = -1> bool write(std::string file, std::string graph_name="Graph", file_type ft = file_type::ASCII)
	{
		// Header for the vtk
		std::string vtk_header;
		// Point list of the VTK
		std::string point_list;
		// Vertex list of the VTK
		std::string cell_list;
		// Graph header
		std::string vtk_binary_or_ascii;
		// Edge list of the GraphML
		std::string edge_list;
		// vertex properties header
		std::string point_prop_header;
		// edge properties header
		std::string cell_prop_header;
		// edge properties header
		std::string edge_prop_header;
		// Data point header
		std::string point_data_header;
		// Data point
		std::string point_data;
		// Cell type header
		std::string cell_types_header;
		// Cell type list
		std::string cell_types_list;
		// Cell data header
		std::string cell_data_header;
		// Cell data list
		std::string cell_data_list;

		// VTK header
		vtk_header = "# vtk DataFile Version 3.0\n"
				     + graph_name + "\n";

		// Choose if binary or ASCII
		if (ft == file_type::ASCII)
		{vtk_header += "ASCII\n";}
		else
		{vtk_header += "BINARY\n";}

		// Data type for graph is DATASET POLYDATA
		vtk_header += "DATASET UNSTRUCTURED_GRID\n";

		// point properties header
		point_prop_header = get_point_properties_list();

		// Get point list
		point_list = get_point_list();

		// cell properties header
		cell_prop_header = get_cell_properties_list();

		// Get cell list
		cell_list = get_cell_list();

		// Get cell types
		cell_types_header = get_cell_types_header();

		// Get cell type list
		cell_types_list = get_cell_types_list();

		// Get cell data header
		cell_data_header = get_cell_data_header();

		// Get cell data list
		cell_data_list = get_cell_data_list();

		// write the file
		std::ofstream ofs(file);

		// Check if the file is open
		if (ofs.is_open() == false)
		{std::cerr << "Error cannot create the VTK file: " + file + "\n";}

		ofs << vtk_header << point_prop_header << point_list <<
				cell_prop_header << cell_list << cell_types_header << cell_types_list << cell_data_header << cell_data_list;

		// Close the file

		ofs.close();

		// Completed succefully
		return true;
	}
};



#endif /* VTKWRITER_VECTOR_BOX_HPP_ */
