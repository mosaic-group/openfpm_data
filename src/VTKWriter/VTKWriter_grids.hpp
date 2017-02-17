/*
 * VTKWriter_grids.hpp
 *
 *  Created on: May 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VTKWRITER_GRIDS_HPP_
#define VTKWRITER_GRIDS_HPP_

#include <boost/mpl/pair.hpp>
#include "VTKWriter_grids_util.hpp"
#include "is_vtk_writable.hpp"

/*! \brief It store one grid
 *
 * \tparam Grid type of grid
 * \tparam St type of space where the grid is defined
 *
 */
template <typename Grid, typename St>
class ele_g
{
public:

	typedef Grid value_type;

	ele_g(const Grid & g, const Point<Grid::dims,St> & offset, const Point<Grid::dims,St> & spacing, const Box<Grid::dims,St> & dom)
	:g(g),offset(offset),spacing(spacing),dom(dom)
	{}

	//! Dataset name
	std::string dataset;
	//! Grid
	const Grid & g;
	//! offset where it start
	Point<Grid::dims,St> offset;
	// spacing of the grid
	Point<Grid::dims,St> spacing;
	// Part of the grid that is real domain
	Box<Grid::dims,size_t> dom;
};



/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to produce at output for each property
 *
 * \tparam ele_g element that store the grid and its attributes
 * \param St type of space where the grid live
 *
 */

template<typename ele_g, typename St>
struct prop_out_g
{
	// property output string
	std::string & v_out;

	// grid that we are processing
	const openfpm::vector_std< ele_g > & vg;

	/*! \brief constructor
	 *
	 * \param v_out string to fill with the vertex properties
	 *
	 */
	prop_out_g(std::string & v_out, const openfpm::vector_std< ele_g > & vg)
	:v_out(v_out),vg(vg)
	{};

	//! It produce an output for each property
    template<typename T>
    void operator()(T& t) const
    {
    	typedef typename boost::mpl::at<typename ele_g::value_type::value_type::type,boost::mpl::int_<T::value>>::type ptype;
    	typedef typename std::remove_all_extents<ptype>::type base_ptype;

    	meta_prop<boost::mpl::int_<T::value> ,ele_g,St, ptype, is_vtk_writable<base_ptype>::value > m(vg,v_out);
    }

    void lastProp()
	{
		// Create point data properties
		v_out += "SCALARS domain float\n";

		// Default lookup table
		v_out += "LOOKUP_TABLE default\n";

		// Produce point data
		for (size_t k = 0 ; k < vg.size() ; k++)
		{
			//! Get a vertex iterator
			auto it = vg.get(k).g.getIterator();

			// if there is the next element
			while (it.isNext())
			{
				if (vg.get(k).dom.isInside(it.get().toPoint()) == true)
					v_out += "1.0\n";
				else
					v_out += "0.0\n";

				// increment the iterator and counter
				++it;
			}
		}
	}
};

/*!
 *
 * It write a VTK format file in case of grids defined on a space
 *
 * \tparam boost::mpl::pair<G,S>
 *
 * where G is the type of grid S is the type of space, float, double ...
 *
 */
template <typename pair>
class VTKWriter<pair,VECTOR_GRIDS>
{
	//! Vector of grids
	openfpm::vector< ele_g<typename pair::first,typename pair::second> > vg;

	/*! \brief Get the total number of points
	 *
	 * \return the total number
	 *
	 */
	size_t get_total()
	{
		size_t tot = 0;

		//! Calculate the full number of vertices
		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			tot += vg.get(i).g.size();
		}
		return tot;
	}

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
		v_out += "VERTICES " + std::to_string(get_total()) + " " + std::to_string(get_total() * 2) + "\n";

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
		v_out += "POINTS " + std::to_string(get_total()) + " float" + "\n";

		// return the vertex properties string
		return v_out;
	}

	/*! \brief Create the VTK point definition
	 *
	 */
	std::string get_point_list()
	{
		//! vertex node output string
		std::stringstream v_out;

		//! For each defined grid

		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			//! Get the iterator
			auto it = vg.get(i).g.getIterator();

			//! Where the grid is defined
			Box<pair::first::dims,typename pair::second> dom;

			// if there is the next element
			while (it.isNext())
			{
				Point<pair::first::dims,typename pair::second> p;
				p = it.get().toPoint();
				p = pmul(p,vg.get(i).spacing) + vg.get(i).offset;

				if (pair::first::dims == 2)
					v_out << p.toString() << " 0.0" << "\n";
				else
					v_out << p.toString() << "\n";

				// increment the iterator and counter
				++it;
			}
		}

		// return the vertex list
		return v_out.str();
	}

	/*! \brief Create the VTK vertex definition
	 *
	 */
	std::string get_vertex_list()
	{
		//! vertex node output string
		std::string v_out;

		size_t k = 0;

		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			//! For each grid point create a vertex
			auto it = vg.get(i).g.getIterator();

			while (it.isNext())
			{
				v_out += "1 " + std::to_string(k) + "\n";

				++k;
				++it;
			}
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

		v_out += "POINT_DATA " + std::to_string(get_total()) + "\n";

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

	/*! \brief Add grid dataset
	 *
	 * \param g Grid to add
	 * \param offset grid offset
	 * \param spacing spacing of the grid
	 * \param dom part of the spacethat is the domain
	 *
	 */
	void add(const typename pair::first & g, const Point<pair::first::dims,typename pair::second> & offset, const Point<pair::first::dims,typename pair::second> & spacing, const Box<pair::first::dims,typename pair::second> & dom)
	{
		ele_g<typename pair::first,typename pair::second> t(g,offset,spacing,dom);

		vg.add(t);
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

	template<int prp = -1> bool write(std::string file, std::string f_name = "grids" , file_type ft = file_type::ASCII)
	{
		// Header for the vtk
		std::string vtk_header;
		// Point list of the VTK
		std::string point_list;
		// Vertex list of the VTK
		std::string vertex_list;
		// Graph header
		std::string vtk_binary_or_ascii;
		// vertex properties header
		std::string point_prop_header;
		// edge properties header
		std::string vertex_prop_header;
		// Data point header
		std::string point_data_header;
		// Data point
		std::string point_data;

		// VTK header
		vtk_header = "# vtk DataFile Version 3.0\n"
				     + f_name + "\n";

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
		point_list = get_point_list();

		// vertex properties header
		vertex_prop_header = get_vertex_properties_list();

		// Get vertex list
		vertex_list = get_vertex_list();

		// Get the point data header
		point_data_header = get_point_data_header();

		// For each property in the vertex type produce a point data

		prop_out_g< ele_g<typename pair::first,typename pair::second>, typename pair::second > pp(point_data, vg);

		if (prp == -1)
			boost::mpl::for_each< boost::mpl::range_c<int,0, pair::first::value_type::max_prop> >(pp);
		else
			boost::mpl::for_each< boost::mpl::range_c<int,prp, prp> >(pp);

		// Add the last property
		pp.lastProp();


		// write the file
		std::ofstream ofs(file);

		// Check if the file is open
		if (ofs.is_open() == false)
		{std::cerr << "Error cannot create the VTK file: " + file + "\n";}

		ofs << vtk_header << point_prop_header << point_list <<
				vertex_prop_header << vertex_list << point_data_header << point_data;

		// Close the file

		ofs.close();

		// Completed succefully
		return true;
	}
};


#endif /* VTKWRITER_GRAPH_HPP_ */
