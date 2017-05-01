/*
 * VTKWriter_point_set.hpp
 *
 *  Created on: Feb 6, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_VTKWRITER_POINT_SET_HPP_
#define OPENFPM_IO_SRC_VTKWRITER_POINT_SET_HPP_

#include <cstddef>
#include <boost/mpl/pair.hpp>
#include "VTKWriter_grids_util.hpp"
#include "is_vtk_writable.hpp"
#include <string>
#include "byteswap_portable.hpp"

/*! \brief Store a reference to the vector position
 *
 * \tparam Vps Type of vector that store the position of the particles
 *
 */
template <typename Vps>
class ele_vps
{
public:

	//! type of vector that store the particle position
	typedef Vps value_type;

	//! particle position vector
	const Vps & g;

	//! ghost marker
	size_t mark;

	//! constructor
	ele_vps(const Vps & g, size_t mark)
	:g(g),mark(mark)
	{}

};

/*! \brief Store a reference to the vector properties
 *
 * \tparam Vpp Type of vector that store the property of the particles
 *
 */
template <typename Vpp>
class ele_vpp
{
public:

	//! type of vector that store the particle properties
	typedef Vpp value_type;


	//! Reference to the particle properties
	const Vpp & g;

	//! ghost marker
	size_t mark;

	//! constructor
	ele_vpp(const Vpp & vpp, size_t mark)
	:g(vpp),mark(mark)
	{}

};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to produce an output for each property
 *
 * \tparam ele_v It is the class ele_v that store the couple vector of position and property
 *
 *
 */
template<typename ele_v, typename St>
struct prop_out_v
{
	//! Binary or ASCII
	file_type ft;

	//! property output string
	std::string & v_out;

	//! vector that we are processing
	const openfpm::vector_std< ele_v > & vv;

	/*! \brief constructor
	 *
	 * \param v_out string to fill with the vertex properties
	 * \param vv vector we are processing
	 * \param ft ASCII or BINARY format
	 *
	 */
	prop_out_v(std::string & v_out, const openfpm::vector_std< ele_v > & vv, file_type ft)
	:ft(ft),v_out(v_out),vv(vv)
	{};

	//! It produce an output for each property
    template<typename T>
    void operator()(T& t) const
    {
    	typedef typename boost::mpl::at<typename ele_v::value_type::value_type::type,boost::mpl::int_<T::value>>::type ptype;
    	typedef typename std::remove_all_extents<ptype>::type base_ptype;

    	meta_prop<boost::mpl::int_<T::value> ,ele_v,St, ptype, is_vtk_writable<base_ptype>::value > m(vv,v_out,ft);
    }

    void lastProp()
	{
		// Create point data properties
		v_out += "SCALARS domain float\n";

		// Default lookup table
		v_out += "LOOKUP_TABLE default\n";

		// Produce point data
		for (size_t k = 0 ; k < vv.size() ; k++)
		{
			//! Get a vertex iterator
			auto it = vv.get(k).g.getIterator();

			// if there is the next element
			while (it.isNext())
			{
				if (ft == file_type::ASCII)
				{
					if (it.get() < vv.get(k).mark)
						v_out += "1.0\n";
					else
						v_out += "0.0\n";
				}
				else
				{
					if (it.get() < vv.get(k).mark)
					{
						int one = 1;
						one = swap_endian_lt(one);
						v_out.append((const char *)&one,sizeof(int));
					}
					else
					{
						int zero = 0;
						zero = swap_endian_lt(zero);
						v_out.append((const char *)&zero,sizeof(int));
					}
				}

				// increment the iterator and counter
				++it;
			}
		}
	}
};

/*!
 *
 * It write a VTK format file for a list of grids defined on a space
 *
 * \tparam boost::mpl::pair<G,S>
 *
 * where G is the type of the vector containing the properties, S is the
 * type of vector containing the particle positions
 *
 */
template <typename pair>
class VTKWriter<pair,VECTOR_POINTS>
{
	//! Vector of position
	openfpm::vector< ele_vps<typename pair::first >> vps;
	//! Vector of properties
	openfpm::vector< ele_vpp<typename pair::second>> vpp;

	/*! \brief Get the total number of points
	 *
	 * \return the total number
	 *
	 */
	size_t get_total()
	{
		size_t tot = 0;

		//! Calculate the full number of vertices
		for (size_t i = 0 ; i < vps.size() ; i++)
		{
			tot += vps.get(i).g.size();
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

	/*! \brief Create the VTK point list
	 *
	 * \param ft file_type
	 *
	 * \return the list of points
	 *
	 */
	std::string get_point_list(file_type ft)
	{
		//! vertex node output string
		std::stringstream v_out;

		//! For each defined grid

		for (size_t i = 0 ; i < vps.size() ; i++)
		{
			//! write the particle position
			auto it = vps.get(i).g.getIterator();

			// if there is the next element
			while (it.isNext())
			{
				Point<pair::first::value_type::dims,typename pair::first::value_type::coord_type> p;
				p = vps.get(i).g.get(it.get());

				output_point<pair::first::value_type::dims,typename pair::first::value_type::coord_type>(p,v_out,ft);

				// increment the iterator and counter
				++it;
			}
		}

		//! In case of binary we have to add a new line at the end of the list
		if (ft == file_type::BINARY)
			v_out << std::endl;

		// return the vertex list
		return v_out.str();
	}

	/*! \brief Create the VTK vertex list
	 *
	 * \param ft file_type
	 *
	 * \return the list of vertices
	 *
	 */
	std::string get_vertex_list(file_type ft)
	{
		// vertex node output string
		std::string v_out;

		size_t k = 0;

		for (size_t i = 0 ; i < vps.size() ; i++)
		{
			//! For each grid point create a vertex
			auto it = vps.get(i).g.getIterator();

			while (it.isNext())
			{
				output_vertex(k,v_out,ft);

				++k;
				++it;
			}
		}

		//! In case of binary we have to add a new line at the end of the list
		if (ft == file_type::BINARY)
			v_out += "\n";

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

	/*! \brief Add a vector dataset
	 *
	 * \param vps vector of positions
	 * \param vpp vector of properties
	 * \param mark additional information that divide the dataset into 2
	 *        (in general is used to mark real from ghost information)
	 *
	 */
	void add(const typename pair::first & vps, const typename pair::second & vpp,size_t mark)
	{
		ele_vps<typename pair::first> t1(vps,mark);
		ele_vpp<typename pair::second> t2(vpp,mark);

		this->vps.add(t1);
		this->vpp.add(t2);
	}

	/*! \brief It write a VTK file from a vector of points
	 *
	 * \tparam prp_out which properties to output [default = -1 (all)]
	 *
	 * \param file path where to write
	 * \param f_name name of the dataset
	 * \param ft specify if it is a VTK BINARY or ASCII file [default = ASCII]
	 *
	 */
	template<int prp = -1> bool write(std::string file, std::string f_name = "points" , file_type ft = file_type::ASCII)
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
		point_list = get_point_list(ft);

		// vertex properties header
		vertex_prop_header = get_vertex_properties_list();

		// Get vertex list
		vertex_list = get_vertex_list(ft);

		// Get the point data header
		point_data_header = get_point_data_header();

		// For each property in the vertex type produce a point data

		prop_out_v< ele_vpp<typename pair::second>, typename pair::first::value_type::coord_type> pp(point_data, vpp, ft);

		if (prp == -1)
			boost::mpl::for_each< boost::mpl::range_c<int,0, pair::second::value_type::max_prop> >(pp);
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


#endif /* OPENFPM_IO_SRC_VTKWRITER_POINT_SET_HPP_ */
