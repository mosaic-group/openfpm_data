/*
 * VTKWriter_grids_st.hpp
 *
 *  Created on: Sep 3, 2015
 *      Author: Pietro Incardona
 */

#ifndef SRC_VTKWRITER_GRIDS_ST_HPP_
#define SRC_VTKWRITER_GRIDS_ST_HPP_


#include <boost/mpl/pair.hpp>
#include "VTKWriter_grids_util.hpp"
#include "util/util_debug.hpp"
#include "util/convert.hpp"

/*! \brief for each combination in the cell grid you can have different grids
 *
 * \tparam Grid type of grid
 *
 */
template <typename Grid>
struct cell_grid
{
	//! vector of fused grids
	openfpm::vector<const Grid *> grids;

	//! combination
	//! (used to calculate the grid shift from the starting point of the cell)
	comb<Grid::dims> cmb;

	cell_grid() {}

	/*! \brief construct a cell grid
	 *
	 * \param cmb in which position this grid live
	 *
	 */
	cell_grid(const comb<Grid::dims> & cmb)
	:cmb(cmb)
	{}

	/*! \brief copy contructor
	 *
	 * \param cg element to copy
	 *
	 */
	cell_grid(const cell_grid<Grid> & cg)
	{
		this->operator=(cg);
	}

	/*! \brief copy constructor
	 *
	 * \param cg element to copy
	 *
	 */
	cell_grid(cell_grid<Grid> && cg)
	{
		this->operator=(cg);
	}

	/*! \brief Copy the cell grid
	 *
	 * \param cg cell_grid to copy
	 *
	 * \return itself
	 *
	 */
	cell_grid<Grid> & operator=(const cell_grid<Grid> & cg)
	{
		cmb = cg.cmb;
		grids = cg.grids;

		return *this;
	}

	/*! \brief Copy the cell grid
	 *
	 * \param cg cell_grid to copy
	 *
	 * \return itself
	 *
	 */
	cell_grid<Grid> & operator=(cell_grid<Grid> && cg)
	{
		cmb = cg.cmb;
		grids = cg.grids;

		return *this;
	}
};

/*! \brief convert a staggered element into a string for vtk write
 *
 * \tparam Grid type of the grid
 * \tparam St space type
 *
 */
template <typename Grid, typename St>
class ele_g_st
{
public:

	//! grid type
	typedef Grid value_type;

	//! constructor
	ele_g_st(){};

	/*! \brief convert a staggered grid property into a string
	 *
	 * \param offset shift of the staggered element
	 * \param spacing of the grid
	 * \param dom Part of the grid that is real domain
	 *
	 */
	ele_g_st(const Point<Grid::dims,St> & offset,
			 const Point<Grid::dims,St> & spacing,
			 const Box<Grid::dims,St> & dom)
	:offset(offset),spacing(spacing),dom(dom)
	{}

	//! output string
	std::string dataset;
	//! fused grids
	openfpm::vector<cell_grid<Grid>> g;
	//! offset where it start the grid
	Point<Grid::dims,St> offset;
	//! spacing of the grid
	Point<Grid::dims,St> spacing;
	//! Part of the grid that is real domain
	Box<Grid::dims,size_t> dom;

	/*! \brief Copy constructor
	 *
	 * \param ele element to copy
	 *
	 */
	inline ele_g_st(const ele_g_st & ele)
	{
		this->operator=(ele);
	}

	/*! \brief Copy constructor
	 *
	 * \param ele element to copy
	 *
	 */
	inline ele_g_st(ele_g_st && ele)
	{
		this->operator=(ele);
	}

	/*! \brief Copy the object
	 *
	 * \param ele ele_g_st to copy
	 *
	 * \return itself
	 *
	 */
	ele_g_st<Grid,St> & operator=(const ele_g_st & ele)
	{
		dataset = ele.dataset;
		g = ele.g;
		offset = ele.offset;
		spacing = ele.spacing;
		dom = ele.dom;

		return *this;
	}

	/*! \brief Copy the object
	 *
	 * \param ele ele_g_st to copy
	 *
	 * \return itself
	 *
	 */
	ele_g_st<Grid,St> & operator=(ele_g_st && ele)
	{
		dataset = ele.dataset;
		g = ele.g;
		offset = ele.offset;
		spacing = ele.spacing;
		dom = ele.dom;

		return *this;
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
class VTKWriter<pair,VECTOR_ST_GRIDS>
{
	//! Vector of grids
	openfpm::vector< ele_g_st<typename pair::first,typename pair::second> > vg;

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
			for (size_t j = 0 ; j < vg.get(i).g.size() ; j++)
			{
				if (vg.get(i).g.get(j).grids.size() != 0)
					tot += vg.get(i).g.get(j).grids.get(0)->size();
			}
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
	 * \return the list of points
	 *
	 */
	std::string get_point_list()
	{
		//! vertex node output string
		std::stringstream v_out;

		//! For each sub-domain
		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			// For each position in the cell
			for (size_t j = 0 ; j < vg.get(i).g.size() ; j++)
			{
				// If there are no grid in this position
				if (vg.get(i).g.get(j).grids.size() == 0)
					continue;

				//! Get the iterator
				auto it = vg.get(i).g.get(j).grids.get(0)->getIterator();

				//! Where the grid is defined
				Box<pair::first::dims,typename pair::second> dom;

				// Calculate the offset of the grid considering its cell position
				Point<pair::first::dims,typename pair::second> middle = vg.get(i).spacing / 2;
				Point<pair::first::dims,typename pair::second> one;
				one.one();
				one = one + toPoint<pair::first::dims,typename pair::second>::convert(vg.get(i).g.get(j).cmb);
				Point<pair::first::dims,typename pair::second> offset = pmul(middle,one) + vg.get(i).offset;

				// if there is the next element
				while (it.isNext())
				{
					Point<pair::first::dims,typename pair::second> p;
					p = it.get().toPoint();
					p = pmul(p,vg.get(i).spacing) + offset;

					if (pair::first::dims == 2)
						v_out << p.toString() << " 0.0" << "\n";
					else
						v_out << p.toString() << "\n";

					// increment the iterator
					++it;
				}
			}
		}

		// return the vertex list
		return v_out.str();
	}

	/*! \brief It generate a name for the property cell component
	 *
	 * \param k component in the cell
	 *
	 * \return property name
	 *
	 */
	std::string get_prop_components(size_t k)
	{
		std::stringstream v_out;

		//! For each sub-domain
		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			// For each position in the cell
			for (size_t j = 0 ; j < vg.get(i).g.size() ; j++)
			{
				if (k < vg.get(i).g.get(j).grids.size())
				{
					// get the combination string
					v_out << vg.get(i).g.get(j).cmb.to_string();
				}
			}
		}

		return v_out.str();
	}

	/*! \brief Create the VTK properties output
	 *
	 * \param k component
	 * \param prop_name property name
	 *
	 * \return the property output string for the grid
	 *
	 */
	std::string get_properties_output(size_t k, std::string prop_name)
	{
		//! vertex node output string
		std::stringstream v_out;

		// Check if T is a supported format
		// for now we support only scalar of native type

		typedef typename boost::mpl::at<typename pair::first::value_type::type,boost::mpl::int_<0>>::type ctype;

		std::string type = getType<ctype>();

		// if the type is not supported return
		if (type.size() == 0)
		{
#ifndef DISABLE_ALL_RTTI
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the type " << demangle(typeid(ctype).name()) << " is not supported by vtk\n";
#endif
			return "";
		}

		std::string prp_cp = get_prop_components(k);

		// Create point data properties
		v_out << "SCALARS " << prop_name << "_" << prp_cp << " " << type + "\n";

		// Default lookup table
		v_out << "LOOKUP_TABLE default\n";

		//! For each sub-domain
		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			// For each position in the cell
			for (size_t j = 0 ; j < vg.get(i).g.size() ; j++)
			{
				// If there are no grid in this position
				if (vg.get(i).g.get(j).grids.size() == 0)
					continue;

				if (k < vg.get(i).g.get(j).grids.size())
				{
					// Grid source
					auto & g_src = *vg.get(i).g.get(j).grids.get(k);

					//! Get the iterator
					auto it = g_src.getIterator();

					//! Where the grid is defined
					Box<pair::first::dims,typename pair::second> dom;

					// if there is the next element
					while (it.isNext())
					{
						auto key = it.get();

						v_out << std::to_string(g_src.template get<0>(key))  << "\n";

						// increment the iterator
						++it;
					}
				}
				else
				{
					// Grid source
					auto & g_src = *vg.get(i).g.get(j).grids.get(0);

					//! Get the iterator
					auto it = g_src.getIterator();

					//! Where the grid is defined
					Box<pair::first::dims,typename pair::second> dom;

					// if there is the next element
					while (it.isNext())
					{
						v_out << "0\n";

						// increment the iterator
						++it;
					}
				}
			}
		}

		// return the vertex list
		return v_out.str();
	}

	/*! \brief Return the output of the domain property
	 *
	 * \return vtk output
	 *
	 */
    std::string lastProp()
	{
		//! vertex node output string
		std::stringstream v_out;

		// Create point data properties
		v_out << "SCALARS domain float\n";

		// Default lookup table
		v_out << "LOOKUP_TABLE default\n";

		//! For each sub-domain
		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			// For each position in the cell
			for (size_t j = 0 ; j < vg.get(i).g.size() ; j++)
			{
				// If there are no grid in this position
				if (vg.get(i).g.get(j).grids.size() == 0)
					continue;

				//! Get the iterator
				auto it = vg.get(i).g.get(j).grids.get(0)->getIterator();

				// if there is the next element
				while (it.isNext())
				{
					if (vg.get(i).dom.isInside(it.get().toPoint()) == true)
						v_out << "1.0\n";
					else
						v_out << "0.0\n";

					// increment the iterator and counter
					++it;
				}
			}
		}

		return v_out.str();
	}

	/*! \brief Get the maximum number of fused grid
	 *
	 * \return the maximum number of fused grids
	 *
	 */
	size_t getMaxFused()
	{
		size_t max = 0;

		//! For each sub-domain
		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			// For each position in the cell
			for (size_t j = 0 ; j < vg.get(i).g.size() ; j++)
			{
				// If there are no grid in this position
				if (vg.get(i).g.get(j).grids.size() > max)
						max = vg.get(i).g.get(j).grids.size();
			}
		}

		return max;
	}

	/*! \brief Create the VTK vertex definition
	 *
	 * \return the string with the vertices as string
	 *
	 */
	std::string get_vertex_list()
	{
		//! vertex node output string
		std::string v_out;

		size_t k = 0;

		//! For each sub-domain
		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			// For each position in the cell
			for (size_t j = 0 ; j < vg.get(i).g.size() ; j++)
			{
				// If there are no grid in this position
				if (vg.get(i).g.get(j).grids.size() == 0)
						continue;
				//! For each grid point create a vertex
				auto it = vg.get(i).g.get(j).grids.get(0)->getIterator();

				while (it.isNext())
				{
					v_out += "1 " + std::to_string(k) + "\n";

					++k;
					++it;
				}
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

	/*! \brief Append the grid to the sub-domain, if for a sub-domain we have a grid that is overlapping
	 *         fuse them, otherwise create a new combination and grid
	 *
	 * \param id sub-domain id
	 * \param g grid to output
	 * \param cmb position of the grid
	 *
	 * \return a valid slot, if does not exist it append the grid at the end with the new combination
	 *
	 */
	void append_grid(size_t id, const typename pair::first & g, const comb<pair::first::dims> & cmb)
	{
		for(size_t i = 0 ; i < vg.get(id).g.size() ; i++)
		{
			// for each defined grid if exist the combination fuse
			if (cmb == vg.get(id).g.get(i).cmb)
			{
				vg.get(id).g.get(i).grids.add(&g);
				return;
			}
		}

		// if the combination does not exist add the grid
		cell_grid< typename pair::first> cg(cmb);
		vg.get(id).g.add(cg);
		vg.get(id).g.last().grids.add(&g);
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
	 * \param i sub-domain id
	 * \param g Grid to add
	 * \param offset grid offset
	 * \param spacing spacing of the grid
	 * \param dom part of the spacethat is the domain
	 * \param cmb position of the grid
	 *
	 */
	void add(size_t i,
			 const typename pair::first & g,
			 const Point<pair::first::dims,typename pair::second> & offset,
			 const Point<pair::first::dims,typename pair::second> & spacing,
			 const Box<pair::first::dims,typename pair::second> & dom,
			 const comb<pair::first::dims> & cmb)
	{
		//! Increase the size
		if (i >= vg.size())
			vg.resize(i+1);

		vg.get(i).offset = offset;
		vg.get(i).spacing = spacing;
		vg.get(i).dom = dom;

		// append the grid
		append_grid(i,g,cmb);
	}

	/*! \brief It write a VTK file from a graph
	 *
	 * \tparam prp_out which properties to output [default = -1 (all)]
	 *
	 * \param file path where to write
	 * \param g_name of the set of grids
	 * \param ft specify if it is a VTK BINARY or ASCII file [default = ASCII]
	 *
	 * \return true if the file is succeful written
	 *
	 */

	template<int prp = -1> bool write(std::string file,
			                          std::string g_name = "grids",
									  file_type ft = file_type::ASCII)
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
				     + g_name + "\n";

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

		// Get the maximum number of fused grids
		size_t mf = getMaxFused();

		// For each property in the vertex type produce a point data
		for (size_t i = 0 ; i < mf ; i++)
			point_data += get_properties_output(i,g_name);

		lastProp();


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


#endif /* SRC_VTKWRITER_GRIDS_ST_HPP_ */
