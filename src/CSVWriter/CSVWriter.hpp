/*
 * CSVWriter.hpp
 *
 *  Created on: Dec 15, 2014
 *      Author: Pietro Incardona
 */

#ifndef CSVWRITER_HPP_
#define CSVWRITER_HPP_

#include <iostream>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <fstream>
#include "util/common.hpp"
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>
#include "csv_multiarray.hpp"
#include "util/util.hpp"
#include "is_csv_writable.hpp"

#define CSV_WRITER 0x30000

/*! \brief this class is a functor for "for_each" algorithm
 *
 * For each element of the boost::vector the operator() is called.
 * Is mainly used to create a string containing all the properties of the object
 *
 * \tparam Tobj object
 *
 */

template<typename Tobj>
struct csv_prp
{
	//! String containing the csv line constructed from an object
	std::stringstream & str;

	//! Object to write
	Tobj & obj;

	/*! \brief Constructor
	 *
	 * Create a vertex properties list
	 *
	 * \param str streamstring
	 * \param obj object to write
	 *
	 */
	csv_prp(std::stringstream & str, Tobj & obj)
	:str(str),obj(obj)
	{
	};

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
		// This is the type of the csv column
		typedef decltype(obj.template get<T::value>()) col_type;

		// Remove the reference from the column type
		typedef typename boost::remove_reference<col_type>::type col_rtype;
		typedef typename std::remove_all_extents<col_rtype>::type base_col_rtype;

    	csv_value_str<col_rtype, is_csv_writable<base_col_rtype>::value >(obj.template get<T::value>(),str);
    }
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * For each element of the boost::vector the operator() is called.
 * Is mainly used to create a string containing all the properties of the object
 *
 * \tparam T object
 *
 */

template<typename Tobj, bool attr>
struct csv_col
{
	//! String containing the colums list as string
	std::stringstream & str;

	/*! \brief Constructor
	 *
	 * \str String where to put the colum list
	 *
	 */
	csv_col(std::stringstream & str)
	:str(str)
	{
	};

	//! It call the functor for each member
    template<typename T>
    inline void operator()(T& t)
    {
		// This is the type of the csv column
		typedef typename boost::mpl::at<typename Tobj::type,boost::mpl::int_<T::value>>::type col_type;

    	csv_col_str<col_type>(std::string(Tobj::attributes::name[T::value]),str);
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

template<typename Tobj>
struct csv_col<Tobj,false>
{
	//! String containing the colums list as string
	std::stringstream & str;

	/*! \brief Constructor
	 *
	 * \str String where to put the colum list
	 *
	 */
	csv_col(std::stringstream & str)
	:str(str)
	{
	};

	//! It call the functor for each member
    template<typename T>
    void operator()(T& t)
    {
		// This is the type of the csv column
		typedef typename boost::fusion::result_of::at_c<typename Tobj::type,T::value>::type col_type;

		// Remove the reference from the column type
		typedef typename boost::remove_reference<col_type>::type col_rtype;

		std::stringstream str2;
		str2 << "column_" << T::value;

		csv_col_str<col_rtype>(str2.str(),str);
    }
};

//#define VECTOR 1

/*! \brief CSV Writer
 *
 * It write in CSV format vector of objects living into an N-dimensional space
 *
 * \tparam v_pos Positional vector
 * \tparam v_prp Property vector
 *
 */
template <typename v_pos, typename v_prp, unsigned int impl = 1>
class CSVWriter
{
	/*! \brief Get the colums name (also the positional name)
	 *
	 */
	std::string get_csv_colums()
	{
		std::stringstream str;

		// write positional columns
		for (size_t i = 0 ; i < v_pos::value_type::dims ; i++)
		{
			if (i == 0)
				str << "x[" << i << "]";
			else
				str << "," << "x[" << i << "]";
		}

		// write positional information

		csv_col<typename v_prp::value_type,has_attributes<typename v_prp::value_type>::value> col(str);

		// Iterate through all the vertex and create the vertex list
		boost::mpl::for_each< boost::mpl::range_c<int,0,v_prp::value_type::max_prop> >(col);

		str << "\n";

		return str.str();
	}

	/*! \brief Get the csv data section
	 *
	 * \param vp vector that contain the positional information
	 * \param vpr vector that contain the property information
	 * \param offset from where to start
	 *
	 */
	std::string get_csv_data(v_pos & vp, v_prp & vpr, size_t offset)
	{
		std::stringstream str;

		// The position and property vector size must match
		if (vp.size() != vpr.size())
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " position vector and property vector must have the same size \n";
			return std::string("");
		}

		// Write the data
		for (size_t i = offset ; i < vp.size() ; i++)
		{
			for (size_t j = 0 ; j < v_pos::value_type::dims ; j++)
			{
				if (j == 0)
					str << vp.template get<0>(i)[j];
				else
					str << "," << vp.template get<0>(i)[j];
			}

			// Object to write
			auto obj = vpr.get(i);

			csv_prp<decltype(obj)> c_prp(str,obj);

			// write the properties to the stream string
			boost::mpl::for_each< boost::mpl::range_c<int,0,v_prp::value_type::max_prop> >(c_prp);

			str << "\n";
		}

		return str.str();
	}

public:

	/*! \brief It write a CSV file
	 *
	 * \tparam prp which properties to output [default = -1 (all)]
	 *
	 * \param file path where to write
	 * \param v positional vector
	 * \param prp properties vector
	 * \param offset from where to start to write
	 *
	 */
	bool write(std::string file, v_pos & v , v_prp & prp, size_t offset=0)
	{
		// Header for csv (colums name)
		std::string csv_header;
		// Data point
		std::string point_data;

		// Get csv columns
		csv_header = get_csv_colums();

		// For each property in the vertex type produce a point data
		point_data = get_csv_data(v,prp,offset);

		// write the file
		std::ofstream ofs(file);

		// Check if the file is open
		if (ofs.is_open() == false)
		{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " cannot create the CSV file: " << file << std::endl;}

		ofs << csv_header << point_data;

		// Close the file

		ofs.close();

		// Completed succefully
		return true;
	}
};


#endif /* CSVWRITER_HPP_ */

