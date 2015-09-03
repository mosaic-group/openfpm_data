/*
 * VTKWriter_grids_util.hpp
 *
 *  Created on: Aug 10, 2015
 *      Author: Pietro Incardona
 */

#ifndef SRC_VTKWRITER_GRIDS_UTIL_HPP_
#define SRC_VTKWRITER_GRIDS_UTIL_HPP_

#include "util/util_debug.hpp"

/*! \brief This class specialize functions in the case the type T
 * has or not defined attributes
 *
 * In C++ partial specialization of a function is not allowed so we have to
 * encapsulate this function in a class
 *
 * \tparam has_attributes parameter that specialize the function in case the grid
 *         define or not attributes name
 *
 * \tparam Grid type we are processing
 * \tparam i the property we are going to write
 *
 */
template<bool has_attributes, typename St, typename ele_g, unsigned int i>
class prop_output_g
{
public:

	/*! \brief For each vertex set the value
	 *
	 * \tparam i vertex property to print
	 *
	 */
	static std::string get_point_data(const openfpm::vector< ele_g > & vg)
	{
		//! vertex node output string
		std::string v_out;

		for (size_t k = 0 ; k < vg.size() ; k++)
		{
			//! Get a vertex iterator
			auto it = vg.get(k).g.getIterator();

			// if there is the next element
			while (it.isNext())
			{
				// Print the property
				v_out += std::to_string(vg.get(k).get_o(it.get()).template get<i>()) + "\n";

				// increment the iterator and counter
				++it;
			}
		}

		return v_out;
	}

	/*! \brief Given a Graph return the point data header for a typename T
	 *
	 * \tparam T type to write
	 * \param n_node number of the node
	 *
	 */

	static std::string get_point_property_header(const std::string & oprp)
	{
		//! vertex node output string
		std::string v_out;

		typedef typename boost::fusion::result_of::at<typename ele_g::value_type::type,boost::mpl::int_<i>>::type ctype;

		// Check if T is a supported format
		// for now we support only scalar of native type

		std::string type = getType<ctype>();

		// if the type is not supported return
		// if the type is not supported return
		if (type.size() == 0)
		{
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the type " << demangle(typeid(ctype).name()) << " is not supported by vtk\n";
			return "";
		}

		// Create point data properties
		v_out += "SCALARS " + get_attributes(oprp) + " " + type + "\n";

		// Default lookup table
		v_out += "LOOKUP_TABLE default\n";

		// return the vertex list
		return v_out;
	}

	/*! \brief Get the attributes name
	 *
	 */

	static std::string get_attributes(const std::string & out)
	{
		return ele_g::value_type::attributes::name[i] + out;
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

template<typename ele_g, typename St, unsigned int i>
class prop_output_g<false,St,ele_g,i>
{
public:

	/*! \brief Get the vtk properties header appending a prefix at the end
	 *
	 * \param oprp prefix
	 *
	 */
	static std::string get_point_property_header(const std::string & oprp)
	{
		//! vertex node output string
		std::string v_out;

		// Check if T is a supported format
		// for now we support only scalar of native type

		typedef typename boost::mpl::at<typename ele_g::value_type::value_type::type,boost::mpl::int_<i>>::type ctype;
		typedef typename std::remove_all_extents<ctype>::type vttype;

		std::string type = getType<vttype>();

		// if the type is not supported return
		if (type.size() == 0)
		{
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the type " << demangle(typeid(ctype).name()) << " is not supported by vtk\n";
			return "";
		}

		// Create point data properties
		v_out += "SCALARS " + get_attributes(oprp) + " " + type + "\n";

		// Default lookup table
		v_out += "LOOKUP_TABLE default\n";

		// return the vertex list
		return v_out;
	}

	/*! \brief Get the attributes name
	 *
	 */
	static std::string get_attributes(const std::string & oprp)
	{
		return std::string("attr" + std::to_string(i) + oprp);
	}
};

/*! \brief This class is an helper to create properties output from scalar and compile-time array elements
 *
 * This class is an helper to copy scalar and compile-time array elements
 *
 */
template<typename I, typename ele_g, typename St, typename T>
struct meta_prop
{
	inline meta_prop(const openfpm::vector< ele_g > & vg, std::string & v_out)
	{
    	// actual string size
    	size_t sz = v_out.size();

		// Produce the point properties header
		v_out += prop_output_g<has_attributes<typename ele_g::value_type::value_type>::value,St ,ele_g,I::value>::get_point_property_header("");

		// If the output has changed, we have to write the properties
		if (v_out.size() != sz)
		{
			// Produce point data

			for (size_t k = 0 ; k < vg.size() ; k++)
			{
				//! Get a vertex iterator
				auto it = vg.get(k).g.getIterator();

				// if there is the next element
				while (it.isNext())
				{
					// Print the property
					v_out += std::to_string(vg.get(k).g.get_o(it.get()).template get<I::value>()) + "\n";

					// increment the iterator and counter
					++it;
				}
			}
		}
	}
};

//! Partial specialization for N=1 1D-Array
template<typename I, typename ele_g, typename St, typename T, size_t N1>
struct meta_prop<I, ele_g,St,T[N1]>
{
	inline meta_prop(const openfpm::vector< ele_g > & vg, std::string & v_out)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
	    	// actual string size
	    	size_t sz = v_out.size();

			// Produce the point properties header
			v_out += prop_output_g<has_attributes<typename ele_g::value_type::value_type>::value,St ,ele_g,I::value>::get_point_property_header("_" + std::to_string(i1));

			// If the output has changed, we have to write the properties
			if (v_out.size() != sz)
			{
				// Produce point data

				for (size_t k = 0 ; k < vg.size() ; k++)
				{
					//! Get a vertex iterator
					auto it = vg.get(k).g.getIterator();

					// if there is the next element
					while (it.isNext())
					{
						// Print the property
						v_out += std::to_string(vg.get(k).g.get_o(it.get()).template get<I::value>()[i1]) + "\n";

						// increment the iterator and counter
						++it;
					}
				}
			}
		}
	}
};

//! Partial specialization for N=2 2D-Array
template<typename I, typename ele_g, typename St ,typename T,size_t N1,size_t N2>
struct meta_prop<I, ele_g,St, T[N1][N2]>
{
	inline meta_prop(const openfpm::vector< ele_g > & vg, std::string & v_out)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
		    	// actual string size
		    	size_t sz = v_out.size();

				// Produce the point properties header
				v_out += prop_output_g<has_attributes<typename ele_g::value_type::value_type>::value,St ,ele_g,I::value>::get_point_property_header("_" + std::to_string(i1) + "_" + std::to_string(i2));

				// If the output has changed, we have to write the properties
				if (v_out.size() != sz)
				{
					// Produce point data

					for (size_t k = 0 ; k < vg.size() ; k++)
					{
						//! Get a vertex iterator
						auto it = vg.get(k).g.getIterator();

						// if there is the next element
						while (it.isNext())
						{
							// Print the property
							v_out += std::to_string(vg.get(k).g.get_o(it.get()).template get<I::value>()[i1][i2]) + "\n";

							// increment the iterator and counter
							++it;
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=3
template<typename I, typename ele_g, typename St,typename T,size_t N1,size_t N2,size_t N3>
struct meta_prop<I,ele_g,St,T[N1][N2][N3]>
{
	inline meta_prop(const openfpm::vector< ele_g > & vg, std::string & v_out)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
			    	// actual string size
			    	size_t sz = v_out.size();

					// Produce the point properties header
					v_out += prop_output_g<has_attributes<typename ele_g::value_type::value_type>::value,St ,ele_g,I::value>::get_point_property_header("_" + std::to_string(i1) + "_" + std::to_string(i2) + "_" + std::to_string(i3));

					// If the output has changed, we have to write the properties
					if (v_out.size() != sz)
					{
						std::string attr = prop_output_g<has_attributes<typename ele_g::value_type::value_type>::value, St,ele_g,I::value>::get_attributes() + "_";

						// Produce point data

						for (size_t k = 0 ; k < vg.size() ; k++)
						{
							//! Get a vertex iterator
							auto it = vg.get(k).g.getIterator();

							// if there is the next element
							while (it.isNext())
							{
								// Print the property
								v_out += std::to_string(vg.get(k).g.get_o(it.get()).template get<I::value>()[i1][i2][i3]) + "\n";

								// increment the iterator and counter
								++it;
							}
						}
					}
				}
			}
		}
	}
};


#endif /* SRC_VTKWRITER_GRIDS_UTIL_HPP_ */
