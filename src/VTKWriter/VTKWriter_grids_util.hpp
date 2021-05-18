/*
 * VTKWriter_grids_util.hpp
 *
 *  Created on: Aug 10, 2015
 *      Author: Pietro Incardona
 */

#ifndef SRC_VTKWRITER_GRIDS_UTIL_HPP_
#define SRC_VTKWRITER_GRIDS_UTIL_HPP_

#include "util/util_debug.hpp"
#include "is_vtk_writable.hpp"
#include "byteswap_portable.hpp"

/*! \brief Return the Attributes name from the type
 *
 *
 */
template<typename ele_g, bool has_attributes>
struct getAttrName
{
	/*! \brief Get attribute name
	 *
	 * \param i id of the attribute
	 * \param oprp post-fix to add
	 * \param prop_names properties name
	 *
	 * \return the string with the property name
	 *
	 */
	static inline std::string get(size_t i, const openfpm::vector<std::string> & prop_names, const std::string & oprp)
	{
		return ele_g::value_type::value_type::attributes::name[i] + oprp;
	}
};

/*! \brief Return the Attributes name from the type
 *
 * This is the case when the type has no name
 *
 */
template<typename ele_g>
struct getAttrName<ele_g,false>
{
	/*! \brief Get attribute name
	 *
	 * \param i id of the attribute
	 * \param oprp post-fix to add
	 * \param  prop_names properties name
	 *
	 * \return return the string with the property name
	 *
	 */
	static inline std::string get(size_t i, const openfpm::vector<std::string> & prop_names, const std::string & oprp)
	{
		if (i >= prop_names.size())
		{
			return std::string("attr" + std::to_string(i) + oprp);
		}
		else
		{
			return prop_names.get(i) + oprp;
		}
	}
};

/*! \brief Get the vtk properties header appending a prefix at the end
 *
 * \tparam has_attributes indicate if the properties have attributes name
 * \param oprp prefix
 *
 */
template<unsigned int i, typename ele_g, bool has_attributes> std::string get_point_property_header_impl(const std::string & oprp, const openfpm::vector<std::string> & prop_names)
{
    //! vertex node output string
    std::string v_out;

    typedef typename boost::mpl::at<typename ele_g::value_type::value_type::type,boost::mpl::int_<i>>::type ctype;

    // Check if T is a supported format
    // for now we support only scalar of native type
    if (std::rank<ctype>::value == 1)
    {
        if (std::extent<ctype>::value <= 3)
        {
            //Get type of the property
            std::string type = getType<typename std::remove_all_extents<ctype>::type>();

            // if the type is not supported skip-it
            if (type.size() == 0)
            {
#ifndef DISABLE_ALL_RTTI
                std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the type " << demangle(typeid(ctype).name()) << " is not supported by vtk\n";
#endif
                return "";
            }

            // Create point data properties
            v_out += "VECTORS " + getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp) + " " + type + "\n";
        }
    }
    else
    {
        std::string type = getType<typename std::remove_all_extents<ctype>::type>();

        // if the type is not supported return
        if (type.size() == 0)
        {
            // We check if is a custom vtk writable object

            if (is_vtk_writable<ctype>::value == true)
            {
                type = getType<typename vtk_type<ctype,is_custom_vtk_writable<ctype>::value>::type >();

                // We check if it is a vector or scalar like type
                if (vtk_dims<ctype>::value == 1)
                    v_out += "SCALARS " + getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp) + " " + type + "\n";
                else
                    v_out += "VECTORS " + getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp) + " " + type + "\n";
            }

            return v_out;
        }

        // Create point data properties
        v_out += "SCALARS " + getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp) + " " + type + "\n";

        // Default lookup table
        v_out += "LOOKUP_TABLE default\n";

    }

    // return the vertex list
    return v_out;
}


/*! \brief Get the vtk properties header appending a prefix at the end
 *
 * \tparam has_attributes indicate if the properties have attributes name
 * \param oprp prefix
 *
 */
/*template<unsigned int i, typename ele_g, bool has_attributes> std::string get_point_property_header_impl_new(const std::string & oprp, const openfpm::vector<std::string> & prop_names,file_type ft)
{
	//! vertex node output string
	std::string v_out;

	typedef typename boost::mpl::at<typename ele_g::value_type::value_type::type,boost::mpl::int_<i>>::type ctype;

	// Check if T is a supported format
	// for now we support only scalar of native type
	if (std::rank<ctype>::value == 1)
	{
		if (std::extent<ctype>::value <= 3)
		{
			//Get type of the property
			std::string type = getTypeNew<typename std::remove_all_extents<ctype>::type>();

			// if the type is not supported skip-it
			if (type.size() == 0)
			{
#ifndef DISABLE_ALL_RTTI
				std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the type " << demangle(typeid(ctype).name()) << " is not supported by vtk\n";
#endif
				return "";
			}

			// Create point data properties
			v_out += "        <DataArray type=\""+type+"\" Name=\""+getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp)+"\""+" NumberOfComponents=\"3\"";
			if(ft==file_type::ASCII){
                v_out+=" format=\"ascii\">\n";
			}
			else{
                v_out+=" format=\"binary\">\n";
            }

			//v_out += "VECTORS " + getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp) + " " + type + "\n";
		}
	}
	else
	{
		std::string type = getTypeNew<typename std::remove_all_extents<ctype>::type>();

		// if the type is not supported return
		if (type.size() == 0)
		{
			// We check if is a custom vtk writable object

			if (is_vtk_writable<ctype>::value == true)
			{
				type = getTypeNew<typename vtk_type<ctype,is_custom_vtk_writable<ctype>::value>::type >();

				// We check if it is a vector or scalar like type
				if (vtk_dims<ctype>::value == 1) {
                    v_out += "        <DataArray type=\"" + type + "\" Name=\"" +
                             getAttrName<ele_g, has_attributes>::get(i, prop_names, oprp) + "\"";
                    if (ft == file_type::ASCII) {
                        v_out += " format=\"ascii\">\n";
                    } else {
                        v_out += " format=\"binary\">\n";
                    }
                    //v_out += "SCALARS " + getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp) + " " + type + "\n";
                }
				else{
                    v_out += "        <DataArray type=\""+type+"\" Name=\""+getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp)+"\""+" NumberOfComponents=\"3\"";
                    if(ft==file_type::ASCII){
                            v_out+=" format=\"ascii\">\n";
                        }
                    else{
                            v_out+=" format=\"binary\">\n";
                        }
					//v_out += "VECTORS " + getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp) + " " + type + "\n";
			    }
			}

			return v_out;
		}

		// Create point data properties
        v_out += "        <DataArray type=\"" + type + "\" Name=\"" +getAttrName<ele_g, has_attributes>::get(i, prop_names, oprp) + "\"";
        if (ft == file_type::ASCII) {
            v_out += " format=\"ascii\">\n";
        }
        else {
            v_out += " format=\"binary\">\n";
        }
		//v_out += "SCALARS " + getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp) + " " + type + "\n";
		// Default lookup table
		//v_out += "LOOKUP_TABLE default\n";

	}

	// return the vertex list
	return v_out;
}*/

/*! \brief Write the vectror property
 *
 * \tparam dim Dimensionality of the property
 *
 */
template<unsigned int dim, typename T>
class prop_write_out
{
public:

	template<typename vector, typename iterator, typename I> static void write(std::ostringstream & v_out, vector & vg, size_t k, iterator & it, file_type ft)
	{
		if (ft == file_type::ASCII)
		{
			// Print the properties
			for (size_t i1 = 0 ; i1 < vtk_dims<T>::value ; i1++)
			{
				v_out << vg.get(k).g.get_o(it.get()).template get<I::value>().get_vtk(i1) << " ";
			}
			if (vtk_dims<T>::value == 2)
			{
				v_out << "0.0";
			}
			v_out << "\n";
		}
		else
		{
			typedef decltype(vg.get(k).g.get_o(it.get()).template get<I::value>().get_vtk(0)) ctype_;
			typedef typename std::remove_reference<ctype_>::type ctype;

			// Print the properties
			for (size_t i1 = 0 ; i1 < vtk_dims<T>::value ; i1++)
			{
				typename is_vtk_writable<ctype>::base tmp = vg.get(k).g.get_o(it.get()).template get<I::value>().get_vtk(i1);
				tmp = swap_endian_lt(tmp);
				v_out.write((const char *)&tmp,sizeof(tmp));
			}
			if (vtk_dims<T>::value == 2)
			{
				typename is_vtk_writable<ctype>::base zero = 0.0;
				zero = swap_endian_lt(zero);
				v_out.write((const char *)&zero,sizeof(zero));
			}
		}
	}
};

/*! \brief Write the scalar property
 *
 *
 */
template<typename T>
class prop_write_out<1,T>
{
public:

	/*! \brief Write the property
	 *
	 *  \param v_out output string of the property
	 *  \param vg vector of properties
	 *  \param k data-set to output
	 *  \param it iterator with the point to output
	 *  \param ft output type BINARY or ASCII
	 *
	 */
	template<typename vector, typename iterator, typename I> static void write(std::ostringstream & v_out, vector & vg, size_t k, iterator & it, file_type ft)
	{
		typedef decltype(vg.get(k).g.template get<I::value>(it.get())) ctype_;
		typedef typename std::remove_const<typename std::remove_reference<ctype_>::type>::type ctype;

		if (ft == file_type::ASCII)
		{
			// Print the property
			v_out << vg.get(k).g.template get<I::value>(it.get()) << "\n";
		}
		else
		{
			typename is_vtk_writable<ctype>::base tmp = vg.get(k).g.template get<I::value>(it.get());
			tmp = swap_endian_lt(tmp);
			v_out.write((const char *)&tmp,sizeof(tmp));
		}
	}
};

/*! \brief Write the vectror property
 *
 * \tparam dim Dimensionality of the property
 *
 */
template<unsigned int dim, typename T>
class prop_write_out_new
{
public:

    template<typename vector, typename iterator, typename I> static void write(std::ostringstream & v_out, vector & vg, size_t k, iterator & it, file_type ft)
    {

        if (ft == file_type::ASCII)
        {
            // Print the properties
            for (size_t i1 = 0 ; i1 < vtk_dims<T>::value ; i1++)
            {
                v_out << vg.get(k).g.get_o(it.get()).template get<I::value>().get_vtk(i1) << " ";
            }
            if (vtk_dims<T>::value == 2)
            {
                v_out << "0.0";
            }
            v_out << "\n";
        }
        else
        {
            typedef decltype(vg.get(k).g.get_o(it.get()).template get<I::value>().get_vtk(0)) ctype_;
            typedef typename std::remove_reference<ctype_>::type ctype;

            // Print the properties
            for (size_t i1 = 0 ; i1 < vtk_dims<T>::value ; i1++)
            {
                typename is_vtk_writable<ctype>::base tmp = vg.get(k).g.get_o(it.get()).template get<I::value>().get_vtk(i1);
                v_out.write((const char *)&tmp,sizeof(tmp));
            }
            if (vtk_dims<T>::value == 2)
            {
                typename is_vtk_writable<ctype>::base zero = 0.0;
                v_out.write((const char *)&zero,sizeof(zero));
            }
        }
    }
};

/*! \brief Write the scalar property
 *
 *
 */
template<typename T>
class prop_write_out_new<1,T>
{
public:

    /*! \brief Write the property
     *
     *  \param v_out output string of the property
     *  \param vg vector of properties
     *  \param k data-set to output
     *  \param it iterator with the point to output
     *  \param ft output type BINARY or ASCII
     *
     */
    template<typename vector, typename iterator, typename I> static void write(std::ostringstream & v_out, vector & vg, size_t k, iterator & it, file_type ft)
    {
        typedef decltype(vg.get(k).g.template get<I::value>(it.get())) ctype_;
        typedef typename std::remove_const<typename std::remove_reference<ctype_>::type>::type ctype;

        if (ft == file_type::ASCII)
        {
            // Print the property
            v_out << vg.get(k).g.template get<I::value>(it.get()) << "\n";
        }
        else
        {
            typename is_vtk_writable<ctype>::base tmp = vg.get(k).g.template get<I::value>(it.get());
            v_out.write((const char *)&tmp,sizeof(tmp));
        }
    }
};


/*! \brief This class is an helper to create properties output from scalar and compile-time array elements
 *
 * \tparam I It is an boost::mpl::int_ that indicate which property we are writing
 * \tparam ele_g element type that store the grid information
 * \tparam St type of space where the grid live
 * \tparam T the type of the property
 * \tparam is_writable flag that indicate if a property is writable
 *
 */
template<typename I, typename ele_g, typename St, typename T, bool is_writable>
struct meta_prop
{
    /*! \brief Write a vtk compatible type into vtk format
     *
     * \param vg array of elements to write
     * \param v_out string containing the string
     * \param prop_names property names
     * \param ft ASCII or BINARY
     *
     */
    inline meta_prop(const openfpm::vector< ele_g > & vg, std::string & v_out, const openfpm::vector<std::string> & prop_names, file_type ft)
    {
        // actual string size
        size_t sz = v_out.size();

        // Produce the point properties header
        v_out += get_point_property_header_impl<I::value,ele_g,has_attributes<typename ele_g::value_type::value_type>::value>("",prop_names);

        // If the output has changed, we have to write the properties
        if (v_out.size() != sz)
        {
            std::ostringstream stream_out;
            if (std::is_same<T,float>::value == true)
            {stream_out << std::setprecision(7);}
            else
            {stream_out << std::setprecision(16);}

            // Produce point data

            for (size_t k = 0 ; k < vg.size() ; k++)
            {
                //! Get a vertex iterator
                auto it = vg.get(k).g.getIterator();

                // if there is the next element
                while (it.isNext())
                {
                    prop_write_out<vtk_dims<T>::value,T>::template write<decltype(vg),decltype(it),I>(stream_out,vg,k,it,ft);

                    // increment the iterator and counter
                    ++it;
                }
            }

            v_out += stream_out.str();

            if (ft == file_type::BINARY)
                v_out += "\n";
        }
    }
};

//! Partial specialization for N=1 1D-Array
template<typename I, typename ele_g, typename St, typename T, size_t N1, bool is_writable>
struct meta_prop<I, ele_g,St,T[N1],is_writable>
{
    /*! \brief Write a vtk compatible type into vtk format
     *
     * \param vg array of elements to write
     * \param v_out string containing the string
     * \param prop_names properties name
     * \param ft ASCII or BINARY
     *
     */
    inline meta_prop(const openfpm::vector< ele_g > & vg, std::string & v_out, const openfpm::vector<std::string> & prop_names , file_type ft)
    {
        // actual string size
        size_t sz = v_out.size();

        // Produce the point properties header
        v_out += get_point_property_header_impl<I::value,ele_g,has_attributes<typename ele_g::value_type::value_type>::value>("",prop_names);

        // If the output has changed, we have to write the properties
        if (v_out.size() != sz)
        {
            std::ostringstream stream_out;
            if (std::is_same<T,float>::value == true)
            {stream_out << std::setprecision(7);}
            else
            {stream_out << std::setprecision(16);}

            // Produce point data

            for (size_t k = 0 ; k < vg.size() ; k++)
            {
                //! Get a vertex iterator
                auto it = vg.get(k).g.getIterator();

                // if there is the next element
                while (it.isNext())
                {
                    if (ft == file_type::ASCII)
                    {
                        // Print the properties
                        for (size_t i1 = 0 ; i1 < N1 ; i1++)
                        {stream_out << vg.get(k).g.template get<I::value>(it.get())[i1] << " ";}

                        if (N1 == 2)
                        {stream_out << (decltype(vg.get(k).g.template get<I::value>(it.get())[0])) 0;}

                        stream_out << "\n";
                    }
                    else
                    {
                        T tmp;

                        // Print the properties
                        for (size_t i1 = 0 ; i1 < N1 ; i1++)
                        {
                            tmp = vg.get(k).g.template get<I::value>(it.get())[i1];
                            tmp = swap_endian_lt(tmp);
                            stream_out.write((const char *)&tmp,sizeof(T));
                        }
                        if (N1 == 2)
                        {
                            tmp = 0.0;
                            tmp = swap_endian_lt(tmp);
                            stream_out.write((const char *)&tmp,sizeof(T));
                        }
                    }

                    // increment the iterator and counter
                    ++it;
                }
            }

            v_out += stream_out.str();

            if (ft == file_type::BINARY)
            {v_out += "\n";}
        }
    }
};

//! Partial specialization for N=2 2D-Array
template<typename I, typename ele_g, typename St ,typename T,size_t N1,size_t N2, bool is_writable>
struct meta_prop<I, ele_g,St, T[N1][N2],is_writable>
{

    /*! \brief Write a vtk compatible type into vtk format
     *
     * \param vg array of elements to write
     * \param v_out string containing the string
     * \param prop_names property names
     * \param ft ASCII or BINARY
     *
     */
    inline meta_prop(const openfpm::vector< ele_g > & vg, std::string & v_out, const openfpm::vector<std::string> & prop_names, file_type ft)
    {
        for (size_t i1 = 0 ; i1 < N1 ; i1++)
        {
            for (size_t i2 = 0 ; i2 < N2 ; i2++)
            {
                std::ostringstream stream_out;
                if (std::is_same<T,float>::value == true)
                {stream_out << std::setprecision(7);}
                else
                {stream_out << std::setprecision(16);}

                // actual string size
                size_t sz = v_out.size();

                // Produce the point properties header
                v_out += get_point_property_header_impl<I::value,ele_g,has_attributes<typename ele_g::value_type::value_type>::value>("_" + std::to_string(i1) + "_" + std::to_string(i2),prop_names);

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
                            T tmp;

                            if (ft == file_type::ASCII)
                            {
                                // Print the property
                                stream_out << vg.get(k).g.template get<I::value>(it.get())[i1][i2] << "\n";
                            }
                            else
                            {
                                tmp = vg.get(k).g.template get<I::value>(it.get())[i1][i2];
                                tmp = swap_endian_lt(tmp);
                                stream_out.write((const char *)&tmp,sizeof(tmp));
                            }

                            // increment the iterator and counter
                            ++it;
                        }
                    }

                    v_out += stream_out.str();

                    if (ft == file_type::BINARY)
                        v_out += "\n";
                }
            }
        }
    }
};


//! Specialication when is not writable
template<typename I, typename ele_g, typename St, typename T>
struct meta_prop<I,ele_g,St,T,false>
{

    /*! \brief Write a vtk compatible type into vtk format
     *
     * \param vg array of elements to write
     * \param v_out string containing the string
     * \param prop_names properties name
     * \param ft ASCII or BINARY
     *
     */
    inline meta_prop(const openfpm::vector< ele_g > & vg, std::string & v_out, const openfpm::vector<std::string> & prop_names, file_type ft)
    {
    }
};


/*! \brief This class is an helper to create properties output from scalar and compile-time array elements
 *
 * \tparam I It is an boost::mpl::int_ that indicate which property we are writing
 * \tparam ele_g element type that store the grid information
 * \tparam St type of space where the grid live
 * \tparam T the type of the property
 * \tparam is_writable flag that indicate if a property is writable
 *
 */
#if 0
template<typename I, typename ele_g, typename St, typename T, bool is_writable>
struct meta_prop_new
{
	/*! \brief Write a vtk compatible type into vtk format
	 *
	 * \param vg array of elements to write
	 * \param v_out string containing the string
	 * \param prop_names property names
	 * \param ft ASCII or BINARY
	 *
	 */
	inline meta_prop_new(const openfpm::vector< ele_g > & vg, std::string & v_out, const openfpm::vector<std::string> & prop_names, file_type ft)
	{
    	// actual string size
    	size_t sz = v_out.size();
    	std::ostringstream v_outToEncode_;
        std::string v_Encoded;

        if (std::is_same<T,float>::value == true)
        {v_outToEncode_ << std::setprecision(7);}
        else
        {v_outToEncode_ << std::setprecision(16);}

		// Produce the point properties header
    	v_out += get_point_property_header_impl_new<I::value,ele_g,has_attributes<typename ele_g::value_type::value_type>::value>("",prop_names,ft);

		// If the output has changed, we have to write the properties
		if (v_out.size() != sz)
		{
            if (ft == file_type::BINARY) {
                size_t pp;
                v_outToEncode_.write((char *)&pp,8);
            }
			// Produce point data
			for (size_t k = 0 ; k < vg.size() ; k++)
			{
				//! Get a vertex iterator
				auto it = vg.get(k).g.getIterator();

				// if there is the next element
				while (it.isNext())
				{
					prop_write_out_new<vtk_dims<T>::value,T>::template write<decltype(vg),decltype(it),I>(v_outToEncode_,vg,k,it,ft);

					// increment the iterator and counter
					++it;
				}
			}

            if (ft == file_type::BINARY)
            {
                std::string v_outToEncode = v_outToEncode_.str();
                *(size_t *) &v_outToEncode[0] = v_outToEncode.size()-sizeof(size_t);
                v_Encoded.resize(v_outToEncode.size()/3*4+4);
                size_t sz=EncodeToBase64((const unsigned char*)&v_outToEncode[0],v_outToEncode.size(),(unsigned char *)&v_Encoded[0],0);
                v_Encoded.resize(sz);
                v_out += v_Encoded + "\n";
            }
            else{
                v_out += v_outToEncode_.str();
            };
			v_out += "        </DataArray>\n";

		}
	}
};

//! Partial specialization for N=1 1D-Array
template<typename I, typename ele_g, typename St, typename T, size_t N1, bool is_writable>
struct meta_prop_new<I, ele_g,St,T[N1],is_writable>
{
	/*! \brief Write a vtk compatible type into vtk format
	 *
	 * \param vg array of elements to write
	 * \param v_out string containing the string
	 * \param prop_names properties name
	 * \param ft ASCII or BINARY
	 *
	 */
	inline meta_prop_new(const openfpm::vector< ele_g > & vg, std::string & v_out, const openfpm::vector<std::string> & prop_names , file_type ft)
	{
	    // actual string size
	    size_t sz = v_out.size();
        std::ostringstream v_outToEncode_;
	    std::string v_Encoded;

        if (std::is_same<T,float>::value == true)
        {v_outToEncode_ << std::setprecision(7);}
        else
        {v_outToEncode_ << std::setprecision(16);}

		// Produce the point properties header
		v_out += get_point_property_header_impl_new<I::value,ele_g,has_attributes<typename ele_g::value_type::value_type>::value>("",prop_names,ft);

		// If the output has changed, we have to write the properties
		if (v_out.size() != sz)
		{
            if (ft == file_type::BINARY) {
                size_t pp = 0;
                v_outToEncode_.write((char *)&pp,8);
            }
			// Produce point data

			for (size_t k = 0 ; k < vg.size() ; k++)
			{
				//! Get a vertex iterator
				auto it = vg.get(k).g.getIterator();

				// if there is the next element
				while (it.isNext())
				{
					if (ft == file_type::ASCII)
					{
						// Print the properties
						for (size_t i1 = 0 ; i1 < N1 ; i1++)
						{v_outToEncode_ << vg.get(k).g.template get<I::value>(it.get())[i1] << " ";}

						if (N1 == 2)
						{v_outToEncode_ << (decltype(vg.get(k).g.template get<I::value>(it.get())[0])) 0;}

						v_outToEncode_ << "\n";
					}
					else
					{
						T tmp;

						// Print the properties
						for (size_t i1 = 0 ; i1 < N1 ; i1++)
						{
							tmp = vg.get(k).g.template get<I::value>(it.get())[i1];
							//tmp = swap_endian_lt(tmp);
							v_outToEncode_.write((const char *)&tmp,sizeof(T));
						}
						if (N1 == 2)
						{
							tmp = 0.0;
							//tmp = swap_endian_lt(tmp);
							v_outToEncode_.write((const char *)&tmp,sizeof(T));
						}
					}

					// increment the iterator and counter
					++it;
				}
			}
            if (ft == file_type::BINARY)
            {
                std::string v_outToEncode = v_outToEncode_.str();
                *(size_t *) &v_outToEncode[0] = v_outToEncode.size()-sizeof(size_t);
                v_Encoded.resize(v_outToEncode.size()/3*4+4);
                size_t sz=EncodeToBase64((const unsigned char*)&v_outToEncode[0],v_outToEncode.size(),(unsigned char *)&v_Encoded[0],0);
                v_Encoded.resize(sz);
                v_out += v_Encoded + "\n";
            }
            else{
                v_out += v_outToEncode_.str();
            };
            v_out += "        </DataArray>\n";
        }
	}
};

//! Partial specialization for N=2 2D-Array
template<typename I, typename ele_g, typename St ,typename T,size_t N1,size_t N2, bool is_writable>
struct meta_prop_new<I, ele_g,St, T[N1][N2],is_writable>
{

	/*! \brief Write a vtk compatible type into vtk format
	 *
	 * \param vg array of elements to write
	 * \param v_out string containing the string
	 * \param prop_names property names
	 * \param ft ASCII or BINARY
	 *
	 */
	inline meta_prop_new(const openfpm::vector< ele_g > & vg, std::string & v_out, const openfpm::vector<std::string> & prop_names, file_type ft)
	{
        std::string v_outToEncode,v_Encoded;

        for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
                v_outToEncode.clear();
		    	// actual string size
		    	size_t sz = v_out.size();

				// Produce the point properties header
				v_out += get_point_property_header_impl_new<I::value,ele_g,has_attributes<typename ele_g::value_type::value_type>::value>("_" + std::to_string(i1) + "_" + std::to_string(i2),prop_names, ft);

				// If the output has changed, we have to write the properties
				if (v_out.size() != sz)
				{
                    if (ft == file_type::BINARY) {
                        v_outToEncode.append(8,0);
                    }
					// Produce point data

					for (size_t k = 0 ; k < vg.size() ; k++)
					{
						//! Get a vertex iterator
						auto it = vg.get(k).g.getIterator();

						// if there is the next element
						while (it.isNext())
						{
							T tmp;

							if (ft == file_type::ASCII)
							{
								// Print the property
								v_outToEncode += std::to_string(vg.get(k).g.template get<I::value>(it.get())[i1][i2]) + "\n";
							}
							else
							{
								tmp = vg.get(k).g.template get<I::value>(it.get())[i1][i2];
								//tmp = swap_endian_lt(tmp);
								v_outToEncode.append((const char *)&tmp,sizeof(tmp));
							}

							// increment the iterator and counter
							++it;
						}
					}
                    if (ft == file_type::BINARY)
                    {
                        *(size_t *) &v_outToEncode[0] = v_outToEncode.size()-sizeof(size_t);
                        v_Encoded.resize(v_outToEncode.size()/3*4+4);
                        size_t sz=EncodeToBase64((const unsigned char*)&v_outToEncode[0],v_outToEncode.size(),(unsigned char *)&v_Encoded[0],0);
                        v_Encoded.resize(sz);
                        v_out += v_Encoded + "\n";
                    }
                    else{
                        v_out += v_outToEncode;
                    };
                    v_out += "        </DataArray>\n";

                }
			}
		}
	}
};


//! Specialication when is not writable
template<typename I, typename ele_g, typename St, typename T>
struct meta_prop_new<I,ele_g,St,T,false>
{

	/*! \brief Write a vtk compatible type into vtk format
	 *
	 * \param vg array of elements to write
	 * \param v_out string containing the string
	 * \param prop_names properties name
	 * \param ft ASCII or BINARY
	 *
	 */
	inline meta_prop_new(const openfpm::vector< ele_g > & vg, std::string & v_out, const openfpm::vector<std::string> & prop_names, file_type ft)
	{
	}
};

#endif

template<unsigned int dims,typename T> inline void output_point(Point<dims,T> & p,std::stringstream & v_out, file_type ft)
{
	if (ft == file_type::ASCII)
	{
		v_out << p.toString();
		size_t i = dims;
		for ( ; i < 3 ; i++)
		{v_out << " 0.0";}
		v_out << "\n";
	}
	else
	{
		size_t i = 0;
		for ( ; i < dims ; i++)
		{
			// we use float so we have to convert to float
			auto tmp = p.get(i);
			//tmp = swap_endian_lt(tmp);
			v_out.write((const char *)&tmp,sizeof(tmp));
		}
		for ( ; i < 3 ; i++)
		{
			// we use float so we have to convert to float

			/* coverity[dead_error_begin] */
			T tmp = 0.0;
			//tmp = swap_endian_lt(tmp);
			v_out.write((const char *)&tmp,sizeof(tmp));
		}
	}
}


template<unsigned int dims,typename T> inline void output_point_new(Point<dims,T> & p,std::stringstream & v_out, file_type ft)
{
    if (ft == file_type::ASCII)
    {
        v_out << p.toString();
        size_t i = dims;
        for ( ; i < 3 ; i++)
        {v_out << " 0.0";}
        v_out << "\n";
    }
    else
    {
        size_t i = 0;
        for ( ; i < dims ; i++)
        {
            // we use float so we have to convert to float
            auto tmp = p.get(i);
            v_out.write((const char *)&tmp,sizeof(tmp));
        }
        for ( ; i < 3 ; i++)
        {
            // we use float so we have to convert to float

            /* coverity[dead_error_begin] */
            T tmp = 0.0;
            v_out.write((const char *)&tmp,sizeof(tmp));
        }
    }
}

inline void output_vertex(size_t k,std::string & v_out, file_type ft)
{
    if (ft == file_type::ASCII)
    {v_out += "1 " + std::to_string(k) + "\n";}
    else
    {
        int tmp;
        tmp = 1;
        tmp = swap_endian_lt(tmp);
        v_out.append((const char *)&tmp,sizeof(int));
        tmp = k;
        tmp = swap_endian_lt(tmp);
        v_out.append((const char *)&tmp,sizeof(int));
    }
}

inline void output_vertex_new(size_t k,std::string & v_out, file_type ft)
{
    if (ft == file_type::ASCII)
    {v_out += std::to_string(k) + "\n";}
    else
    {
        size_t tmp;
        tmp = k;
        v_out.append((const char *)&tmp,sizeof(size_t));
    }
}

#endif /* SRC_VTKWRITER_GRIDS_UTIL_HPP_ */
