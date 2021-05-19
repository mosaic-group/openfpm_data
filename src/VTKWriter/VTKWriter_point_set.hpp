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
#include "MetaParser/MetaParser.hpp"

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

    //! properties names
    const openfpm::vector<std::string> & prop_names;

    /*! \brief constructor
     *
     * \param v_out string to fill with the vertex properties
     * \param vv vector we are processing
     * \param ft ASCII or BINARY format
     *
     */
    prop_out_v(std::string & v_out,
               const openfpm::vector_std< ele_v > & vv,
               const openfpm::vector<std::string> & prop_names,
               file_type ft)
            :ft(ft),v_out(v_out),vv(vv),prop_names(prop_names)
    {};

    /*! \brief It produce an output for each property
     *
     * \param t property id
     *
     */
    template<typename T>
    void operator()(T& t) const
    {
        typedef typename boost::mpl::at<typename ele_v::value_type::value_type::type,boost::mpl::int_<T::value>>::type ptype;
        typedef typename std::remove_all_extents<ptype>::type base_ptype;

        meta_prop_new<boost::mpl::int_<T::value> ,ele_v,St, ptype, is_vtk_writable<base_ptype>::value > m(vv,v_out,prop_names,ft);
    }

    void lastProp()
    {
        std::string v_outToEncode,v_Encoded;
        // Create point data properties
        //v_out += "SCALARS domain float\n";
        // Default lookup table
        //v_out += "LOOKUP_TABLE default\n";
        v_out += "        <DataArray type=\"Float32\" Name=\"domain\"";
        if (ft == file_type::ASCII) {
            v_out += " format=\"ascii\">\n";
        }
        else {
            v_out += " format=\"binary\">\n";
        }

        if (ft == file_type::BINARY) {
            v_outToEncode.append(8,0);
        }
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
                        v_outToEncode += "1.0\n";
                    else
                        v_outToEncode += "0.0\n";
                }
                else
                {
                    if (it.get() < vv.get(k).mark)
                    {
                        float one = 1;
                        //one = swap_endian_lt(one);
                        v_outToEncode.append((const char *)&one,sizeof(int));
                    }
                    else
                    {
                        float zero = 0;
                        //zero = swap_endian_lt(zero);
                        v_outToEncode.append((const char *)&zero,sizeof(int));
                    }
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
        v_out+="        </DataArray>\n";
    }

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
struct prop_out_v_pvtp
{
    //! property output string
    std::string & v_out;

    //! properties names
    const openfpm::vector<std::string> & prop_names;

    /*! \brief constructor
     *
     * \param v_out string to fill with the vertex properties
     * \param vv vector we are processing
     * \param ft ASCII or BINARY format
     *
     */
    prop_out_v_pvtp(std::string & v_out,
                    const openfpm::vector<std::string> & prop_names)
            :v_out(v_out),prop_names(prop_names)
    {
        //meta_prop_new<boost::mpl::int_<T::value> ,ele_v,St, ptype, is_vtk_writable<base_ptype>::value > m(vv,v_out,prop_names,ft);
    };

    /*! \brief It produce an output for each property
     *
     * \param t property id
     *
     */
    template<typename T>
    void operator()(T& t) const
    {
        typedef typename boost::mpl::at<typename ele_v::value_type::value_type::type,boost::mpl::int_<T::value>>::type ptype;
        typedef typename std::remove_all_extents<ptype>::type base_ptype;

        //std::string type = getTypeNew<base_ptype>();
        meta_prop_new<boost::mpl::int_<T::value> ,ele_v,St, ptype, is_vtk_writable<base_ptype>::value >::get_pvtp_out(v_out,prop_names);
        //v_out += "    <PDataArray type=\""+type+"\" Name=\""+getAttrName<ele_g,has_attributes>::get(i,prop_names,oprp)+"\""+" NumberOfComponents=\"3\"";
    }


    void lastProp()
    {
        v_out += "      <PDataArray type=\"Float32\" Name=\"domain\"/>\n    </PPointData>\n";
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
    std::string get_vertex_properties_list(file_type & opt)
    {
        //! vertex property output string
        std::string v_out;

        v_out += "      <Verts>\n";
        if (opt == file_type::ASCII)

        {
            v_out+="        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        }
        else
        {
            v_out+="        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"binary\">\n";
        }

        // write the number of vertex
        //v_out += "VERTICES " + std::to_string(get_total()) + " " + std::to_string(get_total() * 2) + "\n";
        // return the vertex properties string
        return v_out;
    }

    /*! \brief It get the point position header string
     *
     * It get the vertex position header of the vertex defined as a VTK header
     *
     * \return a string that define the vertex position format
     *
     */
    std::string get_point_properties_list(file_type ft)
    {
        //! vertex property output string
        std::string v_out;

        // write the number of vertex

        v_out += "    <Piece NumberOfPoints=\"" + std::to_string(get_total()) + "\" " +"NumberOfVerts=\"" + std::to_string(get_total()) + "\">\n";

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
    std::string get_point_list(file_type & opt)
    {
        //! vertex node output string
        std::stringstream v_out;

       v_out<<"      <Points>\n";

        if (std::is_same<float,typename pair::first::value_type::coord_type>::value == true)
        {
            if (opt == file_type::ASCII)
            {
                v_out<<"        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            }
            else
            {
                v_out<<"        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"binary\">\n";
            }
        }
        else
        {
            if (opt == file_type::ASCII)
            {
                v_out<<"        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            }
            else
            {
                v_out<<"        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"binary\">\n";
            }
        }

        std::stringstream binaryToEncode;
        if (std::is_same<float,typename pair::first::value_type::coord_type>::value == true)
        {
            binaryToEncode << std::setprecision(7);
        }
        else
        {
            binaryToEncode << std::setprecision(16);
        }

        //! For each defined grid
        if (opt == file_type::BINARY)
        {
            size_t tmp=0;
            binaryToEncode.write((const char *)&tmp,sizeof(tmp));
        }

        for (size_t i = 0 ; i < vps.size() ; i++)
        {
            //! write the particle position
            auto it = vps.get(i).g.getIterator();

            // if there is the next element
            while (it.isNext())
            {
                Point<pair::first::value_type::dims,typename pair::first::value_type::coord_type> p;
                p = vps.get(i).g.get(it.get());

                output_point_new<pair::first::value_type::dims,typename pair::first::value_type::coord_type>(p,binaryToEncode,opt);

                // increment the iterator and counter
                ++it;
            }
        }
        //! In case of binary we have to add a new line at the end of the list
        if (opt == file_type::BINARY){
            std::string buffer_out,buffer_bin;
            buffer_bin=binaryToEncode.str();
            *(size_t *)&buffer_bin[0]=buffer_bin.size()-8;
            buffer_out.resize(buffer_bin.size()/3*4+4);
            unsigned long sz = EncodeToBase64((const unsigned char*)&buffer_bin[0],buffer_bin.size(),(unsigned char*)&buffer_out[0],0);
            buffer_out.resize(sz);
            v_out << buffer_out<<std::endl;
        }
        else
        {
            v_out<<binaryToEncode.str();
        }
        v_out<<"        </DataArray>\n";
        v_out<<"      </Points>\n";
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
        std::string v_out,v_outToEncode,v_Encoded;

        size_t k = 0;
        if (ft == file_type::BINARY) {
            v_outToEncode.append(8,0);
        }
        for (size_t i = 0 ; i < vps.size() ; i++)
        {
            //! For each grid point create a vertex
            auto it = vps.get(i).g.getIterator();

            while (it.isNext())
            {
                output_vertex_new(k,v_outToEncode,ft);

                ++k;
                ++it;
            }
        }
        //! In case of binary we have to add a new line at the end of the list
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
        v_out += "                <DataArray type=\"Int64\" Name=\"offsets\" ";

        if (ft == file_type::ASCII)
        {
            v_out += "format=\"ascii\">\n";
        }
        else{
            v_out += "format=\"binary\">\n";
        }

        k=0;
        v_outToEncode.clear();
        if (ft == file_type::BINARY) {
            v_outToEncode.append(8,0);
        }

        for (size_t i = 0 ; i < vps.size() ; i++)
        {
            //! For each grid point create a vertex
            auto it = vps.get(i).g.getIterator();
            while (it.isNext())
            {
                output_vertex_new(k+1,v_outToEncode,ft);

                ++k;
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
        v_out += "      </Verts>\n";
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

        v_out += "      <PointData>\n";

        return v_out;
    }
    struct doubleint{
        long int i;
        double d;
    };
    /*! \brief return the meta data string
     *
     * \param meta_data string with the meta-data to add
     *
     */
    std::string add_meta_data(std::string & meta_data, file_type & opt)
    {
        std::string meta_string;

        // check for time metadata
        MetaParser_options opts;
        opts.add_options()
                ("time", MetaParser_def::value<double>());

        MetaParser mp(opts);
        mp.parse(meta_data);

        double time = 0.0;
        bool exist = mp.getOption("time",time);

        if (exist == true)
        {
            //<DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii" RangeMin="2" RangeMax="2">
            //meta_string += "";
            meta_string += "    <FieldData>\n";

            if (opt == file_type::ASCII)
            {   meta_string += "        <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">\n";
                meta_string += std::to_string(time);
            }
            else
            {
                meta_string += "        <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"binary\">\n";

                //time = swap_endian_lt(time);
                unsigned char time_string[24];//= base64_encode((const unsigned char*)&time,6);
                //const unsigned char Time=(const unsigned char)time;
                doubleint timeInit;
                timeInit.i=8;
                timeInit.d=time;
                size_t sz=EncodeToBase64((const unsigned char*)&timeInit,16,time_string,0);
                //meta_string.append((const char *)&time,sizeof(double));
                //meta_string += time_string;
                meta_string.append((const char *)time_string,sz);
            }
            meta_string += "\n";
            meta_string += "      </DataArray>\n";
            meta_string += "    </FieldData>\n";
        }


        return meta_string;
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
     * \param mark additional information that divide the dataset into 2 (in general is used to mark real from ghost information)
     * \param opt_names optional parameter that indicate the names of the properties
     *
     */
    void add(const typename pair::first & vps,
             const typename pair::second & vpp,
             size_t mark)
    {
        ele_vps<typename pair::first> t1(vps,mark);
        ele_vpp<typename pair::second> t2(vpp,mark);

        this->vps.add(t1);
        this->vpp.add(t2);
    }

/*! \brief It write a Merged VTP type file from a vector of points
	 *
	 * \tparam prp_out which properties to output [default = -1 (all)]
	 *
	 * \return true if the write complete successfully
	 *
	 */
    bool write_pvtp(std::string file,const openfpm::vector<std::string> & prop_names,size_t n,size_t timestamp=-1)
    {
        //openfpm::vector< ele_vpp<typename pair::second>> vpp;
        // Header for the vtk
        std::string vtk_header;
        std::string Name_data;
        std::string PpointEnd;
        std::string Piece;

        vtk_header = "<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PPolyData>\n    <PPointData>\n";
        prop_out_v_pvtp< ele_vpp<typename pair::second>, typename pair::first::value_type::coord_type> pp(Name_data,prop_names);
        boost::mpl::for_each< boost::mpl::range_c<int,0, pair::second::value_type::max_prop> >(pp);
        pp.lastProp();
        PpointEnd += "    <PPoints>\n      <PDataArray type=\""+getTypeNew<typename decltype(vps)::value_type::value_type::value_type::coord_type>()+"\" Name=\"Points\" NumberOfComponents=\"3\"/>\n    </PPoints>\n";


        if (timestamp==-1) {
            for (int i = 0; i < n; i++)
            { Piece += "    <Piece Source=\"" + file.substr(0, file.size()) + "_" +std::to_string(i) + ".vtp\"/>\n";}
            file +=  ".pvtp";
        }
        else{
            for (int i = 0; i < n; i++)
            { Piece += "    <Piece Source=\"" + file.substr(0, file.size()) + "_" +std::to_string(i) + "_" + std::to_string(timestamp) + ".vtp\"/>\n";}
            file += "_" + std::to_string(timestamp) + ".pvtp";
        }
        std::string closingFile="  </PPolyData>\n</VTKFile>";

        // write the file
        std::ofstream ofs(file);

        // Check if the file is open
        if (ofs.is_open() == false)
        {std::cerr << "Error cannot create the PVTP file: " + file + "\n";}

        ofs << vtk_header << Name_data <<PpointEnd<< Piece << closingFile;

        // Close the file

        ofs.close();

        return true;
    }

    /*! \brief It write a VTK file from a vector of points
     *
     * \tparam prp_out which properties to output [default = -1 (all)]
     *
     * \param file path where to write
     * \param f_name name of the dataset
     * \param prop_names properties names
     * \param ft specify if it is a VTK BINARY or ASCII file [default = ASCII]
     *
     * \return true if the write complete successfully
     *
     */
    template<int prp = -1> bool write(std::string file,
                                      const openfpm::vector<std::string> & prop_names,
                                      std::string f_name = "points" ,
                                      std::string meta_data = "",
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
        vtk_header = "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";

        vtk_header +="  <PolyData>\n";

        // Choose if binary or ASCII
/*		if (ft == file_type::ASCII)
		{vtk_header += "ASCII\n";}
		else
		{vtk_header += "BINARY\n";}*/

        // Data type for graph is DATASET POLYDATA
        //vtk_header += "DATASET POLYDATA\n";

        vtk_header += add_meta_data(meta_data,ft);

        // point properties header
        point_prop_header = get_point_properties_list(ft);

        // Get point list
        point_list = get_point_list(ft);

        // vertex properties header
        vertex_prop_header = get_vertex_properties_list(ft);

        // Get vertex list
        vertex_list = get_vertex_list(ft);

        // Get the point data header
        point_data_header = get_point_data_header();

        // For each property in the vertex type produce a point data

        prop_out_v< ele_vpp<typename pair::second>, typename pair::first::value_type::coord_type> pp(point_data, vpp, prop_names,ft);

        if (prp == -1)
        {boost::mpl::for_each< boost::mpl::range_c<int,0, pair::second::value_type::max_prop> >(pp);}
        else
        {boost::mpl::for_each< boost::mpl::range_c<int,prp, prp> >(pp);}

        // Add the last property
        pp.lastProp();

        std::string closingFile="      </PointData>\n    </Piece>\n  </PolyData>\n</VTKFile>";

        // write the file
        std::ofstream ofs(file);

        // Check if the file is open
        if (ofs.is_open() == false)
        {std::cerr << "Error cannot create the VTK file: " + file + "\n";}

        ofs << vtk_header << point_prop_header << point_list <<
            vertex_prop_header << vertex_list << point_data_header << point_data << closingFile;

        // Close the file

        ofs.close();

        // Completed succefully
        return true;
    }
};


#endif /* OPENFPM_IO_SRC_VTKWRITER_POINT_SET_HPP_ */
