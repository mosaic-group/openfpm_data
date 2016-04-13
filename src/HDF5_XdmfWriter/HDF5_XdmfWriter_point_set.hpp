/*
 * H5PartWriter_point_set.hpp
 *
 *  Created on: Feb 7, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_HDF5_XDMFWRITER_HDF5_XDMFWRITER_POINT_SET_HPP_
#define OPENFPM_IO_SRC_HDF5_XDMFWRITER_HDF5_XDMFWRITER_POINT_SET_HPP_

#include "HDF5_XdmfWriter_util.hpp"
#include "Vector/map_vector.hpp"
#include "VCluster.hpp"

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to produce write each property in H5Part
 *
 * \tparam ele_v is the vector of properties
 * \tparam seq, sequence of property to output
 * \tparam has_name define if the structure define names for the properties
 *
 */

template<typename ele_v, bool has_name>
struct H5_prop_out
{
	// HDF5 file
	hid_t file_id;

	// vector that we are processing
	ele_v & vv;

	// Up to which element to write
	size_t stop;

	/*! \brief constructor
	 *
	 * \param v_out string to fill with the vertex properties
	 *
	 */
	H5_prop_out(hid_t file_id, ele_v & vv, size_t stop)
	:file_id(file_id),vv(vv),stop(stop)
	{};

	//! It produce an output for each property
    template<typename T>
    void operator()(T& t) const
    {
    	typedef typename boost::mpl::at<typename ele_v::value_type::value_type::type,boost::mpl::int_<T::value>>::type ptype;

    	H5_write<ptype,T::value,ele_v>::write(file_id,std::string(ele_v::value_type::attributes::names[T::value]),vv,stop);
    }
};



/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to produce an output for each property
 *
 * \tparam ele_v is the vector of properties
 * \tparam seq, sequence of property to output
 * \tparam has_name define if the structure define names
 *
 */
template<typename ele_v>
struct H5_prop_out<ele_v,false>
{
	// HDF5 file
	hid_t file_id;

	// vector that we are processing
	ele_v & vv;

	// Up to which element to write
	size_t stop;

	/*! \brief constructor
	 *
	 * \param v_out string to fill with the vertex properties
	 *
	 */
	H5_prop_out(hid_t file_id, ele_v & vv, size_t stop)
	:file_id(file_id),vv(vv),stop(stop)
	{};

	//! It produce an output for each property
    template<typename T>
    void operator()(T& t) const
    {
    	typedef typename boost::mpl::at<typename ele_v::value_type::type,boost::mpl::int_<T::value>>::type ptype;

    	H5_write<ptype,T::value,ele_v>::write(file_id,std::string("attr") + std::to_string(T::value),vv,stop);
    }
};

template <>
class HDF5_XdmfWriter<H5_POINTSET>
{
	// Time step
	int t;

	//! HDF5 file
	hid_t file_id;

public:

	/*!
	 *
	 * H5PartWriter constructor
	 *
	 */
	HDF5_XdmfWriter()
	:t(0)
	{}


	/*!
	 *
	 * \brief Write a set of particle position and properties into HDF5
	 *
	 * \tparam Pos Vector of positions type
	 * \taparam Prp Vector of properties type
	 * \tparam prp list of properties to output
	 *
	 * \param pos Vector with the positions
	 * \param prp Vector with the properties
	 * \param stop size of the vector to output
	 *
	 */
	template<typename VPos, typename VPrp, int ... prp > bool write(const std::string & file, openfpm::vector<VPos> & v_pos, openfpm::vector<VPrp> & v_prp, size_t stop)
	{
		Vcluster & v_cl = create_vcluster();

		// Open and HDF5 file in parallel

		hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_id, v_cl.getMPIComm(), MPI_INFO_NULL);
		file_id = H5Fcreate(file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		H5Pclose(plist_id);

		// Single coordinate positional vector
		openfpm::vector<typename VPos::coord_type> x_n;
		x_n.resize(stop);

		//for each component, fill x_n
		for (size_t i = 0 ; i < VPos::dims ; i++)
		{
			//
			for (size_t j = 0 ; j < stop ; j++)
				x_n.get(j) = v_pos.template get<0>(j)[i];

			std::stringstream str;
			str << "x" << i;

			HDF5CreateDataSet<typename VPos::coord_type>(file_id,str.str(),x_n.getPointer(),stop*sizeof(typename VPos::coord_type));
		}

		// Now we write the properties

		typedef typename to_boost_vmpl<prp ... >::type v_prp_seq;
		H5_prop_out<openfpm::vector<VPrp>,has_attributes<VPrp>::value> f(file_id,v_prp,stop);

		boost::mpl::for_each_ref<v_prp_seq>(f);

		H5Fclose(file_id);

		return true;
	}
};


#endif /* OPENFPM_IO_SRC_HDF5_XDMFWRITER_HDF5_XDMFWRITER_POINT_SET_HPP_ */
