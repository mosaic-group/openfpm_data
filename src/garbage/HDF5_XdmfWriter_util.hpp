/*
 * H5PartWriteData_meta.hpp
 *
 *  Created on: Feb 22, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_HDF5_XDMFWRITER_HDF5_XDMFWRITER_UTIL_HPP_
#define OPENFPM_IO_SRC_HDF5_XDMFWRITER_HDF5_XDMFWRITER_UTIL_HPP_

#include "hdf5.h"
#include "Vector/map_vector.hpp"

/*! \brief HDF5 Create the data-set in the file
 *
 * \tparam type Type to write
 *
 * \param file_id Id of the file
 * \param filespace id where to write
 * \param str dataset to write
 * \param ptr pointer with the data to write
 * \param sz size of the data to write
 *
 * \return true if the function succeed
 *
 */
template<typename type> bool HDF5CreateDataSet(hid_t file_id, const std::string & str ,void * ptr, size_t sz)
{
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	if (plist_id < 0)
		return false;

	/* Create the dataspace for the position dataset.  */
	hsize_t dimsf[1] = {sz};
	hid_t filespace = H5Screate_simple(1, dimsf, NULL);
	if (filespace < 0)
		return false;

	if (std::is_same<type,char>::value == true)
	{
		hid_t dset_id = H5Dcreate(file_id, str.c_str(), H5T_NATIVE_CHAR, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
			return false;

		herr_t status = H5Dwrite(dset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, plist_id, ptr);
		if (status < 0)
			return false;

	    H5Dclose(dset_id);
	    H5Dclose(filespace);
	    return true;
	}
/*	else if (std::is_same<type,signed char>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_SCHAR, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dclose(dset_id);
		return status;
	}
	else if (std::is_same<type,unsigned char>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_UCHAR, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else if (std::is_same<type,short>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_SHORT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else if (std::is_same<type,unsigned short>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_USHORT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else if (std::is_same<type,int>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else if (std::is_same<type,unsigned int>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_UINT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else if (std::is_same<type,long>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_LONG, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else if (std::is_same<type,unsigned long>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_ULONG, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else if (std::is_same<type,long long>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_LLONG, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else if (std::is_same<type,unsigned long long>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_ULLONG, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else if (std::is_same<type,float>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}*/
	else if (std::is_same<type,double>::value == true)
	{
		hid_t dset_id = H5Dcreate(file_id, str.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
			return false;

		herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, ptr);
		if (status < 0)
			return false;

	    H5Dclose(dset_id);
	    H5Dclose(filespace);
	    return true;
	}
	/*else if (std::is_same<type,long double>::value == true)
	{
		dset_id = H5Dcreate(file_id, std.c_str(), H5T_NATIVE_LDOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}*/

	return true;
}


/*! \brief Write an HDF5 dataset in case of scalars and vectors
 *
 * \tparam T type to write
 * \tparam pid Property id
 * \tparam V Vector containing the information
 *
 */
template<typename T,size_t pid, typename V>
struct H5_write
{
	/*! \brief write
	 *
	 * \param file_id HDF5 file
	 * \param str dataset name
	 * \param v Vector containing the information
	 * \param stop size to store
	 *
	 */
	static inline void write(hid_t file_id, const std::string & str, V & v, size_t stop)
	{
		// Create the buffer
		openfpm::vector<T> buffer;
		buffer.resize(stop);

		for (size_t j = 0 ; j < stop ; j++)
			buffer.get(j) = v.template get<pid>(j);

		HDF5CreateDataSet<T>(file_id,str.c_str(),buffer.getPointer(),stop*sizeof(T));
	}
};

//! Partial specialization for N=1 1D-Array
template<typename T,size_t pid, typename V,size_t N1>
struct H5_write<T[N1],pid,V>
{

	/*! \brief write
	 *
	 * \param file_id HDF5 file
	 * \param str dataset name
	 * \param v Vector containing the information
	 * \param stop size to store
	 *
	 */
	static inline void write(hid_t file_id, const std::string & str, V & v, size_t stop)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			// Create the buffer
			openfpm::vector<T> buffer;
			buffer.resize(stop);

			for (size_t j = 0 ; j < stop ; j++)
				buffer.get(j) = v.template get<pid>(j)[i1];

			std::stringstream sstr;
			sstr << "_" << i1;

			HDF5CreateDataSet<T>(file_id,std::string(str) + sstr.str(),v.getPointer(),stop*sizeof(T));
		}
	}
};

//! Partial specialization for N=2 2D-Array
template<typename T, size_t pid, typename V,size_t N1,size_t N2>
struct H5_write<T[N1][N2],pid,V>
{

	/*! \brief write
	 *
	 * \param file_id HDF5 file
	 * \param str dataset name
	 * \param v Vector containing the information
	 * \param stop size to store
	 *
	 */
	static inline void write(hid_t file_id, const std::string & str, V & v, size_t stop)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				// Create the buffer
				openfpm::vector<T> buffer;
				buffer.resize(stop);

				for (size_t j = 0 ; j < stop ; j++)
					buffer.get(j) = v.template get<pid>(j)[i1][i2];

				std::stringstream sstr;
				sstr << "_" << i1 << "_" << i2;

				HDF5CreateDataSet<T>(file_id,std::string(str) + sstr.str(),v.getPointer(),stop*sizeof(T));
			}
		}
	}
};

//! Partial specialization for N=3
template<typename T, size_t pid, typename V,size_t N1,size_t N2,size_t N3>
struct H5_write<T[N1][N2][N3],pid,V>
{

	/*! \brief write
	 *
	 * \param file_id HDF5 file
	 * \param str dataset name
	 * \param v Vector containing the information
	 * \param stop size to store
	 *
	 */
	static inline void write(hid_t file_id, const std::string & str, V & v, size_t stop)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					// Create the buffer
					openfpm::vector<T> buffer;
					buffer.resize(stop);

					for (size_t j = 0 ; j < stop ; j++)
						buffer.get(j) = v.template get<pid>(j)[i1][i2][i3];

					std::stringstream sstr;
					sstr << "_" << i1 << "_" << i2 << "_" << i3;

					HDF5CreateDataSet<T>(file_id,std::string(str) + sstr.str(),v.getPointer(),stop*sizeof(T));
				}
			}
		}
	}
};

//! Partial specialization for N=4
template<typename T, size_t pid, typename V ,size_t N1,size_t N2,size_t N3,size_t N4>
struct H5_write<T[N1][N2][N3][N4],pid,V>
{

	/*! \brief write
	 *
	 * \param file_id HDF5 file
	 * \param str dataset name
	 * \param v Vector containing the information
	 * \param stop size to store
	 *
	 */
	static inline void write(hid_t file_id, const std::string & str, V & v, size_t stop)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						// Create the buffer
						openfpm::vector<T> buffer;
						buffer.resize(stop);

						for (size_t j = 0 ; j < stop ; j++)
							buffer.get(j) = v.template get<pid>(j)[i1][i2][i3][i4];


						std::stringstream sstr;
						sstr << "_" << i1 << "_" << i2 << "_" << i3 << "_" << i4;

						HDF5CreateDataSet<T>(file_id,std::string(str) + sstr.str(),v.getPointer(),stop*sizeof(T));
					}
				}
			}
		}
	}
};

#endif /* OPENFPM_IO_SRC_HDF5_XDMFWRITER_HDF5_XDMFWRITER_UTIL_HPP_ */
