/*
 * HDF5_writer_gd.hpp
 *
 *  Created on: May 2, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_HDF5_WR_HDF5_WRITER_GD_HPP_
#define OPENFPM_IO_SRC_HDF5_WR_HDF5_WRITER_GD_HPP_


#include "Packer_Unpacker/Pack_selector.hpp"
#include "Packer_Unpacker/Packer.hpp"
#include "Packer_Unpacker/Unpacker.hpp"
#include "util/GBoxes.hpp"

template <>
class HDF5_writer<GRID_DIST>
{
public:

	template<typename device_grid>
	inline void save(const std::string & filename,
			         const openfpm::vector<device_grid> & loc_grid,
					 const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext) const
	{
		Vcluster & v_cl = create_vcluster();

		//Pack_request vector
		size_t req = 0;

		//Pack request
		Packer<typename std::remove_reference<decltype(loc_grid)>::type,HeapMemory>::packRequest(loc_grid,req);
		Packer<typename std::remove_reference<decltype(gdb_ext)>::type,HeapMemory>::packRequest(gdb_ext,req);

		//std::cout << "Req: " << req << std::endl;

		// allocate the memory
		HeapMemory pmem;
		//pmem.allocate(req);
		ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(req,pmem));
		mem.incRef();

		//Packing

		Pack_stat sts;

		Packer<typename std::remove_reference<decltype(loc_grid)>::type,HeapMemory>::pack(mem,loc_grid,sts);
		Packer<typename std::remove_reference<decltype(gdb_ext)>::type,HeapMemory>::pack(mem,gdb_ext,sts);

	    /*****************************************************************
	     * Create a new file with default creation and access properties.*
	     * Then create a dataset and write data to it and close the file *
	     * and dataset.                                                  *
	     *****************************************************************/

		int mpi_rank = v_cl.getProcessUnitID();
		int mpi_size = v_cl.getProcessingUnits();

		MPI_Comm comm = v_cl.getMPIComm();
		MPI_Info info  = MPI_INFO_NULL;

	    // Set up file access property list with parallel I/O access

		hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_id, comm, info);

		// Create a new file collectively and release property list identifier.
		hid_t file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		H5Pclose(plist_id);

		size_t sz = pmem.size();
		//std::cout << "Pmem.size: " << pmem.size() << std::endl;
		openfpm::vector<size_t> sz_others;
		v_cl.allGather(sz,sz_others);
		v_cl.execute();

		size_t sum = 0;

		for (size_t i = 0; i < sz_others.size(); i++)
			sum += sz_others.get(i);

		//Size for data space in file
		hsize_t fdim[1] = {sum};

		//Size for data space in file
		hsize_t fdim2[1] = {(size_t)mpi_size};

		//Create data space in file
		hid_t file_dataspace_id = H5Screate_simple(1, fdim, NULL);

		//Create data space in file
		hid_t file_dataspace_id_2 = H5Screate_simple(1, fdim2, NULL);

		//Size for data space in memory
		hsize_t mdim[1] = {pmem.size()};

		//Create data space in memory
		hid_t mem_dataspace_id = H5Screate_simple(1, mdim, NULL);

		//if (mpi_rank == 0)
			//std::cout << "Total object size: " << sum << std::endl;

		//Create data set in file
		hid_t file_dataset = H5Dcreate (file, "grid_dist", H5T_NATIVE_CHAR, file_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		//Create data set 2 in file
		hid_t file_dataset_2 = H5Dcreate (file, "metadata", H5T_NATIVE_INT, file_dataspace_id_2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	    //H5Pclose(plist_id);
	    H5Sclose(file_dataspace_id);
	    H5Sclose(file_dataspace_id_2);

	    hsize_t block[1] = {pmem.size()};

	    //hsize_t stride[1] = {1};

		hsize_t count[1] = {1};

	    hsize_t offset[1] = {0};

	    for (int i = 0; i < mpi_rank; i++)
	    {
	    	if (mpi_rank == 0)
				offset[0] = 0;
	    	else
	    		offset[0] += sz_others.get(i);
	    }

	//    std::cout << "MPI rank: " << mpi_rank << ", MPI size: " << mpi_size << ", Offset: " << offset[0] << ", Block: " << block[0] << std::endl;

	    int metadata[mpi_size];

	    for (int i = 0; i < mpi_size; i++)
	    	metadata[i] = sz_others.get(i);

	    //Select hyperslab in the file.
	    file_dataspace_id = H5Dget_space(file_dataset);
	    H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, block);

	    file_dataspace_id_2 = H5Dget_space(file_dataset_2);


	    //Create property list for collective dataset write.
	    plist_id = H5Pcreate(H5P_DATASET_XFER);
	    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

		//Write a data set to a file
		H5Dwrite(file_dataset, H5T_NATIVE_CHAR, mem_dataspace_id, file_dataspace_id, plist_id, (const char *)pmem.getPointer());

		//Write a data set 2 to a file
		H5Dwrite(file_dataset_2, H5T_NATIVE_INT, H5S_ALL, file_dataspace_id_2, plist_id, metadata);
/*
		for (size_t i = 0; i < gdb_ext.size(); i++)
		{
			Box<dim,long int> box = gdb_ext.get(i).Dbox;
			std::cout << "Dboxes saved: (" << box.getLow(0) << "; " << box.getLow(1) << "); (" << box.getHigh(0) << "; " << box.getHigh(1) << ")" << std::endl;
		}

		for (size_t i = 0; i < loc_grid.size(); i++)
		{
			std::cout << "loc_grids saved: (" << loc_grid.get(i).getGrid().getBox().getLow(0) << "; " << loc_grid.get(i).getGrid().getBox().getLow(1) << "); (" << loc_grid.get(i).getGrid().getBox().getHigh(0) << "; " << loc_grid.get(i).getGrid().getBox().getHigh(1) << ")" << std::endl;
		}
*/
	    //Close/release resources.
	    H5Dclose(file_dataset);
	    H5Sclose(file_dataspace_id);
	    H5Dclose(file_dataset_2);
	    H5Sclose(file_dataspace_id_2);
	    H5Sclose(mem_dataspace_id);
	    H5Pclose(plist_id);
	    H5Fclose(file);
	}

};

#endif /* OPENFPM_IO_SRC_HDF5_WR_HDF5_WRITER_GD_HPP_ */
