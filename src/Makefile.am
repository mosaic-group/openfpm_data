
LINKLIBS = $(PTHREAD_LIBS) $(HDF5_LDFLAGS) $(HDF5_LIBS)  $(OPT_LIBS) $(BOOST_IOSTREAMS_LIB)  $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB)  $(BOOST_THREAD_LIB) 

noinst_PROGRAMS = io
io_SOURCES = main.cpp ../../openfpm_vcluster/src/VCluster/VCluster.cpp ../../openfpm_devices/src/memory/HeapMemory.cpp ../../openfpm_devices/src/Memleak_check.cpp ../../openfpm_devices/src/memory/PtrMemory.cpp
io_CXXFLAGS = $(AM_CXXFLAGS) $(CUDA_CFLAGS) $(INCLUDES_PATH) $(HDF5_CPPFLAGS) $(BOOST_CPPFLAGS)  -I/usr/local/include
io_CFLAGS = $(CUDA_CFLAGS)
io_LDADD = $(LINKLIBS) 

nobase_include_HEADERS = RawReader/RawReader.hpp RawReader/RawReader_unit_tests.hpp CSVWriter/csv_multiarray.hpp CSVWriter/CSVWriter.hpp CSVWriter/is_csv_writable.hpp \
GraphMLWriter/GraphMLWriter.hpp util/util.hpp util/GBoxes.hpp \
VTKWriter/VTKWriter.hpp VTKWriter/byteswap_portable.hpp VTKWriter/VTKWriter_dist_graph.hpp VTKWriter/VTKWriter_graph.hpp VTKWriter/VTKWriter_point_set.hpp  VTKWriter/VTKWriter_grids.hpp VTKWriter/VTKWriter_grids_st.hpp VTKWriter/VTKWriter_grids_util.hpp VTKWriter/VTKWriter_vector_box.hpp VTKWriter/is_vtk_writable.hpp HDF5_wr/HDF5_wr.hpp HDF5_wr/HDF5_writer.hpp HDF5_wr/HDF5_writer_vd.hpp HDF5_wr/HDF5_writer_gd.hpp HDF5_wr/HDF5_reader_gd.hpp  HDF5_wr/HDF5_reader.hpp HDF5_wr/HDF5_reader_vd.hpp  \
Plot/GoogleChart.hpp Plot/util.hpp 


.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

test: io
	source $(HOME)/openfpm_vars &&  cd .. && ./src/io

