
LINKLIBS = $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_IOSTREAMS_LIB)  $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_THREAD_LIB)

noinst_PROGRAMS = io
io_SOURCES = main.cpp ../../openfpm_devices/src/memory/HeapMemory.cpp ../../openfpm_devices/src/Memleak_check.cpp
io_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH) $(BOOST_CPPFLAGS) -I/usr/local/include
io_CFLAGS = $(CUDA_CFLAGS)
io_LDADD = $(LINKLIBS)

nobase_include_HEADERS = CSVWriter/csv_multiarray.hpp CSVWriter/CSVWriter.hpp GraphMLWriter/GraphMLWriter.hpp util.hpp VTKWriter/VTKWriter.hpp VTKWriter/VTKWriter_graph.hpp VTKWriter/VTKWriter_grids.hpp VTKWriter/VTKWriter_grids_st.hpp VTKWriter/VTKWriter_grids_util.hpp VTKWriter/VTKWriter_vector_box.hpp


.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
