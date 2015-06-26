
LINKLIBS = $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_IOSTREAMS_LIB)  $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_THREAD_LIB)

bin_PROGRAMS = io
io_SOURCES = main.cpp ../../OpenFPM_devices/src/memory/HeapMemory.cpp ../../OpenFPM_data/src/Memleak_check.cpp
io_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH) -fopenmp
io_CFLAGS = $(CUDA_CFLAGS)
io_LDADD = $(LINKLIBS)

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
