
LINKLIBS = $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(CUDA_LIBS) $(BOOST_THREAD_LIB)

bin_PROGRAMS = mem_map
mem_map_SOURCES = main.cpp
mem_map_CXXFLAGS = $(CUDA_CFLAGS)
mem_map_CFLAGS = $(CUDA_CFLAGS)
mem_map_LDADD = $(LINKLIBS) -L/usr/lib64/nvidia-bumblebee/

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
