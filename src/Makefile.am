
LINKLIBS = $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(CUDA_LIBS) $(BOOST_THREAD_LIB)

if BUILDCUDA
 CUDA_SOURCES=../../OpenFPM_devices/src/memory/CudaMemory.cu
else
 CUDA_SOURCES=
endif

bin_PROGRAMS = mem_map
mem_map_SOURCES = main.cpp Memleak_check.cpp $(CUDA_SOURCES) ../../OpenFPM_devices/src/memory/HeapMemory.cpp ../../OpenFPM_devices/src/memory/PtrMemory.cpp
mem_map_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH) -I/usr/local/include
mem_map_CFLAGS = $(CUDA_CFLAGS)
mem_map_LDADD = $(LINKLIBS) -L/usr/lib64/nvidia-bumblebee/ -L/usr/local/lib -L/usr/local/lib64 -lboost_unit_test_framework 

.cu.o :
	$(NVCC) $(NVCCFLAGS) $(INCLUDES_PATH) -o $@ -c $<
