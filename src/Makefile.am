
LINKLIBS = $(PTHREAD_LIBS) $(OPT_LIBS)  $(BOOST_PROGRAM_OPTIONS_LIB) $(CUDA_LIBS) $(BOOST_IOSTREAMS_LIB)

if BUILDCUDA
 CUDA_SOURCES=../../openfpm_devices/src/memory/CudaMemory.cu
else
 CUDA_SOURCES=
endif

noinst_PROGRAMS = mem_map
mem_map_SOURCES = main.cpp Memleak_check.cpp $(CUDA_SOURCES) ../../openfpm_devices/src/memory/HeapMemory.cpp ../../openfpm_devices/src/memory/PtrMemory.cpp
mem_map_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH) $(BOOST_CPPFLAGS) -I/usr/local/include
mem_map_CFLAGS = $(CUDA_CFLAGS)
mem_map_LDADD = $(LINKLIBS) $(CUDA_LIBS) -L/usr/local/lib -L/usr/local/lib64 

nobase_include_HEADERS= data_type/scalar.hpp \
Graph/graph_unit_tests.hpp Graph/map_graph.hpp \
Grid/comb.hpp Grid/Encap.hpp Grid/grid_key.hpp Grid/grid_key_dx_expression_unit_tests.hpp Grid/grid_key_expression.hpp Grid/grid_sm.hpp Grid/grid_unit_tests.hpp Grid/grid_util_test.hpp Grid/map_grid.hpp Grid/se_grid.hpp Grid/util.hpp \
Point_test.hpp \
Point_orig.hpp \
memory_ly/memory_array.hpp memory_ly/memory_c.hpp memory_ly/memory_conf.hpp memory_ly/t_to_memory_c.hpp \
NN/CellList/CellDecomposer.hpp NN/CellList/CellKey.hpp NN/CellList/CellList.hpp NN/CellList/CellList_test.hpp NN/CellList/CellListBal.hpp NN/CellList/CellListFast.hpp NN/CellList/CellListMem.hpp NN/CellList/CellNNIterator.hpp \
Space/Ghost.hpp Space/Matrix.hpp Space/SpaceBox.hpp Space/SpaceBox_unit_tests.hpp \
Space/Shape/AdaptiveCylinderCone.hpp Space/Shape/Box.hpp Space/Shape/Box_unit_tests.hpp Space/Shape/HyperCube.hpp Space/Shape/HyperCube_unit_test.hpp Space/Shape/Point.hpp Space/Shape/Point_unit_test.hpp Space/Shape/Sphere.hpp \
util/check_no_pointers.hpp util/common.hpp util/convert.hpp util/create_vmpl_sequence.hpp util/ct_array.hpp util/for_each_ref.hpp util/mathutil.hpp util/meta_copy.hpp util/object_creator.hpp util/object_s_di.hpp util/object_si_d.hpp util/object_util.hpp util/util_debug.hpp util/util_test.hpp util/variadic_to_vmpl.hpp util/variadic_to_vmpl_unit_test.hpp \
Vector/map_vector.hpp Vector/map_vector_std.hpp Vector/util.hpp Vector/vect_isel.hpp Vector/vector_test_util.hpp Vector/vector_unit_tests.hpp \
Memleak_check.hpp \
timer.hpp

.cu.o :
	$(NVCC) $(NVCCFLAGS) $(INCLUDES_PATH) -o $@ -c $<
