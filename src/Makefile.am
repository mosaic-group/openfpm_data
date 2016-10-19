
LINKLIBS = $(LIBHILBERT_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(CUDA_LIBS) $(BOOST_IOSTREAMS_LIB) $(LOCAL_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) $(BOOST_CHRONO_LIB) $(BOOST_TIMER_LIB) $(BOOST_SYSTEM_LIB)

if BUILDCUDA
 CUDA_SOURCES=../../openfpm_devices/src/memory/CudaMemory.cu
else
 CUDA_SOURCES=
endif

noinst_PROGRAMS = mem_map
mem_map_SOURCES = ../../openfpm_devices/src/Memleak_check.cpp main.cpp $(CUDA_SOURCES) ../../openfpm_devices/src/memory/HeapMemory.cpp ../../openfpm_devices/src/memory/PtrMemory.cpp
mem_map_CXXFLAGS = $(AM_CXXFLAGS) $(LIBHILBERT_INCLUDE) $(CUDA_CFLAGS) $(INCLUDES_PATH) $(BOOST_CPPFLAGS) -I/usr/local/include -I/usr/local/libhilbert/include
mem_map_CFLAGS = $(CUDA_CFLAGS)
mem_map_LDADD = $(LINKLIBS)

nobase_include_HEADERS= data_type/scalar.hpp data_type/aggregate.hpp \
Graph/graph_unit_tests.hpp Graph/map_graph.hpp \
Grid/comb.hpp Grid/grid_base_implementation.hpp Grid/grid_pack_unpack.ipp  Grid/grid_base_impl_layout.hpp Grid/grid_common.hpp Grid/grid_gpu.hpp Grid/Encap.hpp Grid/grid_key.hpp Grid/grid_key_dx_expression_unit_tests.hpp Grid/grid_key_expression.hpp Grid/grid_sm.hpp Grid/grid_unit_tests.hpp Grid/grid_util_test.hpp Grid/map_grid.hpp Grid/se_grid.hpp Grid/util.hpp \
Grid/iterators/grid_key_dx_iterator_sp.hpp Grid/grid_key_dx_iterator_hilbert.hpp  Grid/iterators/grid_key_dx_iterator_sub_bc.hpp Grid/iterators/grid_key_dx_iterator_sub.hpp Grid/iterators/grid_key_dx_iterator.hpp Grid/iterators/grid_skin_iterator.hpp  \
Point_test.hpp \
Point_orig.hpp \
memory_ly/memory_array.hpp memory_ly/memory_c.hpp memory_ly/memory_conf.hpp memory_ly/t_to_memory_c.hpp \
NN/CellList/CellDecomposer.hpp NN/VerletList/VerletList.hpp NN/VerletList/VerletListFast.hpp NN/VerletList/VerletNNIterator.hpp NN/CellList/CellListNNIteratorRadius.hpp NN/CellList/CellListIterator.hpp NN/CellList/CellListM.hpp NN/CellList/CellNNIteratorM.hpp  NN/CellList/CellListFast_gen.hpp NN/CellList/CellList.hpp NN/CellList/CellList_test.hpp NN/CellList/CellListBal.hpp NN/CellList/CellListFast.hpp NN/CellList/CellListMem.hpp NN/CellList/CellNNIterator.hpp \
Space/Ghost.hpp Space/Matrix.hpp Space/SpaceBox.hpp Space/SpaceBox_unit_tests.hpp \
Space/Shape/AdaptiveCylinderCone.hpp Space/Shape/Box.hpp Space/Shape/Box_unit_tests.hpp Space/Shape/HyperCube.hpp Space/Shape/HyperCube_unit_test.hpp Space/Shape/Point.hpp Space/Shape/Point_unit_test.hpp Space/Shape/Point_operators_functions.hpp Space/Shape/Point_operators.hpp  Space/Shape/Sphere.hpp \
util/check_no_pointers.hpp util/common.hpp util/convert.hpp util/create_vmpl_sequence.hpp util/ct_array.hpp util/for_each_ref.hpp util/mathutil.hpp util/object_creator.hpp util/object_s_di.hpp util/object_si_d.hpp util/object_util.hpp util/util_debug.hpp util/util_test.hpp util/variadic_to_vmpl.hpp util/variadic_to_vmpl_unit_test.hpp util/Pack_stat.hpp \
Vector/map_vector.hpp Vector/map_vector_std_util.hpp Vector/map_vector_std_ptr.hpp Vector/map_vector_std.hpp Vector/util.hpp Vector/vect_isel.hpp Vector/vector_test_util.hpp Vector/vector_unit_tests.hpp Vector/se_vector.hpp Vector/map_vector_grow_p.hpp Vector/vector_std_pack_unpack.ipp Vector/vector_pack_unpack.ipp  \
timer.hpp \
util/copy_compare/compare_fusion_vector.hpp util/SimpleRNG.hpp util/copy_compare/compare_general.hpp util/copy_compare/copy_compare_aggregates.hpp util/copy_compare/copy_fusion_vector.hpp util/copy_compare/copy_general.hpp util/copy_compare/meta_compare.hpp util/copy_compare/meta_copy.hpp util/mul_array_extents.hpp \
Packer_Unpacker/Pack_selector.hpp Packer_Unpacker/Packer_nested_tests.hpp Packer_Unpacker/Packer_unit_tests.hpp Packer_Unpacker/Packer.hpp Packer_Unpacker/Unpacker.hpp Packer_Unpacker/Packer_util.hpp Packer_Unpacker/prp_all_zero.hpp Packer_Unpacker/has_pack_encap.hpp Packer_Unpacker/has_pack_agg.hpp Packer_Unpacker/has_max_prop.hpp

.cu.o :
	$(NVCC) $(NVCCFLAGS) $(INCLUDES_PATH) -o $@ -c $<

test: mem_map
	source $(HOME)/openfpm_vars && cd .. && ./src/mem_map

.PHONY: test

