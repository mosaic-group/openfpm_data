cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

########################### Executables

if (CUDA_FOUND)
	set(CUDA_SOURCES Vector/vector_gpu_unit_tests.cu Grid/cuda/cuda_grid_gpu_tests.cu  Vector/cuda/map_vector_cuda_funcs_tests.cu  ../../openfpm_devices/src/memory/CudaMemory.cu NN/CellList/CellList_gpu_test.cu util/cuda/scan_cuda_unit_tests.cu  Grid/cuda/cuda_grid_unit_tests_func.cu  util/cuda/modern_gpu_tests.cu)
else()
 	set(CUDA_SOURCES )
endif()

add_executable(mem_map ../../openfpm_devices/src/Memleak_check.cpp main.cpp util/multi_array_openfpm/multi_array_ref_openfpm_unit_test.cpp memory_ly/memory_conf_unit_tests.cpp  Space/tests/SpaceBox_unit_tests.cpp  Grid/copy_grid_unit_test.cpp NN/Mem_type/Mem_type_unit_tests.cpp  ${CUDA_SOURCES} ../../openfpm_devices/src/memory/HeapMemory.cpp ../../openfpm_devices/src/memory/PtrMemory.cpp)

if ( CMAKE_COMPILER_IS_GNUCC )
    target_compile_options(mem_map PRIVATE "-Wno-deprecated-declarations")
endif()

###########################

if (CUDA_FOUND)
	target_compile_options(mem_map PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe "--display_error_number --diag_suppress=611 --diag_suppress=2885 --diag_suppress=2886  --diag_suppress=2887  --diag_suppress=2888 --diag_suppress=186 --diag_suppress=111" --expt-extended-lambda>)
endif()

target_include_directories (mem_map PUBLIC ${CUDA_INCLUDE_DIRS})
target_include_directories (mem_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (mem_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_devices/src/)
target_include_directories (mem_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_io/src/)
target_include_directories (mem_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (mem_map PUBLIC ${LIBHILBERT_INCLUDE_DIRS})
target_include_directories (mem_map PUBLIC ${Boost_INCLUDE_DIRS})

target_link_libraries(mem_map ${Boost_LIBRARIES})
target_link_libraries(mem_map -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})


# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(mem_map PUBLIC cxx_std_11)
target_link_libraries(mem_map ${MPI_C_LIBRARIES})

install(FILES Grid/comb.hpp 
	      Grid/copy_grid_fast.hpp 
	      Grid/grid_base_implementation.hpp 
	      Grid/grid_pack_unpack.ipp  
	      Grid/grid_base_impl_layout.hpp 
	      Grid/grid_common.hpp 
	      Grid/grid_gpu.hpp Grid/Encap.hpp 
	      Grid/grid_key.hpp Grid/grid_key_dx_expression_unit_tests.hpp 
	      Grid/grid_key_expression.hpp Grid/grid_sm.hpp 
	      Grid/grid_unit_tests.hpp Grid/grid_util_test.hpp 
	      Grid/map_grid.hpp Grid/se_grid.hpp Grid/util.hpp
	      Grid/iterators/grid_key_dx_iterator_sp.hpp 
	      Grid/grid_key_dx_iterator_hilbert.hpp
	      Grid/grid_key_dx_iterator_hilbert.hpp 
              DESTINATION openfpm_data/include/Grid/ )

install(FILES Grid/cuda/cuda_grid_gpu_funcs.cuh
	      Grid/cuda/map_grid_cuda_ker.cuh
	      DESTINATION openfpm_data/include/Grid/cuda/ )

install(FILES data_type/aggregate.hpp
	DESTINATION openfpm_data/include/data_type )

install(FILES Graph/map_graph.hpp
	DESTINATION openfpm_data/include/Graph)

install(FILES Point_test.hpp
	      Point_orig.hpp
	      timer.hpp
	DESTINATION openfpm_data/include/ )

install(FILES memory_ly/memory_array.hpp 
	      memory_ly/memory_c.hpp 
	      memory_ly/memory_conf.hpp 
	      memory_ly/t_to_memory_c.hpp
	DESTINATION openfpm_data/include/memory_ly )

install(FILES NN/Mem_type/MemBalanced.hpp
              NN/Mem_type/MemFast.hpp 
              NN/Mem_type/MemMemoryWise.hpp
	      DESTINATION openfpm_data/include/NN/Mem_type )

install(FILES NN/CellList/CellListNNIteratorRadius.hpp 
	      NN/CellList/CellListIterator.hpp 
	      NN/CellList/CellListM.hpp 
	      NN/CellList/CellNNIteratorM.hpp 
	      NN/CellList/CellList.hpp 
	      NN/CellList/CellList_test.hpp 
	      NN/CellList/CellListFast_gen.hpp
	      NN/CellList/CellList_util.hpp
	      NN/CellList/CellNNIterator.hpp 
	      NN/CellList/ProcKeys.hpp
              NN/CellList/CellNNIteratorRuntime.hpp 
	      NN/CellList/NNc_array.hpp 
	      NN/CellList/ParticleItCRS_Cells.hpp 
	      NN/CellList/ParticleIt_Cells.hpp 
	      NN/CellList/CellDecomposer.hpp
	      NN/CellList/CellNNIteratorRuntime.hpp
	      NN/CellList/CellNNIteratorRuntimeM.hpp
	      DESTINATION openfpm_data/include/NN/CellList )

      install(FILES NN/CellList/cuda/CellList_cpu_ker.cuh
		    NN/CellList/cuda/CellList_gpu_ker.cuh
		    NN/CellList/cuda/CellList_gpu.hpp
		    NN/CellList/cuda/Cuda_cell_list_util_func.hpp
	            DESTINATION openfpm_data/include/NN/CellList/cuda )


install(FILES NN/VerletList/VerletList.hpp 
	      NN/VerletList/VerletListFast.hpp 
	      NN/VerletList/VerletNNIterator.hpp
	      NN/VerletList/VerletListM.hpp 
	      NN/VerletList/VerletNNIteratorM.hpp
	DESTINATION openfpm_data/include/NN/VerletList/)

install(FILES Space/Ghost.hpp Space/Matrix.hpp Space/SpaceBox.hpp
	DESTINATION openfpm_data/include/Space/ )

install(FILES Space/Shape/AdaptiveCylinderCone.hpp 
	      Space/Shape/Box.hpp 
	      Space/Shape/Box_unit_tests.hpp 
	      Space/Shape/HyperCube.hpp 
	      Space/Shape/HyperCube_unit_test.hpp 
	      Space/Shape/Point.hpp 
	      Space/Shape/Point_unit_test.hpp 
	      Space/Shape/Point_operators_functions.hpp 
	      Space/Shape/Point_operators.hpp  
	      Space/Shape/Sphere.hpp
	      DESTINATION openfpm_data/include/Space/Shape )



install(FILES util/check_no_pointers.hpp 
	      util/common.hpp util/convert.hpp 
	      util/create_vmpl_sequence.hpp 
	      util/ct_array.hpp 
	      util/for_each_ref.hpp 
	      util/mathutil.hpp util/object_creator.hpp 
	      util/object_s_di.hpp 
	      util/object_si_d.hpp 
	      util/object_util.hpp 
	      util/util_debug.hpp 
	      util/variadic_to_vmpl.hpp 
	      util/Pack_stat.hpp
	      util/cuda/cuda_launch.hpp
	      util/tokernel_transformation.hpp
	      util/SimpleRNG.hpp
	      util/math_util_complex.hpp
	      util/mul_array_extents.hpp
	      DESTINATION openfpm_data/include/util )

install(FILES util/copy_compare/compare_fusion_vector.hpp 
	      util/copy_compare/compare_general.hpp 
	      util/copy_compare/copy_compare_aggregates.hpp 
	      util/copy_compare/copy_fusion_vector.hpp 
	      util/copy_compare/copy_general.hpp util/copy_compare/meta_compare.hpp 
	      util/copy_compare/meta_copy.hpp
	      DESTINATION openfpm_data/include/util/copy_compare/ )

install(FILES Packer_Unpacker/Pack_selector.hpp
	      Packer_Unpacker/Packer_nested_tests.hpp 
	      Packer_Unpacker/Packer_unit_tests.hpp 
	      Packer_Unpacker/Packer.hpp 
	      Packer_Unpacker/Unpacker.hpp 
	      Packer_Unpacker/Packer_util.hpp 
	      Packer_Unpacker/prp_all_zero.hpp 
	      Packer_Unpacker/has_pack_encap.hpp 
	      Packer_Unpacker/has_pack_agg.hpp 
	      Packer_Unpacker/has_max_prop.hpp
	      DESTINATION openfpm_data/include/Packer_Unpacker )

install(FILES util/boost/boost_multi_array_base_openfpm.hpp
	      util/boost/boost_multi_array_openfpm.hpp
	      util/boost/boost_multi_array_subarray_openfpm.hpp
	      util/boost/boost_multi_array_ref_openfpm.hpp
	      util/boost/boost_multi_array_view_openfpm.hpp
	      util/boost/boost_array_openfpm.hpp
	      util/boost/boost_multi_array_iterator_openfpm.hpp
	      DESTINATION openfpm_data/include/util/boost/ )

install(FILES Grid/iterators/grid_key_dx_iterator_sp.hpp 
	      Grid/iterators/stencil_type.hpp  
	      Grid/iterators/grid_key_dx_iterator_sub_bc.hpp 
	      Grid/iterators/grid_key_dx_iterator_sub.hpp 
	      Grid/iterators/grid_key_dx_iterator.hpp 
	      Grid/iterators/grid_skin_iterator.hpp
	      DESTINATION openfpm_data/include/Grid/iterators)

install(FILES Vector/map_vector.hpp 
	      Vector/vector_def.hpp 
	      Vector/map_vector_std_util.hpp 
	      Vector/map_vector_std_ptr.hpp 
	      Vector/map_vector_std.hpp 
	      Vector/util.hpp 
	      Vector/vect_isel.hpp 
	      Vector/vector_test_util.hpp 
	      Vector/vector_unit_tests.hpp 
	      Vector/se_vector.hpp 
	      Vector/map_vector_grow_p.hpp 
	      Vector/vector_std_pack_unpack.ipp 
	      Vector/vector_pack_unpack.ipp 
	      Vector/vector_map_iterator.hpp
	      DESTINATION openfpm_data/include/Vector )

install(FILES Vector/cuda/map_vector_cuda_ker.cuh
	      DESTINATION openfpm_data/include/Vector/cuda/ )

install(FILES util/multi_array_openfpm/array_openfpm.hpp
	      util/multi_array_openfpm/multi_array_iterator_openfpm.hpp 
	      util/multi_array_openfpm/multi_array_ref_base_openfpm.hpp
	      util/multi_array_openfpm/multi_array_ref_openfpm.hpp
	      util/multi_array_openfpm/multi_array_ref_subarray_openfpm.hpp
	      util/multi_array_openfpm/multi_array_view_openfpm.hpp
	      util/multi_array_openfpm/storage_order.hpp
	      util/multi_array_openfpm/types.hpp 
	      DESTINATION openfpm_data/include/util/multi_array_openfpm )

install(FILES util/cuda/scan_cuda.cuh
	      util/cuda/ofp_context.hxx
	      DESTINATION openfpm_data/include/util/cuda )

install(FILES util/cuda/moderngpu/context.hxx
	      util/cuda/moderngpu/cpp11.hxx
	      util/cuda/moderngpu/cta_load_balance.hxx
	      util/cuda/moderngpu/cta_merge.hxx
	      util/cuda/moderngpu/cta_mergesort.hxx
	      util/cuda/moderngpu/cta_reduce.hxx
	      util/cuda/moderngpu/cta_scan.hxx
	      util/cuda/moderngpu/cta_search.hxx
	      util/cuda/moderngpu/cta_segscan.hxx
	      util/cuda/moderngpu/intrinsics.hxx
	      util/cuda/moderngpu/kernel_bulkinsert.hxx
	      util/cuda/moderngpu/kernel_bulkremove.hxx
	      util/cuda/moderngpu/kernel_compact.hxx
	      util/cuda/moderngpu/kernel_intervalmove.hxx
	      util/cuda/moderngpu/kernel_join.hxx
	      util/cuda/moderngpu/kernel_load_balance.hxx
	      util/cuda/moderngpu/kernel_merge.hxx
	      util/cuda/moderngpu/kernel_mergesort.hxx
	      util/cuda/moderngpu/kernel_reduce.hxx
	      util/cuda/moderngpu/kernel_scan.hxx
	      util/cuda/moderngpu/kernel_segreduce.hxx
	      util/cuda/moderngpu/kernel_segsort.hxx
	      util/cuda/moderngpu/kernel_sortedsearch.hxx
	      util/cuda/moderngpu/kernel_workcreate.hxx
	      util/cuda/moderngpu/launch_box.hxx
	      util/cuda/moderngpu/launch_params.hxx
	      util/cuda/moderngpu/loadstore.hxx
	      util/cuda/moderngpu/memory.hxx
	      util/cuda/moderngpu/meta.hxx
	      util/cuda/moderngpu/operators.hxx
	      util/cuda/moderngpu/search.hxx
	      util/cuda/moderngpu/sort_networks.hxx
	      util/cuda/moderngpu/transform.hxx
	      util/cuda/moderngpu/tuple.hxx
	      util/cuda/moderngpu/types.hxx
	      util/cuda/moderngpu/util.hxx
              DESTINATION openfpm_data/include/util/cuda/moderngpu/ )

#if(BUILD_TESTING)

#  add_executable(particle_test test.cu)

#  set_target_properties(particle_test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#  target_link_libraries(particle_test PRIVATE particles)

#  add_test(NAME particles_10k COMMAND particle_test 10000 )
#  add_test(NAME particles_256k COMMAND particle_test 256000 )

#  if(APPLE)
    # We need to add the default path to the driver (libcuda.dylib) as an rpath,
    # so that the static cuda runtime can find it at runtime.
    #    set_property(TARGET particle_test PROPERTY BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    #  endif()
    #endif()

