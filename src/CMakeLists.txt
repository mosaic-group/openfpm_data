cmake_minimum_required(VERSION 3.8 FATAL_ERROR)


########################### Executables

add_definitions(-DSCAN_WITH_CUB)


if (TEST_PERFORMANCE)
	set(CUDA_SOURCES SparseGridGpu/performance/SparseGridGpu_performance_heat_stencil_sparse.cu
			 SparseGridGpu/performance/SparseGridGpu_performance_insert_stencil.cu
			 SparseGridGpu/performance/SparseGridGpu_performance_heat_stencil.cu 
			 SparseGridGpu/performance/SparseGridGpu_performance_stencil_heat_host.cu  
			 SparseGridGpu/performance/SparseGridGpu_performance_get_nn.cu  
			 SparseGridGpu/performance/SparseGridGpu_performance_get_single.cu 
			 SparseGridGpu/performance/SparseGridGpu_performance_insert_single.cu
			 SparseGridGpu/performance/SparseGridGpu_performance_tests.cu
			 SparseGridGpu/performance/SparseGridGpu_performance_insert_block.cu
                         SparseGridGpu/performance/SparseGridGpu_performance_heat_stencil_3d.cu
                         SparseGridGpu/performance/performancePlots.cpp
                         Vector/performance/vector_performance_test.cu)
endif ()


if (NOT CUDA_ON_BACKEND STREQUAL "None")
	set(CUDA_SOURCES ${CUDA_SOURCES}
	    Vector/map_vector_sparse_unit_tests.cu
            Vector/vector_gpu_unit_tests.cu
            Grid/cuda/cuda_grid_gpu_tests.cu
            Vector/cuda/map_vector_cuda_funcs_tests.cu
            NN/CellList/CellList_gpu_test.cu 
	    util/cuda/scan_sort_cuda_unit_tests.cu
            Grid/cuda/cuda_grid_unit_tests_func.cu
            util/cuda/modern_gpu_tests.cu
            Vector/cuda/map_vector_sparse_cuda_ker_unit_tests.cu
            NN/CellList/tests/CellDecomposer_gpu_ker_unit_test.cu
	    SparseGridGpu/tests/BlockMapGpu_tests.cu
	    util/test/zmorton_unit_tests.cpp
	    util/cuda/test/segreduce_block_cuda_tests.cu 
	    SparseGridGpu/tests/BlockMapGpu_kernels_tests.cu 
	    SparseGridGpu/tests/SparseGridGpu_tests.cu
	    )
endif ()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	add_definitions("-DBOOST_MPL_CFG_HAS_TYPEOF")
endif()

if (CUDA_ON_BACKEND STREQUAL "OpenMP" OR CUDA_ON_BACKEND STREQUAL "SEQUENTIAL")
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES COMPILE_FLAGS "-D__NVCC__ -DCUDART_VERSION=11000")
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CXX)
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
                add_definitions("-x c++")
        endif()
endif()

if ( CUDA_ON_BACKEND STREQUAL "HIP" AND HIP_FOUND )

	list(APPEND HIP_HIPCC_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})

	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		list(APPEND HIP_HIPCC_FLAGS -O0)
	endif()

	list(APPEND HIP_HIPCC_FLAGS -D__NVCC__ -D__HIP__  -DCUDART_VERSION=11000 -D__CUDACC__ -D__CUDACC_VER_MAJOR__=11 -D__CUDACC_VER_MINOR__=0 -D__CUDACC_VER_BUILD__=0 --std=c++14)
	set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CXX)


        hip_add_executable(mem_map ${CUDA_SOURCES}
                main.cpp
                data_type/aggregate_unit_tests.cpp
                util/multi_array_openfpm/multi_array_ref_openfpm_unit_test.cpp
                memory_ly/memory_conf_unit_tests.cpp
                Space/tests/SpaceBox_unit_tests.cpp
                Space/Shape/Sphere_unit_test.cpp
                SparseGrid/SparseGrid_unit_tests.cpp
                SparseGrid/SparseGrid_chunk_copy_unit_tests.cpp
                Grid/copy_grid_unit_test.cpp NN/Mem_type/Mem_type_unit_tests.cpp
                Grid/Geometry/tests/grid_smb_tests.cpp
                Space/Shape/Plane_unit_test.cpp)


else()

	add_executable(mem_map ${CUDA_SOURCES}
        	main.cpp
        	data_type/aggregate_unit_tests.cpp
        	util/multi_array_openfpm/multi_array_ref_openfpm_unit_test.cpp
        	memory_ly/memory_conf_unit_tests.cpp
        	Space/tests/SpaceBox_unit_tests.cpp
        	Space/Shape/Sphere_unit_test.cpp
		SparseGrid/SparseGrid_unit_tests.cpp
		SparseGrid/SparseGrid_chunk_copy_unit_tests.cpp
        	Grid/copy_grid_unit_test.cpp NN/Mem_type/Mem_type_unit_tests.cpp
                Grid/Geometry/tests/grid_smb_tests.cpp
                Space/Shape/Plane_unit_test.cpp)

endif()

set_property(TARGET mem_map PROPERTY CUDA_ARCHITECTURES 60 75)

if (CUDA_FOUND)

	target_include_directories(mem_map PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

	add_executable(isolation
		isolation.cu
        	)

	add_dependencies(isolation ofpmmemory)

	set_property(TARGET isolation PROPERTY CUDA_ARCHITECTURES OFF)

endif()

add_dependencies(mem_map ofpmmemory)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	target_compile_options(mem_map PRIVATE "-Wno-undefined-var-template")
	target_compile_options(mem_map PRIVATE "-Wno-macro-redefined")
endif()

if (CMAKE_COMPILER_IS_GNUCC)
    target_compile_options(mem_map PRIVATE "-Wno-deprecated-declarations")

    if (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64")
        target_compile_options(mem_map PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -mavx>)
    endif()

   if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0 )
        target_compile_options(mem_map PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -fabi-version=6>)
    endif()

    if (TEST_COVERAGE)
        target_compile_options(mem_map PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -fprofile-arcs -ftest-coverage>)
    endif ()
endif ()

if (CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
        add_definitions(-D__STRICT_ANSI__)
endif()


###########################

if (TEST_PERFORMANCE)
    target_include_directories(mem_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_io/src/)
    target_include_directories(isolation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_io/src/)
endif ()

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_include_directories(mem_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_io/src/)
    if (CUDA_FOUND)
    	target_include_directories(isolation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_io/src/)
    endif()
ENDIF()


if (CUDA_FOUND)
    if (CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    	target_compile_options(mem_map PUBLIC $<$<COMPILE_LANGUAGE:CUDA>: ${WARNING_SUPPRESSION_AND_OPTION_NVCC} -use_fast_math -gencode arch=compute_61,code=compute_61 -g -lineinfo >)
    	target_compile_options(isolation PUBLIC $<$<COMPILE_LANGUAGE:CUDA>: ${WARNING_SUPPRESSION_AND_OPTION_NVCC} -use_fast_math -gencode arch=compute_61,code=compute_61 -g -lineinfo >)
    else()
	target_compile_options(mem_map PUBLIC $<$<COMPILE_LANGUAGE:CUDA>: -D__NVCC__ -D__CUDACC_VER_MAJOR__=${CUDA_VERSION_MAJOR} -D__CUDACC_VER_MINOR__=${CUDA_VERSION_MINOR} -D__CUDACC_VER_BUILD__=0 >)
        target_compile_options(isolation PUBLIC $<$<COMPILE_LANGUAGE:CUDA>: -D__NVCC__ -D__CUDACC_VER_MAJOR__=${CUDA_VERSION_MAJOR} -D__CUDACC_VER_MINOR__=${CUDA_VERSION_MINOR} -D__CUDACC_VER_BUILD__=0 >)
    endif()
    if (TEST_COVERAGE)
        target_compile_options(mem_map PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler "-fprofile-arcs -ftest-coverage" >)
    endif ()
    if (TEST_COVERAGE)
        target_compile_options(mem_map PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler "-fprofile-arcs -ftest-coverage" >)
    endif ()
endif ()

target_include_directories(mem_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(mem_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_devices/src/)
target_include_directories(mem_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories(mem_map PUBLIC ${LIBHILBERT_INCLUDE_DIRS})
target_include_directories(mem_map PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(mem_map PUBLIC ${ALPAKA_ROOT}/include)

if (CUDA_FOUND)
	target_include_directories(isolation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
	target_include_directories(isolation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_devices/src/)
	target_include_directories(isolation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config)
	target_include_directories(isolation PUBLIC ${LIBHILBERT_INCLUDE_DIRS})
	target_include_directories(isolation PUBLIC ${Boost_INCLUDE_DIRS})
	target_include_directories(isolation PUBLIC ${ALPAKA_ROOT}/include)
endif()

target_include_directories(mem_map PUBLIC ${Vc_INCLUDE_DIR})

if (HIP_FOUND)
	target_link_libraries(mem_map hip::host)
endif()
target_link_libraries(mem_map ${Boost_LIBRARIES})
target_link_libraries(mem_map -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(mem_map ofpmmemory)
target_link_libraries(mem_map ${Vc_LIBRARIES})
if (OPENMP_FOUND)
	target_link_libraries(mem_map OpenMP::OpenMP_CXX)
endif()

if (CUDA_FOUND)
	target_link_libraries(isolation ${Boost_LIBRARIES})
	target_link_libraries(isolation -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
	target_link_libraries(isolation ofpmmemory)
endif()

if (TEST_COVERAGE)
    target_link_libraries(mem_map -lgcov)
endif ()

# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(mem_map PUBLIC cxx_std_11)
target_link_libraries(mem_map ${MPI_C_LIBRARIES})
target_link_libraries(mem_map m)
if (NOT APPLE)
    target_link_libraries(mem_map rt)
endif ()

install(FILES Grid/comb.hpp
        Grid/copy_grid_fast.hpp
        Grid/grid_base_implementation.hpp
        Grid/grid_pack_unpack.ipp
        Grid/grid_base_impl_layout.hpp
        Grid/grid_common.hpp
        Grid/grid_gpu.hpp
        Grid/grid_key.hpp Grid/grid_key_dx_expression_unit_tests.hpp
        Grid/grid_key_expression.hpp 
	Grid/grid_sm.hpp
	Grid/grid_zm.hpp
        Grid/grid_unit_tests.hpp Grid/grid_util_test.hpp
        Grid/map_grid.hpp Grid/se_grid.hpp Grid/util.hpp
        Grid/iterators/grid_key_dx_iterator_sp.hpp
        Grid/grid_key_dx_iterator_hilbert.hpp
        Grid/grid_key_dx_iterator_hilbert.hpp
        DESTINATION openfpm_data/include/Grid/
	COMPONENT OpenFPM)

install(FILES Grid/cuda/cuda_grid_gpu_funcs.cuh
        Grid/cuda/map_grid_cuda_ker.cuh
        DESTINATION openfpm_data/include/Grid/cuda/
	COMPONENT OpenFPM)

install(FILES data_type/aggregate.hpp
        DESTINATION openfpm_data/include/data_type
	COMPONENT OpenFPM)

install(FILES Graph/CartesianGraphFactory.hpp
        Graph/map_graph.hpp
        DESTINATION openfpm_data/include/Graph
	COMPONENT OpenFPM)

install(FILES Point_test.hpp
        Point_orig.hpp
        timer.hpp
        DESTINATION openfpm_data/include/
	COMPONENT OpenFPM)

install(FILES SparseGrid/SparseGrid_iterator.hpp
              SparseGrid/SparseGrid.hpp
              SparseGrid/SparseGridUtil.hpp
	      SparseGrid/SparseGrid_iterator_block.hpp
	      SparseGrid/SparseGrid_chunk_copy.hpp
	      SparseGrid/SparseGrid_conv_opt.hpp
	      SparseGrid/cp_block.hpp
        DESTINATION openfpm_data/include/SparseGrid
	COMPONENT OpenFPM)

install(FILES hash_map/hopscotch_hash.h
              hash_map/hopscotch_map.h
              hash_map/hopscotch_sc_map.h
              hash_map/hopscotch_sc_set.h
              hash_map/hopscotch_set.h
        DESTINATION openfpm_data/include/hash_map
	COMPONENT OpenFPM)

install(FILES memory_ly/Encap.hpp
	memory_ly/memory_array.hpp
        memory_ly/memory_c.hpp
        memory_ly/memory_conf.hpp
        memory_ly/t_to_memory_c.hpp
        DESTINATION openfpm_data/include/memory_ly
	COMPONENT OpenFPM)

install(FILES NN/Mem_type/MemBalanced.hpp
        NN/Mem_type/MemFast.hpp
        NN/Mem_type/MemMemoryWise.hpp
        DESTINATION openfpm_data/include/NN/Mem_type
	COMPONENT OpenFPM)

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
	NN/CellList/CellList_def.hpp
        DESTINATION openfpm_data/include/NN/CellList
	COMPONENT OpenFPM)

install(FILES NN/CellList/cuda/CellList_cpu_ker.cuh
        NN/CellList/cuda/CellList_gpu_ker.cuh
        NN/CellList/cuda/CellList_gpu.hpp
        NN/CellList/cuda/Cuda_cell_list_util_func.hpp
	NN/CellList/cuda/CellDecomposer_gpu_ker.cuh
        DESTINATION openfpm_data/include/NN/CellList/cuda
	COMPONENT OpenFPM)


install(FILES NN/VerletList/VerletList.hpp
        NN/VerletList/VerletListFast.hpp
        NN/VerletList/VerletNNIterator.hpp
        NN/VerletList/VerletListM.hpp
        NN/VerletList/VerletNNIteratorM.hpp
        DESTINATION openfpm_data/include/NN/VerletList/
	COMPONENT OpenFPM)

install(FILES Space/Ghost.hpp Space/Matrix.hpp Space/SpaceBox.hpp
        DESTINATION openfpm_data/include/Space/
	COMPONENT OpenFPM)

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
        DESTINATION openfpm_data/include/Space/Shape
	COMPONENT OpenFPM)


install(FILES util/check_no_pointers.hpp
        util/common.hpp util/convert.hpp
        util/create_vmpl_sequence.hpp
        util/ct_array.hpp
        util/for_each_ref.hpp
	util/for_each_ref_host.hpp
        util/mathutil.hpp util/object_creator.hpp
	util/object_si_di.hpp
        util/object_s_di.hpp
	util/zmorton.hpp
        util/object_si_d.hpp
        util/object_util.hpp
        util/util_debug.hpp
        util/variadic_to_vmpl.hpp
        util/Pack_stat.hpp
        util/tokernel_transformation.hpp
        util/SimpleRNG.hpp
        util/math_util_complex.hpp
        util/mul_array_extents.hpp
	util/hostDevice_util_funcs.hpp
	util/sparsegrid_util_common.hpp
        DESTINATION openfpm_data/include/util
	COMPONENT OpenFPM)

install(FILES util/copy_compare/compare_fusion_vector.hpp
        util/copy_compare/compare_general.hpp
        util/copy_compare/copy_compare_aggregates.hpp
        util/copy_compare/copy_fusion_vector.hpp
        util/copy_compare/copy_general.hpp util/copy_compare/meta_compare.hpp
        util/copy_compare/meta_copy.hpp
        DESTINATION openfpm_data/include/util/copy_compare/
	COMPONENT OpenFPM)

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
        DESTINATION openfpm_data/include/Packer_Unpacker
	COMPONENT OpenFPM)

install(FILES util/boost/boost_array_openfpm.hpp
        DESTINATION openfpm_data/include/util/boost/
	COMPONENT OpenFPM)

install(FILES Grid/iterators/grid_key_dx_iterator_sp.hpp
        Grid/iterators/stencil_type.hpp
        Grid/iterators/grid_key_dx_iterator_sub_bc.hpp
        Grid/iterators/grid_key_dx_iterator_sub.hpp
        Grid/iterators/grid_key_dx_iterator.hpp
        Grid/iterators/grid_skin_iterator.hpp
        DESTINATION openfpm_data/include/Grid/iterators
	COMPONENT OpenFPM)

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
        Vector/map_vector_printers.hpp
        Vector/map_vector_sparse.hpp
        DESTINATION openfpm_data/include/Vector
	COMPONENT OpenFPM)

install(FILES Vector/cuda/map_vector_cuda_ker.cuh
        Vector/cuda/map_vector_std_cuda_ker.cuh
        Vector/cuda/map_vector_std_cuda.hpp
        Vector/cuda/map_vector_sparse_cuda_ker.cuh
        Vector/cuda/map_vector_sparse_cuda_kernels.cuh
        DESTINATION openfpm_data/include/Vector/cuda/
	COMPONENT OpenFPM)

install(DIRECTORY util/cuda/cub_old/
        DESTINATION openfpm_data/include/util/cuda/cub_old
	COMPONENT OpenFPM)

install(FILES util/multi_array_openfpm/array_openfpm.hpp
        util/multi_array_openfpm/multi_array_iterator_openfpm.hpp
        util/multi_array_openfpm/multi_array_ref_base_openfpm.hpp
        util/multi_array_openfpm/multi_array_ref_openfpm.hpp
        util/multi_array_openfpm/multi_array_ref_subarray_openfpm.hpp
        util/multi_array_openfpm/multi_array_view_openfpm.hpp
        util/multi_array_openfpm/storage_order.hpp
        util/multi_array_openfpm/types.hpp
        DESTINATION openfpm_data/include/util/multi_array_openfpm
	COMPONENT OpenFPM)


install(FILES util/cuda/ofp_context.hxx
        util/cuda/kernels.cuh
        util/cuda/scan_ofp.cuh
	util/cuda/sort_ofp.cuh
	util/cuda/reduce_ofp.cuh
	util/cuda/segreduce_ofp.cuh
	util/cuda/merge_ofp.cuh
        DESTINATION openfpm_data/include/util/cuda
	COMPONENT OpenFPM)

install(FILES util/stat/common_statistics.hpp
        DESTINATION openfpm_data/include/util/stat
	COMPONENT OpenFPM)

install(FILES SparseGridGpu/SparseGridGpu.hpp
	      SparseGridGpu/SparseGridGpu_kernels.cuh
	      SparseGridGpu/SparseGridGpu_ker.cuh
	      SparseGridGpu/SparseGridGpu_ker_util.hpp
	      SparseGridGpu/DataBlock.cuh
	      SparseGridGpu/BlockMapGpu.hpp
	      SparseGridGpu/BlockMapGpu_kernels.cuh
	      SparseGridGpu/BlockMapGpu_ker.cuh
	      SparseGridGpu/BlockMapGpu_dimensionalityWrappers.cuh
	      SparseGridGpu/encap_num.hpp
              DESTINATION openfpm_data/include/SparseGridGpu
	      COMPONENT OpenFPM)

install (FILES Grid/Geometry/grid_smb.hpp
	      Grid/Geometry/grid_zmb.hpp
              DESTINATION openfpm_data/include/Grid/Geometry/
	      COMPONENT OpenFPM)

install (FILES SparseGridGpu/Iterators/SparseGridGpu_iterator_sub.hpp
	       SparseGridGpu/Iterators/SparseGridGpu_iterator.hpp
               DESTINATION openfpm_data/include/SparseGridGpu/Iterators/ 
	       COMPONENT OpenFPM)

install (FILES  SparseGridGpu/TemplateUtils/mathUtils.hpp
	        SparseGridGpu/TemplateUtils/encap_shmem.hpp
		DESTINATION openfpm_data/include/SparseGridGpu/TemplateUtils/ 
		COMPONENT OpenFPM)

install(FILES util/cuda/moderngpu/context.hxx
        util/cuda/moderngpu/context_reduced.hxx
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
        DESTINATION openfpm_data/include/util/cuda/moderngpu/
	COMPONENT OpenFPM)

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

