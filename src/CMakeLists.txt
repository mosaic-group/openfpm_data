cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

########################### Executables

add_executable(io main.cpp HDF5_wr/HDF5_writer_cuda.cu ObjReader/ObjReader_unit_test.cpp  ../../openfpm_vcluster/src/VCluster/VCluster.cpp ../../openfpm_devices/src/memory/CudaMemory.cu ../../openfpm_devices/src/memory/HeapMemory.cpp ../../openfpm_devices/src/Memleak_check.cpp ../../openfpm_devices/src/memory/PtrMemory.cpp)


###########################

if (CUDA_FOUND)
	target_compile_options(io PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe "--display_error_number --diag_suppress=2885 --diag_suppress=2887  --diag_suppress=2888 --diag_suppress=186 --diag_suppress=111" --expt-extended-lambda>)
endif()

target_include_directories (io PUBLIC ${CUDA_INCLUDE_DIRS})
target_include_directories (io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_devices/src/)
target_include_directories (io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_vcluster/src/)
target_include_directories (io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_data/src/)
target_include_directories (io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../src/)
target_include_directories (io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (io PUBLIC ${HDF5_ROOT}/include)
target_include_directories (io PUBLIC ${TINYOBJLOADER_INCLUDE_DIRS} )

target_link_libraries(io ${Boost_LIBRARIES})
target_link_libraries(io -L${HDF5_ROOT}/lib hdf5 hdf5_hl)
target_link_libraries(io -L${TINYOBJLOADER_LIBRARY_DIRS} -ltinyobjloader )

# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(io PUBLIC cxx_std_11)
target_link_libraries(io ${MPI_C_LIBRARIES})

install(FILES RawReader/RawReader.hpp 
	      RawReader/RawReader_unit_tests.hpp 
	      DESTINATION openfpm_io/include/RawReader )

install(FILES CSVWriter/csv_multiarray.hpp 
	      CSVWriter/CSVWriter.hpp 
	      CSVWriter/is_csv_writable.hpp
	      DESTINATION openfpm_io/include/CSVWriter/ )
      
install(FILES GraphMLWriter/GraphMLWriter.hpp 
	DESTINATION openfpm_io/include/GraphMLWriter )

install(FILES util/util.hpp util/GBoxes.hpp
	DESTINATION openfpm_io/include/util )

install(FILES VTKWriter/VTKWriter.hpp 
	      VTKWriter/byteswap_portable.hpp 
	      VTKWriter/VTKWriter_dist_graph.hpp 
	      VTKWriter/VTKWriter_graph.hpp 
	      VTKWriter/VTKWriter_point_set.hpp  
	      VTKWriter/VTKWriter_grids.hpp 
	      VTKWriter/VTKWriter_grids_st.hpp 
	      VTKWriter/VTKWriter_grids_util.hpp 
	      VTKWriter/VTKWriter_vector_box.hpp 
	      VTKWriter/is_vtk_writable.hpp
	DESTINATION openfpm_io/include/VTKWriter/ )

install(FILES HDF5_wr/HDF5_wr.hpp
	HDF5_wr/HDF5_writer.hpp 
	HDF5_wr/HDF5_writer_vd.hpp 
	HDF5_wr/HDF5_writer_gd.hpp 
	HDF5_wr/HDF5_reader_gd.hpp  
	HDF5_wr/HDF5_reader.hpp 
	HDF5_wr/HDF5_reader_vd.hpp
	DESTINATION openfpm_io/include/HDF5_wr)

install(FILES Plot/GoogleChart.hpp Plot/util.hpp
	DESTINATION openfpm_io/include/Plot )


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

