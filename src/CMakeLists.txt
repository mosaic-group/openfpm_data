cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

########################### Executables


add_executable(io main.cpp 
				  MetaParser/MetaParser_unit_test.cpp 
				  HDF5_wr/HDF5_writer_cuda.cu 
				  ObjReader/ObjReader_unit_test.cpp)

add_dependencies(io ofpmmemory)
add_dependencies(io vcluster)

if ( CMAKE_COMPILER_IS_GNUCC )
    target_compile_options(io PRIVATE "-Wno-deprecated-declarations")
    if (TEST_COVERAGE)
        target_compile_options(io PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -fprofile-arcs -ftest-coverage>)
    endif()
endif()

if (CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
        add_definitions(-D__STRICT_ANSI__)
endif()

###########################

if (CUDA_FOUND)
	target_include_directories (io PUBLIC ${MPI_C_INCLUDE_DIRS})

   	if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0 )
    	add_compile_options( "-fabi-version=6" )
    endif()

    if (TEST_COVERAGE)
    	target_compile_options(io PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler "-fprofile-arcs -ftest-coverage" >)
	endif()
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
target_include_directories (io PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories (io PUBLIC ${OPENFPM_PDATA_DIR}/src)
target_include_directories (io PUBLIC ${ALPAKA_ROOT}/include)

if (PETSC_FOUND)
	target_include_directories(io PUBLIC ${PETSC_INCLUDES})
endif()

target_link_libraries(io ${Boost_LIBRARIES})
target_link_libraries(io ${HDF5_LIBRARIES})
target_link_libraries(io ${TINYOBJLOADER_LIBRARIES} )
target_link_libraries(io vcluster)
target_link_libraries(io ofpmmemory)
target_link_libraries(io OpenMP::OpenMP_CXX)
if (PETSC_FOUND)
	target_link_libraries(io ${PETSC_LIBRARIES})
endif()

# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(io PUBLIC cxx_std_11)
target_link_libraries(io ${MPI_C_LIBRARIES})

if (TEST_COVERAGE)
    target_link_libraries(io -lgcov)
endif()

install(FILES RawReader/RawReader.hpp 
	      RawReader/RawReader_unit_tests.hpp 
	      DESTINATION openfpm_io/include/RawReader 
	      COMPONENT OpenFPM)

install(FILES CSVWriter/csv_multiarray.hpp 
	      CSVWriter/CSVWriter.hpp 
	      CSVWriter/is_csv_writable.hpp
	      DESTINATION openfpm_io/include/CSVWriter/ 
	      COMPONENT OpenFPM)
      
install(FILES GraphMLWriter/GraphMLWriter.hpp 
	DESTINATION openfpm_io/include/GraphMLWriter 
	COMPONENT OpenFPM)

install(FILES util/util.hpp util/GBoxes.hpp
	DESTINATION openfpm_io/include/util 
	COMPONENT OpenFPM)

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
	DESTINATION openfpm_io/include/VTKWriter/ 
	COMPONENT OpenFPM)

install(FILES MetaParser/MetaParser.hpp
	DESTINATION openfpm_io/include/MetaParser/ 
	COMPONENT OpenFPM)

install(FILES HDF5_wr/HDF5_wr.hpp
	HDF5_wr/HDF5_writer.hpp 
	HDF5_wr/HDF5_writer_vd.hpp 
	HDF5_wr/HDF5_writer_gd.hpp 
	HDF5_wr/HDF5_reader_gd.hpp  
	HDF5_wr/HDF5_reader.hpp 
	HDF5_wr/HDF5_reader_vd.hpp
	DESTINATION openfpm_io/include/HDF5_wr
	COMPONENT OpenFPM)

install(FILES RawReader/InitGridWithPixel.hpp
	DESTINATION openfpm_io/include/RawReader/
	COMPONENT OpenFPM)

install(FILES util/PathsAndFiles.hpp
	DESTINATION openfpm_io/include/util
	COMPONENT OpenFPM)

install(FILES Plot/GoogleChart.hpp Plot/util.hpp
	DESTINATION openfpm_io/include/Plot 
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

