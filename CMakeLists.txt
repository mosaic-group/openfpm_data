cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

if (POLICY CMP0074)
        cmake_policy(SET CMP0074 OLD)
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake_modules/)

if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

set(BOOST_INCLUDE ${Boost_INCLUDE_DIR} CACHE PATH "Include directory for BOOST")
set(HDF5_ROOT CACHE PATH "HDF5 root path")
set(SE_CLASS1 CACHE BOOL "Activate compilation with SE_CLASS1")
set(SE_CLASS3 CACHE BOOL "Activate compilation with SE_CLASS3")
set(ENV{HDF5_ROOT} ${HDF5_ROOT})
set(OPENFPM_PDATA_DIR CACHE PATH "OpenFPM_pdata dirs")
set(PETSC_ROOT CACHE PATH "If compiling with linear algebra indicate the PETSC root directory")
set(ENV{PETSC_DIR} ${PETSC_ROOT})
set(HIP_ENABLE CACHE BOOL "Enable HIP compiler")
set(AMD_ARCH_COMPILE "gfx900" CACHE STRING "AMD gpu architecture used to compile kernels")

add_subdirectory (../openfpm_vcluster openfpm_vcluster)

if (CUDA_ON_BACKEND STREQUAL "CUDA")
        enable_language(CUDA)
endif()

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CUDA_STANDARD 14)

if(NOT TINYOBJLOADER_FOUND)
	set(TINYOBJLOADER_ROOT CACHE PATH "TinyObjLoader library path")
endif()

set(ENV{PATH} "$ENV{PATH}:${HDF5_ROOT}/bin")
set(HDF5_PREFER_PARALLEL TRUE)

set (ENV{BOOST_ROOT} ${BOOST_ROOT})
set (Boost_NO_BOOST_CMAKE OFF)
find_package(Boost 1.72.0 REQUIRED COMPONENTS unit_test_framework iostreams program_options system filesystem
		OPTIONAL_COMPONENTS fiber context)
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)
find_package(TinyObjLoader)
find_package(PETSc)
find_package(OpenMP)
if (CUDA_ON_BACKEND STREQUAL "HIP" AND NOT HIP_FOUND)
        find_package(HIP)
endif()

if(HIP_FOUND)
        set(DEFINE_HIP_GPU "#define HIP_GPU")
        set(DEFINE_CUDIFY_USE_HIP "#define CUDIFY_USE_HIP")
endif()

if(HIP_FOUND)
        set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if(CUDA_FOUND)
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cuda.cu")
else()
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cuda.cpp")
endif()

###### CONFIG.h FILE ######

if(SE_CLASS1)
	set(DEFINE_SE_CLASS1 "#define SE_CLASS1")
endif()

if(SE_CLASS3)
	set(DEFINE_SE_CLASS3 "#define SE_CLASS3")
endif()

if(CUDA_FOUND)
	set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if(MPI_FOUND)
	set(DEFINE_HAVE_MPI "#define HAVE_MPI")
else()
	message( FATAL_ERROR "MPI is required in order to install OpenFPM" )
endif()

if(PETSC_FOUND)
        set(DEFINE_HAVE_PETSC "#define HAVE_PETSC")
endif()

if(HDF5_FOUND)
	if (HDF5_IS_PARALLEL)
		set(DEFINE_HAVE_HDF5 "#define HAVE_HDF5")
	else()
		message( STATUS "HDF5 found ${HDF5_INCLUDE_DIRS} does not have parallel support, OpenFPM require it" )
	endif()
else()
	message( FATAL_ERROR "HDF5 with parallel support is required in order to install OpenFPM" )
endif()

if(TINYOBJLOADER_FOUND)
        set(DEFINE_HAVE_TINYOBJLOADER "#define HAVE_TINYOBJLOADER 1")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/config/config_cmake.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/config/config.h)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})


add_subdirectory (src)

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
	set(DEFINE_HAVE_BOOST ${DEFINE_HAVE_BOOST} PARENT_SCOPE)
	set(DEFINE_HAVE_BOOST_IOSTREAMS ${DEFINE_HAVE_BOOST_IOSTREAMS} PARENT_SCOPE)
	set(DEFINE_HAVE_BOOST_PROGRAM_OPTIONS ${DEFINE_HAVE_BOOST_PROGRAM_OPTIONS} PARENT_SCOPE)
	set(DEFINE_HAVE_BOOST_UNIT_TEST_FRAMEWORK ${DEFINE_HAVE_BOOST_UNIT_TEST_FRAMEWORK} PARENT_SCOPE)
	set(DEFINE_HAVE_BOOST_CONTEXT ${DEFINE_HAVE_BOOST_CONTEXT} PARENT_SCOPE)
	set(DEFINE_HAVE_BOOST_FIBER ${DEFINE_HAVE_BOOST_FIBER} PARENT_SCOPE)
	set(DEFINE_HAVE_OPENMP ${DEFINE_HAVE_OPENMP} PARENT_SCOPE)
	set(DEFINE_HAVE_ALPAKA ${DEFINE_HAVE_ALPAKA} PARENT_SCOPE)
	set(DEFINE_CUDA_GPU ${DEFINE_CUDA_GPU} PARENT_SCOPE)
	set(DEFINE_CUDIFY_BACKEND ${DEFINE_CUDIFY_BACKEND} PARENT_SCOPE)
	set(OPTIONAL_BOOST_LIBS ${OPTIONAL_BOOST_LIBS} PARENT_SCOPE)
	set(WARNING_SUPPRESSION_AND_OPTION_NVCC ${WARNING_SUPPRESSION_AND_OPTION_NVCC} PARENT_SCOPE)
	set(WARNING_SUPPRESSION_AND_OPTION_NVCC_TEXT ${WARNING_SUPPRESSION_AND_OPTION_NVCC_TEXT} PARENT_SCOPE)
        set(CUDA_FOUND ${CUDA_FOUND} PARENT_SCOPE)
        set(CUDA_VERSION_MAJOR ${CUDA_VERSION_MAJOR} PARENT_SCOPE)
        set(CUDA_VERSION_MINOR ${CUDA_VERSION_MINOR} PARENT_SCOPE)
        set(CUDA_VERSION_PATCH ${CUDA_VERSION_PATCH} PARENT_SCOPE)
endif()

