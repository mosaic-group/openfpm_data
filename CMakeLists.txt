cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)


list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake_modules/)

set(BOOST_INCLUDE ${Boost_INCLUDE_DIR} CACHE PATH "Include directory for BOOST")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(SE_CLASS1 CACHE BOOL "Activate compilation with SE_CLASS1")
set(SE_CLASS2 CACHE BOOL "Activate compilation with SE_CLASS2")
set(SE_CLASS3 CACHE BOOL "Activate compilation with SE_CLASS3")
set(ENABLE_GPU CACHE BOOL "Disable the GPU code independently that a cuda compiler is found")

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CUDA_STANDARD 11)

if (ENABLE_GPU)
	enable_language(CUDA)
	find_package(CUDA)
endif()

find_package(Boost 1.66.0 REQUIRED COMPONENTS unit_test_framework iostreams program_options)
find_package(LibHilbert REQUIRED)


###### CONFIG.h FILE ######

if(SE_CLASS1)
	set(DEFINE_SE_CLASS1 "#define SE_CLASS1")
endif()

if(SE_CLASS2)
	set(DEFINE_SE_CLASS2 "#define SE_CLASS2")
endif()

if(SE_CLASS3)
	set(DEFINE_SE_CLASS3 "#define SE_CLASS3")
endif()

if(CUDA_FOUND)
	set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if (Boost_FOUND)
	set(DEFINE_HAVE_BOOST "#define HAVE_BOOST")
	set(DEFINE_HAVE_BOOST_IOSTREAMS "#define HAVE_BOOST_IOSTREAMS")
	set(DEFINE_HAVE_BOOST_PROGRAM_OPTIONS "#define HAVE_BOOST_PROGRAM_OPTIONS")
	set(DEFINE_HAVE_BOOST_UNIT_TEST_FRAMEWORK "#define HAVE_BOOST_UNIT_TEST_FRAMEWORK")
else()
	message( FATAL_ERROR "BOOST is required in order to install OpenFPM" )
endif()

if(LIBHILBERT_FOUND)
        set(DEFINE_HAVE_LIBHILBERT "#define HAVE_LIBHILBERT 1")
else()
        message( FATAL_ERROR "LibHilbert is required in order to install OpenFPM")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/config/config_cmake.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/config/config.h)


add_subdirectory (src)

