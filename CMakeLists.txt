cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

if (POLICY CMP0074)
        cmake_policy(SET CMP0074 OLD)
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake_modules/ /opt/rocm/hip/cmake)

set(BOOST_INCLUDE ${Boost_INCLUDE_DIR} CACHE PATH "Include directory for BOOST")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(SE_CLASS1 CACHE BOOL "Activate compilation with SE_CLASS1")
set(SE_CLASS3 CACHE BOOL "Activate compilation with SE_CLASS3")
set(ENABLE_GPU CACHE BOOL "Disable the GPU code independently that a cuda compiler is found")
set(TEST_PERFORMANCE CACHE BOOL "Enable test performance")
set(ALPAKA_ROOT CACHE PATH "Alpaka root path")
set(CUDA_ON_CPU CACHE BOOL "Make Cuda work on heap")
set(HIP_ENABLE CACHE BOOL "Enable HIP compiler")
set(AMD_ARCH_COMPILE "gfx900" CACHE STRING "AMD gpu architecture used to compile kernels")

if (ENABLE_GPU)
	set(CUDA_ON_CPU OFF)
	if (NOT HIP_ENABLE)
		enable_language(CUDA)
		find_package(CUDA)
	endif()
endif()

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CUDA_STANDARD 14)

set(Vc_DIR "${Vc_ROOT}/lib/cmake/Vc/")
message("Searching Vc in ${Vc_DIR}")

find_package(Boost 1.72.0 REQUIRED COMPONENTS unit_test_framework iostreams program_options system filesystem OPTIONAL_COMPONENTS fiber context)
find_package(LibHilbert REQUIRED)
find_package(Vc REQUIRED)
find_package(HIP)

###### CONFIG.h FILE ######

if(SE_CLASS1)
	set(DEFINE_SE_CLASS1 "#define SE_CLASS1")
endif()

if(SE_CLASS3)
	set(DEFINE_SE_CLASS3 "#define SE_CLASS3")
endif()

if(TEST_PERFORMANCE)
        set(DEFINE_PERFORMANCE_TEST "#define PERFORMANCE_TEST")
endif()

if(HIP_FOUND)
        set(DEFINE_HIP_GPU "#define HIP_GPU")
        set(DEFINE_CUDIFY_USE_HIP "#define CUDIFY_USE_HIP")
endif()

if (Boost_FOUND)
        set(DEFINE_HAVE_BOOST "#define HAVE_BOOST")
        set(DEFINE_HAVE_BOOST_IOSTREAMS "#define HAVE_BOOST_IOSTREAMS")
        set(DEFINE_HAVE_BOOST_PROGRAM_OPTIONS "#define HAVE_BOOST_PROGRAM_OPTIONS")
        set(DEFINE_HAVE_BOOST_UNIT_TEST_FRAMEWORK "#define HAVE_BOOST_UNIT_TEST_FRAMEWORK")
        if (Boost_CONTEXT_FOUND)
                set(DEFINE_HAVE_BOOST_CONTEXT "#define HAVE_BOOST_CONTEXT")
        else()
                #if context is not there CUDA_ON_CPU cannot be activated
                set(CUDA_ON_CPU OFF)
        endif()
        if (Boost_FIBER_FOUND)
                set(DEFINE_HAVE_BOOST_FIBER "#define HAVE_BOOST_FIBER")
        endif()
else()
	message( FATAL_ERROR "BOOST is required in order to install OpenFPM" )
endif()

if(ENABLE_GPU AND CUDA_FOUND)
        set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if(CUDA_ON_CPU)
        set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if(HIP_FOUND)
	set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if(LIBHILBERT_FOUND)
        set(DEFINE_HAVE_LIBHILBERT "#define HAVE_LIBHILBERT 1")
else()
        message( FATAL_ERROR "LibHilbert is required in order to install OpenFPM")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/config/config_cmake.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/config/config.h)


add_subdirectory (../openfpm_devices/ openfpm_devices)

set(WARNING_SUPPRESSION_AND_OPTION_NVCC ${WARNING_SUPPRESSION_AND_OPTION_NVCC} PARENT_SCOPE)
set(WARNING_SUPPRESSION_AND_OPTION_NVCC_TEXT ${WARNING_SUPPRESSION_AND_OPTION_NVCC_TEXT} PARENT_SCOPE)

add_subdirectory (src)
