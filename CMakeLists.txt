cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)


list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake_modules/)

set(BOOST_INCLUDE ${Boost_INCLUDE_DIR} CACHE PATH "Include directory for BOOST")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(SE_CLASS1 CACHE BOOL "Activate compilation with SE_CLASS1")
set(SE_CLASS2 CACHE BOOL "Activate compilation with SE_CLASS2")
set(SE_CLASS3 CACHE BOOL "Activate compilation with SE_CLASS3")
set(ENABLE_GPU CACHE BOOL "Disable the GPU code independently that a cuda compiler is found")
set(TEST_PERFORMANCE CACHE BOOL "Enable test performance")
set(HIP_ENABLE CACHE BOOL "Enable HIP compiler")

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CUDA_STANDARD 11)

if (ENABLE_GPU)
	enable_language(CUDA)
	find_package(CUDA)

        if (CUDA_VERSION_MAJOR EQUAL 9 AND CUDA_VERSION_MINOR EQUAL 2)
                message("CUDA is compatible")
                set(WARNING_SUPPRESSION_AND_OPTION_NVCC  -Xcudafe "--display_error_number --diag_suppress=611 --diag_suppress=2885 --diag_suppress=2886  --diag_suppress=2887  --diag_suppress=2888 --diag_suppress=186 --diag_suppress=111" --expt-extended-lambda)
                set(WARNING_SUPPRESSION_AND_OPTION_NVCC_STR  "-Xcudafe --display_error_number -Xcudafe --diag_suppress=611 -Xcudafe --diag_suppress=2885 -Xcudafe --diag_suppress=2886 -Xcudafe --diag_suppress=2887 -Xcudafe  --diag_suppress=2888 -Xcudafe --diag_suppress=186 -Xcudafe --diag_suppress=111 --expt-extended-lambda")
                FILE(WRITE cuda_options " -Xcudafe \"--display_error_number --diag_suppress=611 --diag_suppress=2885 --diag_suppress=2886  --diag_suppress=2887  --diag_suppress=2888 --diag_suppress=186 --diag_suppress=111\" --expt-extended-lambda ")
        elseif ( CUDA_VERSION_MAJOR EQUAL 10 AND CUDA_VERSION_MINOR EQUAL 1 )
                message("CUDA is compatible")
                set(WARNING_SUPPRESSION_AND_OPTION_NVCC  -Xcudafe "--display_error_number --diag_suppress=2915 --diag_suppress=2912 --diag_suppress=2913 --diag_suppress=111 --diag_suppress=186 --diag_suppress=611 " --expt-extended-lambda )
                set(WARNING_SUPPRESSION_AND_OPTION_NVCC_STR " -Xcudafe --display_error_number -Xcudafe --diag_suppress=2915 -Xcudafe --diag_suppress=2912 -Xcudafe --diag_suppress=2913 -Xcudafe --diag_suppress=111 -Xcudafe --diag_suppress=186 -Xcudafe --diag_suppress=611  --expt-extended-lambda ")
                FILE(WRITE cuda_options "-Xcudafe \"--display_error_number --diag_suppress=2915 --diag_suppress=2914  --diag_suppress=2912 --diag_suppress=2913 --diag_suppress=111 --diag_suppress=186 --diag_suppress=611 \" --expt-extended-lambda")
        else()
                message(FATAL_ERROR "CUDA is incompatible, version 9.2 and 10.1 is only supported")
        endif()

endif()

find_package(Boost 1.66.0 REQUIRED COMPONENTS unit_test_framework iostreams program_options system filesystem)
find_package(LibHilbert REQUIRED)
find_package(HIP)
if(HIP_FOUND)
    message(STATUS "Found HIP: " ${HIP_VERSION})
endif()

if (HIP_ENABLE)
        if(NOT HIP_FOUND)
                message( FATAL_ERROR "HIP has not been found" )
        endif()
        if (HIP_ENABLE)
                set(CMAKE_CXX_EXTENSIONS OFF)
                string(APPEND CMAKE_CXX_FLAGS ${WARNING_SUPPRESSION_AND_OPTION_NVCC_STR})
        endif()
endif()


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

if(TEST_PERFORMANCE)
        set(DEFINE_PERFORMANCE_TEST "#define PERFORMANCE_TEST")
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

