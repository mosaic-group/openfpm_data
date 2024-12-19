cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_data LANGUAGES C CXX)

if (POLICY CMP0074)
	cmake_policy(SET CMP0074 OLD)
endif ()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)

###### CONFIG.h FILE ######

if(LIBHILBERT_FOUND)
	set(DEFINE_HAVE_LIBHILBERT "#define HAVE_LIBHILBERT 1")
else()
	message( FATAL_ERROR "LibHilbert is required in order to install OpenFPM")
endif()

if (Vc_FOUND)
	set(DEFINE_HAVE_VCDEVEL)
else()
	message( FATAL_ERROR "Vc is required in roder to install OpenFPM")
endif()

add_subdirectory(src)

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
	set(DEFINE_HAVE_VCDEVEL ${DEFINE_HAVE_VCDEVEL} CACHE INTERNAL "")
	set(DEFINE_HAVE_LIBHILBERT ${DEFINE_HAVE_LIBHILBERT} CACHE INTERNAL "")
endif()

