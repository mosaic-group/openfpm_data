#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_PREREQ(2.59)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CANONICAL_SYSTEM
AC_CONFIG_SRCDIR([src/main.cpp])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([src/config/config.h])
m4_ifdef([ACX_PTHREAD],,[m4_include([m4/acx_pthread.m4])])
m4_ifdef([ACX_MPI],,[m4_include([m4/acx_mpi.m4])])
m4_ifdef([AX_OPENMP],,[m4_include([m4/ax_openmp.m4])])
m4_ifdef([AX_CUDA],,[m4_include([m4/ax_cuda.m4])])
m4_ifdef([AX_BOOST_BASE],,[m4_include([m4/ax_boost_base.m4])])
m4_ifdef([AX_BOOST_IOSTREAMS],,[m4_include([m4/ax_boost_iostreams.m4])])
m4_ifdef([AX_BOOST_PROGRAM_OPTIONS],,[m4_include([m4/ax_boost_program_options.m4])])
m4_ifdef([AX_BOOST_UNIT_TEST_FRAMEWORK],,[m4_include([m4/ax_boost_unit_test_framework.m4])])

CXXFLAGS+=" --std=c++11 -march=native -mtune=native -Wno-unused-local-typedefs -Wextra -Wno-unused-parameter "
NVCCFLAGS=" "
INCLUDES_PATH=" "

# set default prefix
if test "$prefix" = "NONE"; then
   prefix="/usr/local"
fi

# Add openfpm_io to prefix forlder (if prefix end if openfpm_pdata remove that path)
if [ echo "$prefix" | grep -qE "/openfpm_pdata$" ]; then
  prefix=$(echo "$prefix" | sed 's/\(.*\)\/openfpm_pdata$/\1/g')
fi
prefix="$prefix/openfpm_io"
echo "Installation dir is: $prefix"

# Checks for programs.
AC_PROG_CXX

# Checks g++ flags

AC_CANONICAL_HOST

###### Check for debug compilation

AC_MSG_CHECKING(whether to build with debug information)
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [enable debug data generation (def=no)]
	),
  	debuger="$enableval"
)



AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
	AC_DEFINE([DEBUG_MODE],[],[Debug])
	AC_DEFINE([DEBUG],[],[Debug])
	CXXFLAGS="$CXXFLAGS -g3  -Wall -O0 "
	NVCCFLAGS+="$NVCCFLAGS -g -O0 "
else
	CXXFLAGS="$CXXFLAGS -Wall -O3 -g3 -funroll-loops "
	NVCCFLAGS+="$NVCCFLAGS -O3 "
fi

####### Set location of openfpm_pdata

AC_ARG_WITH([pdata],
             AS_HELP_STRING([--with-pdata=directory],
             [specify where is located the pdata project]),
             [pdata_dir="$withval"],
             [pdata_dir="../../openfpm_pdata/src"])

####### include openfpm_devices include path

INCLUDES_PATH+="-I. -Iconfig -I../../openfpm_data/src -I../../openfpm_devices/src -I$pdata_dir "

###### Check for memleak check compilation

AC_MSG_CHECKING(whether to build with memory check)
mem_check=no
AC_ARG_ENABLE(memcheck,
	AC_HELP_STRING(
	    [--enable-memcheck],
	    [enable memory checking (def=no)]
	),
  	mem_check="$enableval"
)



AC_MSG_RESULT($mem_check)
if test x"$mem_check" = x"yes"; then
	AC_DEFINE([MEMLEAK_CHECK],[],[Memory check])
fi

AC_DEFINE([ON_IO_UNIT_TESTS],[],[IO unit tests])

##### CHECK FOR BOOST ##############

AX_BOOST_BASE([1.52],[],[echo "boost not found"
                    exit 202])
AX_BOOST_UNIT_TEST_FRAMEWORK
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_IOSTREAMS

if test x"$ax_cv_boost_unit_test_framework" = x"no"; then
  echo "Notify boost not usable"
  exit 202
fi

if test x"$ax_cv_boost_iostreams" = x"no"; then
  echo "Notify boost not usable"
  exit 202
fi

if test x"$ax_cv_boost_programs_options" = x"no"; then
  echo "Notify boost not usable"
  exit 202
fi

### Unfortunately a lot of linux distros install a pretty old MPI in the system wide folder,
### override such MPI with the installed one is extremely difficult and tricky, because we want
### to include "some" system library but exclude mpi. One possibility is to give specifically
### the wanted libmpi.so directly to the linker. But this is not possible because this lib is
### given by mpic++ in the form -L/path/to/mpi -lmpi, the other is completely eliminate every
### -L with a system default library
###

# eliminate any /usr/lib and and /usr/include from $BOOST_LIB and $BOOST_INCLUDE

BOOST_LDFLAGS=$(echo "$BOOST_LDFLAGS" | sed -e 's/-L\/usr\/lib64//g' | sed -e 's/-L\/usr\/lib//g')
BOOST_CPPFLAGS=$(echo "$BOOST_CPPFLAGS" | sed -e 's/-I\/usr\/include//g')

AC_SUBST(BOOST_LDFLAGS)
AC_SUBST(BOOST_CPPFLAGS)


####### Checking for GPU support

AX_CUDA

if test x"$NVCC_EXIST" = x"yes"; then
  AC_MSG_CHECKING(whether to build with GPU support)
  gpu_support=yes
  AC_ARG_ENABLE(gpu,
        AC_HELP_STRING(
            [--enable-gpu],
            [enable gpu support]
        ),
        gpu_support="$enableval"
  )

  AC_MSG_RESULT($gpu_support)
  if test x"$gpu_support" = x"yes"; then
        AC_DEFINE([CUDA_GPU],[],[CUDA GPU support])
  fi
else
  gpu_support=no
fi


# Set this conditional if cuda is wanted

#

AM_CONDITIONAL(BUILDCUDA, test x$gpu_support = x"yes")

##########################

AC_SUBST(NVCCFLAGS)
AC_SUBST(INCLUDES_PATH)

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
echo ""
echo   "***********************************"
echo   "*                                 *"
if [ test x"$profiler" = x"yes" ]; then
  echo "*    profiler: yes                *"
else
  echo "*    profiler: no                 *"
fi
if [ test x"$memcheck" = x"yes" ]; then
  echo "*    memcheck: yes                *"
else
  echo "*    memcheck: no                 *"
fi
if [ test x"$debuger" = x"yes" ]; then
  echo "*    debug: yes                   *"
else
  echo "*    debug: no                    *"
fi
if [ test x"$gpu_support" = x"no" ]; then
  echo "*    gpu: no                      *"
else
  echo "*    gpu: yes                     *"
fi
echo   "*                                 *"
echo   "***********************************"

