#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_PREREQ(2.59)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CANONICAL_SYSTEM
AC_CONFIG_SRCDIR([src/main.cpp])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([src/config/config.h])
m4_ifdef([ACX_PTHREAD],,[m4_include([m4/acx_pthread.m4])])
m4_ifdef([AX_BOOST],,[m4_include([m4/ax_boost.m4])])
m4_ifdef([ACX_MPI],,[m4_include([m4/acx_mpi.m4])])
m4_ifdef([AX_OPENMP],,[m4_include([m4/ax_openmp.m4])])
m4_ifdef([AX_CUDA],,[m4_include([m4/ax_cuda.m4])])

CXXFLAGS+=" --std=c++11 -march=native -mtune=native -Wno-unused-local-typedefs -Wextra -Wno-unused-parameter "
NVCCFLAGS=" "
INCLUDES_PATH=" "

# Checks for programs.
AC_PROG_CXX

# Checks g++ flags

AC_CANONICAL_HOST

###### Check for test coverage

AC_MSG_CHECKING(whether to build with test coverage)
test_cov=no
AC_ARG_ENABLE(test-coverage,
	AC_HELP_STRING(
	    [--enable-test-coverage],
	    [enable test coverage]
	),
  	test_cov="$enableval"
)



AC_MSG_RESULT($test_cov)
if test x"$test_cov" = x"yes"; then
	AC_DEFINE([TEST_COVERAGE_MODE],[],[Test coverage mode])
	CXXFLAGS="$CXXFLAGS  -fprofile-arcs -ftest-coverage "
fi

###### Check for se-class1

AC_MSG_CHECKING(whether to build with security enhancement class1)
se_class1=no
AC_ARG_ENABLE(se-class1,
        AC_HELP_STRING(
            [--enable-se-class1],
            [enable security enhancement class1]
        ),
        se_class1="$enableval"
)



AC_MSG_RESULT($se_class1)
if test x"$se_class1" = x"yes"; then
        AC_DEFINE([SE_CLASS1],[],[Security enhancement class 1])
fi


###### Check for se-class 2

AC_MSG_CHECKING(whether to build with security enhancement class 2)
se_class2=no
AC_ARG_ENABLE(se-class2,
        AC_HELP_STRING(
            [--enable-se-class2],
            [enable security enhancement class 2]
        ),
        se_class2="$enableval"
)


AC_MSG_RESULT($se_class2)
if test x"$se_class2" = x"yes"; then
        AC_DEFINE([SE_CLASS2],[],[Security enhancement class 2])
fi

###### Check for se-class 3

AC_MSG_CHECKING(whether to build with security enhancement class 3)
se_class3=no
AC_ARG_ENABLE(se-class3,
        AC_HELP_STRING(
            [--enable-se-class3],
            [enable security enhancement class 3]
        ),
        se_class3="$enableval"
)



AC_MSG_RESULT($se_class3)
if test x"$se_class3" = x"yes"; then
        AC_DEFINE([SE_CLASS3],[],[Security enhancement class 3])
fi

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

####### include OpenFPM_devices include path

INCLUDES_PATH+="-I. -Iconfig -I../../OpenFPM_devices/src"

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

##### CHECK FOR BOOST ##############

AX_BOOST([1.52],[],[])

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

##### CHECK FOR BOOST ###

AX_BOOST([1.52],[],[])

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
if [ test x"$test_cov" = x"yes" ]; then
  echo "*    test coverage: yes           *"
else
  echo "*    test coverage: no            *"
fi
if [ test x"$debuger" = x"yes" ]; then
  echo "*    debug: yes                   *"
else
  echo "*    debug: no                    *"
fi
if [ test x"$se_class1" = x"yes" ]; then
  echo "*    se-class1: yes               *"
else
  echo "*    se-class1: no                *"
fi
if [ test x"$se_class2" = x"yes" ]; then
  echo "*    se-class2: yes               *"
else
  echo "*    se-class2: no                *"
fi
if [ test x"$se_class3" = x"yes" ]; then
  echo "*    se-class3: yes               *"
else
  echo "*    se-class3: no                *"
fi
if [ test x"$gpu_support" = x"no" ]; then
  echo "*    gpu: no                      *"
else
  echo "*    gpu: yes                     *"
fi
echo   "*                                 *"
echo   "***********************************"

