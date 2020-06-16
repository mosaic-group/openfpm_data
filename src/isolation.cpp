#define BOOST_DISABLE_ASSERTS
#define FUSION_MAX_VECTOR_SIZE 20

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "util/math_util_complex.hpp"

#define DISABLE_MPI_WRITTERS

// initialization function:
bool init_unit_test()
{
  return true;
}

std::vector<int> sieve_spf;

// entry point:
int main(int argc, char* argv[])
{
	openfpm::math::init_getFactorization();
	return boost::unit_test::unit_test_main( &init_unit_test, argc, argv );
}

#include <boost/fusion/include/mpl.hpp>

#include <iostream>
#include <typeinfo>

// Include tests


