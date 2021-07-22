#define BOOST_DISABLE_ASSERTS
#define FUSION_MAX_VECTOR_SIZE 20

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "util/math_util_complex.hpp"
#include "Vector/map_vector.hpp"

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

BOOST_AUTO_TEST_SUITE( util_test )

template<typename vtype>
__global__ void test(vtype v)
{
	v.template get<1>(100);
}

BOOST_AUTO_TEST_CASE( map_vector_std_util )
{
	openfpm::vector_gpu<aggregate<int,int,int>> v;

	v.resize(1000);

	CUDA_LAUNCH_DIM3(test,1,1,v.toKernel());
	

	typedef boost::fusion::vector<int,int,int> bfv;
}

BOOST_AUTO_TEST_SUITE_END()

#include <boost/fusion/include/mpl.hpp>

#include <iostream>
#include <typeinfo>

// Include tests


