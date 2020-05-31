#define BOOST_TEST_DYN_LINK
#include "util/cuda_util.hpp"
#include <boost/test/unit_test.hpp>
#include "memory_ly/memory_conf.hpp"
#include "Vector/map_vector.hpp"

BOOST_AUTO_TEST_SUITE( memory_conf_test )


BOOST_AUTO_TEST_CASE( memory_conf_use )
{
	bool test = is_multiple_buffer_each_prp<openfpm::vector<float>>::value;

	BOOST_REQUIRE_EQUAL(test,false);

	test = is_multiple_buffer_each_prp<float>::value;

	BOOST_REQUIRE_EQUAL(test,false);

	test = is_multiple_buffer_each_prp<openfpm::vector<aggregate<float,float>>>::value;

	BOOST_REQUIRE_EQUAL(test,false);

	test = is_multiple_buffer_each_prp<openfpm::vector_gpu<aggregate<float,float>>>::value;

	BOOST_REQUIRE_EQUAL(test,true);
}

BOOST_AUTO_TEST_SUITE_END()
