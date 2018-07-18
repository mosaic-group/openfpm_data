#define BOOST_GPU_ENABLED __host__ __device__

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "util/cuda_util.hpp"
#include "Vector/map_vector.hpp"
#include "util/cuda/moderngpu/kernel_load_balance.hxx"

BOOST_AUTO_TEST_SUITE( modern_gpu_tests )

BOOST_AUTO_TEST_CASE( make_vector_compatible_with_modern_gpu_transform_lbs )
{
	std::cout << "Test modern gpu test tansform_lbs" << "\n";

	mgpu::standard_context_t context(false);

	int count = 200030;
	int spacing = 100;

	int num_segments = mgpu::div_up(count, spacing);
	openfpm::vector_gpu<aggregate<int>> segments(num_segments);
	for(int i = 0; i < num_segments; ++i)
	{segments.template get<0>(i) = i * spacing;}

	openfpm::vector_gpu<aggregate<int>>  lbs(count);

	segments.template hostToDevice<0>();

	load_balance_search(count, (int *)segments.template getDeviceBuffer<0>(), num_segments, (int *)lbs.template getDeviceBuffer<0>(),context);

	lbs.deviceToHost<0>();

	bool check = true;
	for(size_t i = 0; i < lbs.size(); ++i)
	{
	    check &= lbs.template get<0>(i) == i / spacing;
	}

	BOOST_REQUIRE_EQUAL(check,true);

	std::cout << "End test modern gpu test tansform_lbs" << "\n";

	// Test the cell list
}


BOOST_AUTO_TEST_SUITE_END()

