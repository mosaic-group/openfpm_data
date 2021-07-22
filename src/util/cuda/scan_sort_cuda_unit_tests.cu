#define BOOST_GPU_ENABLED __host__ __device__

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include "util/cuda_launch.hpp"
#include <boost/test/unit_test.hpp>

#include "util/cuda_util.hpp"
#include "Vector/map_vector.hpp"

#define SORT_WITH_CUB

#include "sort_ofp.cuh"
#include "scan_ofp.cuh"
#include "segreduce_ofp.cuh"

BOOST_AUTO_TEST_SUITE( scan_tests )

BOOST_AUTO_TEST_CASE( test_scan_cub_wrapper )
{
	std::cout << "Test scan CUB" << "\n";

	openfpm::vector_gpu<aggregate<unsigned int>> input;
	openfpm::vector_gpu<aggregate<unsigned int>> output;

	openfpm::vector_gpu<aggregate<unsigned char>> temporal;

	input.resize(10000);
	output.resize(10000);

	// fill input

	for (size_t i = 0 ; i < 10000; i++)
	{
		input.template get<0>(i) = 10.0*(float)rand() / RAND_MAX;
	}

	input.template hostToDevice<0>();

	mgpu::ofp_context_t context;
	openfpm::scan((unsigned int *)input.template getDeviceBuffer<0>(),input.size(),(unsigned int *)output.template getDeviceBuffer<0>(),context);

    output.template deviceToHost<0>();

    size_t cnt = 0;
    for (size_t i = 0 ; i < input.size() ; i++)
    {
    	BOOST_REQUIRE_EQUAL(cnt,output.template get<0>(i));
    	cnt += input.template get<0>(i);
    }

	std::cout << "End scan CUB" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( test_sort_cub_wrapper )
{
	std::cout << "Test sort CUB" << "\n";

	openfpm::vector_gpu<aggregate<unsigned int>> input;
	openfpm::vector_gpu<aggregate<unsigned int>> input_id;

	openfpm::vector_gpu<aggregate<unsigned char>> temporal;

	input.resize(10000);
	input_id.resize(10000);


	// fill input

	for (size_t i = 0 ; i < 10000; i++)
	{
		input.template get<0>(i) = 10000.0*(float)rand() / RAND_MAX;
		input_id.template get<0>(i) = i;
	}

	input.template hostToDevice<0>();
	input_id.template hostToDevice<0>();

	mgpu::ofp_context_t context;

	openfpm::sort((unsigned int *)input.template getDeviceBuffer<0>(),
				  (unsigned int *)input_id.template getDeviceBuffer<0>(),
			      input.size(),mgpu::template less_t<unsigned int>(),context);

	input.template deviceToHost<0>();
	input_id.template deviceToHost<0>();

    for (size_t i = 0 ; i < input.size() - 1 ; i++)
    {
    	BOOST_REQUIRE(input.template get<0>(i) <= input.template get<0>(i+1));
    }

	openfpm::sort((unsigned int *)input.template getDeviceBuffer<0>(),
				  (unsigned int *)input_id.template getDeviceBuffer<0>(),
			      input.size(),mgpu::template greater_t<unsigned int>(),context);

	input.template deviceToHost<0>();
	input_id.template deviceToHost<0>();

    for (size_t i = 0 ; i < input.size() - 1 ; i++)
    {
    	BOOST_REQUIRE(input.template get<0>(i) >= input.template get<0>(i+1));
    }

	std::cout << "End sort CUB" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( test_seg_reduce_wrapper )
{
	std::cout << "Test gpu segmented reduce" << "\n";

	mgpu::ofp_context_t context;

	int count = 130;

	openfpm::vector_gpu<aggregate<int>> vgpu;
	openfpm::vector_gpu<aggregate<int>> segment_offset;
	openfpm::vector_gpu<aggregate<int>> output;
	int init = 0;

	vgpu.resize(count);

	for (size_t i = 0 ; i < count ; i++)
	{
		vgpu.template get<0>(i) = ((float)rand() / (float)RAND_MAX) * 17;
	}

	segment_offset.add();
	segment_offset.template get<0>(0) = 0;
	size_t base = 0;
	while (1)
	{
		int c = ((float)rand() / (float)RAND_MAX) * 17;

		if (c + base >= count)
		{break;}

		segment_offset.add();
		segment_offset.template get<0>(segment_offset.size() - 1) = c + segment_offset.template get<0>(segment_offset.size() - 2);

		base += c;
	}

	vgpu.hostToDevice<0>();

	segment_offset.hostToDevice<0>();
	output.resize(segment_offset.size());

	openfpm::segreduce((int *)vgpu.template getDeviceBuffer<0>(), vgpu.size(),
					(int *)segment_offset.template getDeviceBuffer<0>(), segment_offset.size(),
					(int *)output.template getDeviceBuffer<0>(),
					mgpu::plus_t<int>(), init, context);


	output.template deviceToHost<0>();

	bool match = true;
	size_t i = 0;
	for ( ; i < segment_offset.size()-1 ; i++)
	{
		size_t red = 0;
		for (size_t j = 0 ; j < segment_offset.template get<0>(i+1) - segment_offset.template get<0>(i)  ; j++)
		{
			red += vgpu.template get<0>(segment_offset.template get<0>(i) + j);
		}
		match &= red == output.template get<0>(i);
	}

	BOOST_REQUIRE_EQUAL(match,true);

	size_t red2 = 0;
	for (size_t j = 0 ; j < vgpu.size() - segment_offset.template get<0>(i)  ; j++)
	{
		red2 += vgpu.template get<0>(segment_offset.template get<0>(i) + j);
	}
	match &= red2 == output.template get<0>(i);

	BOOST_REQUIRE_EQUAL(match,true);

	std::cout << "End test modern gpu test reduce" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_SUITE_END()
