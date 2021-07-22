#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "util/cuda_util.hpp"
#include "Vector/map_vector.hpp"

#ifndef CUDA_ON_CPU

#ifndef __HIP__
#include "util/cuda/moderngpu/kernel_load_balance.hxx"
#include "util/cuda/moderngpu/kernel_mergesort.hxx"
#include "util/cuda/moderngpu/kernel_reduce.hxx"
#include "util/cuda/moderngpu/kernel_segreduce.hxx"


BOOST_AUTO_TEST_SUITE( modern_gpu_tests )

BOOST_AUTO_TEST_CASE( modern_gpu_loadbalance_lbs )
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

	mgpu::load_balance_search(count, (int *)segments.template getDeviceBuffer<0>(), num_segments, (int *)lbs.template getDeviceBuffer<0>(),context);

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

BOOST_AUTO_TEST_CASE( modern_gpu_sort )
{
	std::cout << "Test modern gpu test tansform_lbs" << "\n";

	mgpu::standard_context_t context(false);

	int count = 200030;

	openfpm::vector_gpu<aggregate<unsigned int,unsigned int>> vgpu;
	openfpm::vector_gpu<aggregate<unsigned int>> gpu_ns;

	vgpu.resize(count);
	gpu_ns.resize(count);

	for (size_t i = 0 ; i < count ; i++)
	{
		vgpu.template get<0>(i) = ((float)rand() / (float)RAND_MAX) * 17;
		vgpu.template get<1>(i) = i;

		gpu_ns.template get<0>(i) = vgpu.template get<0>(i);
	}

	vgpu.hostToDevice<0,1>();

    mergesort((unsigned int *)vgpu.getDeviceBuffer<0>(),(unsigned int *)vgpu.getDeviceBuffer<1>(), count, mgpu::less_t<unsigned int>(), context);

    vgpu.deviceToHost<0,1>();

    // print

    bool match = true;
    for (int i = 0 ; i < count - 1 ; i++)
    {
    	match &= vgpu.template get<0>(i) <= vgpu.template get<0>(i+1);
    	match &= gpu_ns.template get<0>(vgpu.template get<1>(i)) == vgpu.template get<0>(i);
    }

    BOOST_REQUIRE_EQUAL(match,true);

	std::cout << "End test modern gpu test tansform_lbs" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( modern_gpu_reduce )
{
	std::cout << "Test modern gpu reduce" << "\n";

	mgpu::standard_context_t context(false);

	int count = 200030;

	openfpm::vector_gpu<aggregate<int>> vgpu;

	vgpu.resize(count);

	for (size_t i = 0 ; i < count ; i++)
	{
		vgpu.template get<0>(i) = ((float)rand() / (float)RAND_MAX) * 17;
	}

	vgpu.hostToDevice<0>();

	CudaMemory mem;
	mem.allocate(sizeof(int));
	mgpu::reduce((int *)vgpu.template getDeviceBuffer<0>(), count, (int *)mem.getDevicePointer(), mgpu::plus_t<int>(), context);

    mem.deviceToHost();
    int red_p = *(int *)mem.getPointer();

    // print

    int red = 0;
    for (int i = 0 ; i < count ; i++)
    {
    	red += vgpu.template get<0>(i);
    }

    BOOST_REQUIRE_EQUAL(red,red_p);

	std::cout << "End test modern gpu test reduce" << "\n";

	// Test the cell list
}


BOOST_AUTO_TEST_CASE( modern_gpu_seg_reduce )
{
	std::cout << "Test modern gpu segmented reduce" << "\n";

	mgpu::standard_context_t context(false);

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

	mgpu::segreduce((int *)vgpu.template getDeviceBuffer<0>(), vgpu.size(),
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

#endif

#endif

