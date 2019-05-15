/*
 * cub_gpu_tests.cu
 *
 *  Created on: May 15, 2019
 *      Author: i-bird
 */

#ifndef CUB_GPU_TESTS_CU_
#define CUB_GPU_TESTS_CU_

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "cub/cub.cuh"
#include "Vector/map_vector.hpp"

BOOST_AUTO_TEST_SUITE( cub_gpu_tests )

BOOST_AUTO_TEST_CASE( cub_gpu_scan_test )
{
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

    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;
    cub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes,(unsigned int *)input.template getDeviceBuffer<0>(),
    		                                                    (unsigned int *)output.template getDeviceBuffer<0>(),
    		                                                    input.size());

    temporal.resize(temp_storage_bytes);

    // Run
    cub::DeviceScan::ExclusiveSum(temporal.template getDeviceBuffer<0>(), temp_storage_bytes,(unsigned int *)input.template getDeviceBuffer<0>(),
            (unsigned int *)output.template getDeviceBuffer<0>(),
            input.size());

    // Check

    output.template deviceToHost<0>();

    size_t cnt = 0;
    for (size_t i = 0 ; i < input.size() ; i++)
    {
    	BOOST_REQUIRE_EQUAL(cnt,output.template get<0>(i));
    	cnt += input.template get<0>(i);
    }
}


BOOST_AUTO_TEST_SUITE_END()

#endif /* CUB_GPU_TESTS_CU_ */
