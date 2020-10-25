/*
 * SparseGridGpu_performance_host.cu
 *
 *  Created on: Sep 10, 2019
 *      Author: i-bird
 */

#define BOOST_TEST_DYN_LINK
#define DISABLE_MPI_WRITTERS

#include <boost/test/unit_test.hpp>
#include "performancePlots.hpp"
#include <iostream>
#include "SparseGridGpu/SparseGridGpu.hpp"
#include "SparseGridGpu/tests/utils/SparseGridGpu_util_test.cuh"

extern std::string suiteURI;
extern report_sparse_grid_tests report_sparsegrid_funcs;
extern std::set<std::string> testSet;

BOOST_AUTO_TEST_SUITE(performance)

BOOST_AUTO_TEST_SUITE(SparseGridGpu_test)

BOOST_AUTO_TEST_CASE(testStencilHeatHost_gridScaling)
{

}
BOOST_AUTO_TEST_CASE(testStencilHeatHost_blockScaling)
{

}

BOOST_AUTO_TEST_CASE(testStencilHeatSparseHost_gridScaling)
{

}
BOOST_AUTO_TEST_CASE(testStencilHeatSparseHost_blockScaling)
{

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

