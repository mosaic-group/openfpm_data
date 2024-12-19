/*
 * SparseGridGpu_performance_get_single.cu
 *
 *  Created on: Sep 9, 2019
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


template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void testGetSingle(std::string testURI, unsigned int i)
{
    auto testName = "Get single";
    constexpr unsigned int dim = 2;
//  constexpr unsigned int blockEdgeSize = 8;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
    typedef aggregate<float> AggregateT;

    unsigned int iterations = 10;

//    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").getSingle");
    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","Get");

    report_sparsegrid_funcs.graphs.put(base + ".dim",dim);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*blockEdgeSize);

    dim3 gridSize(gridEdgeSize, gridEdgeSize);
    dim3 blockSize(blockEdgeSize, blockEdgeSize);
    dim3 blockSizeBlockedInsert(1, 1);
    grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
    SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize> sparseGrid(blockGeometry);
    gpu::ofp_context_t gpuContext;
    sparseGrid.template setBackgroundValue<0>(0);

    // Now fill the grid once
    auto offset = 0;
    sparseGrid.setGPUInsertBuffer(gridSize, blockSizeBlockedInsert);
    insertValues2DBlocked<0, 1, blockEdgeSize> << < gridSize, blockSize >> >
            (sparseGrid.toKernel(), offset, offset);
    sparseGrid.template flush < smax_ < 0 >> (gpuContext, flush_type::FLUSH_ON_DEVICE);

    unsigned long long numElements = gridEdgeSize*blockEdgeSize*gridEdgeSize*blockEdgeSize;
    openfpm::vector<double> measures;

    for (unsigned int iter=0; iter<iterations; ++iter)
    {
        auto offset = 0;

        cudaDeviceSynchronize();

        timer ts;
        ts.start();

        getValues2D<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), offset, offset);
        cudaDeviceSynchronize();

        ts.stop();

        float gElemS = numElements / (1e9 * ts.getwct());
        measures.add(gElemS);
    }

    double mean = 0;
    double deviation = 0;
    standard_deviation(measures,mean,deviation);

    report_sparsegrid_funcs.graphs.put(base + ".Gget.mean",mean);
    report_sparsegrid_funcs.graphs.put(base +".Gget.dev",deviation);

    // All times above are in ms

    std::cout << "Test: " << testName << "\n";
    std::cout << "Block: " << blockEdgeSize << "x" << blockEdgeSize << "\n";
    std::cout << "Grid: " << gridEdgeSize*blockEdgeSize << "x" << gridEdgeSize*blockEdgeSize << "\n";
    double dataOccupancyMean, dataOccupancyDev;
    sparseGrid.deviceToHost();
    sparseGrid.measureBlockOccupancy(dataOccupancyMean, dataOccupancyDev);std::cout << "Data Occupancy: " << dataOccupancyMean << " dev:" << dataOccupancyDev << std::endl;
    report_sparsegrid_funcs.graphs.put(base + ".dataOccupancy.mean",dataOccupancyMean);
    report_sparsegrid_funcs.graphs.put(base +".dataOccupancy.dev",dataOccupancyDev);
    std::cout << "Iterations: " << iterations << "\n";
    std::cout << "Throughput:\n\t" << mean << "GElem/s" << "\n";
}

BOOST_AUTO_TEST_SUITE(performance)

BOOST_AUTO_TEST_SUITE(SparseGridGpu_test)

BOOST_AUTO_TEST_CASE(testGet_gridScaling_2)
{
    std::string testURI = suiteURI + ".device.get.dense.single.2D.2.gridScaling";
    unsigned int counter = 0;
    testGetSingle<2, 128>(testURI, counter++);
    testGetSingle<2, 256>(testURI, counter++);
    testGetSingle<2, 512>(testURI, counter++);
    testGetSingle<2, 1024>(testURI, counter++);
    testGetSingle<2, 2048>(testURI, counter++);
    testGetSingle<2, 4096>(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testGet_gridScaling_4)
{
    std::string testURI = suiteURI + ".device.get.dense.single.2D.4.gridScaling";
    unsigned int counter = 0;
    testGetSingle<4, 64>(testURI, counter++);
    testGetSingle<4, 128>(testURI, counter++);
    testGetSingle<4, 256>(testURI, counter++);
    testGetSingle<4, 512>(testURI, counter++);
    testGetSingle<4, 1024>(testURI, counter++);
    testGetSingle<4, 2048>(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testGet_gridScaling_8)
{
    std::string testURI = suiteURI + ".device.get.dense.single.2D.8.gridScaling";
    unsigned int counter = 0;
    testGetSingle<8, 32>(testURI, counter++);
    testGetSingle<8, 64>(testURI, counter++);
    testGetSingle<8, 128>(testURI, counter++);
    testGetSingle<8, 256>(testURI, counter++);
    testGetSingle<8, 512>(testURI, counter++);
    testGetSingle<8, 1024>(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testGet_gridScaling_16)
{
    std::string testURI = suiteURI + ".device.get.dense.single.2D.16.gridScaling";
    unsigned int counter = 0;
    testGetSingle<16, 16>(testURI, counter++);
    testGetSingle<16, 32>(testURI, counter++);
    testGetSingle<16, 64>(testURI, counter++);
    testGetSingle<16, 128>(testURI, counter++);
    testGetSingle<16, 256>(testURI, counter++);
    testGetSingle<16, 512>(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testGet_gridScaling_32)
{
    std::string testURI = suiteURI + ".device.get.dense.single.2D.32.gridScaling";
    unsigned int counter = 0;
    testGetSingle<32, 8>(testURI, counter++);
    testGetSingle<32, 16>(testURI, counter++);
    testGetSingle<32, 32>(testURI, counter++);
    testGetSingle<32, 64>(testURI, counter++);
    testGetSingle<32, 128>(testURI, counter++);
    testGetSingle<32, 256>(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testGet_blockScaling)
{
    std::string testURI = suiteURI + ".device.get.dense.single.2D.blockScaling";
    unsigned int counter = 0;
    testGetSingle<2, 1024>(testURI, counter++);
    testGetSingle<4, 512>(testURI, counter++);
    testGetSingle<8, 256>(testURI, counter++);
    testGetSingle<16, 128>(testURI, counter++);
    testGetSingle<32, 64>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
