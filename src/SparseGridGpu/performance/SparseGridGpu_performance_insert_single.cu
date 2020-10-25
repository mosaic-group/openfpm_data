/*
 * SparseGridGpu_performance_insert_single.cu
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


template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void testInsertSingle(std::string testURI, unsigned int i)
{
	auto testName = "Insert single (one chunk per element)";
	constexpr unsigned int dim = 2;
//	constexpr unsigned int blockEdgeSize = 8;
	constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
	typedef aggregate<float> AggregateT;

	unsigned int iterations = 10;
	bool prePopulateGrid = true;

//    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").insertSingle");
    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","InsertSingle");

    report_sparsegrid_funcs.graphs.put(base + ".dim",dim);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*blockEdgeSize);

	dim3 gridSize(gridEdgeSize, gridEdgeSize);
	dim3 blockSize(blockEdgeSize, blockEdgeSize);
	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	if (prePopulateGrid)
	{
		// Pre-populate grid
		sparseGrid.setGPUInsertBuffer(gridSize, blockSize);
		insertValues2D<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), 0, 0);
		sparseGrid.template flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
		cudaDeviceSynchronize();
		///
	}

	for (unsigned int iter=0; iter<5; ++iter)
	{
		auto offset = 0;
		sparseGrid.setGPUInsertBuffer(gridSize, blockSize);
		insertValues2D<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), offset, offset);
		sparseGrid.template flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
		cudaDeviceSynchronize();
	}

	unsigned long long numElements = gridEdgeSize*blockEdgeSize*gridEdgeSize*blockEdgeSize;
	openfpm::vector<double> measures;

	for (unsigned int iter=0; iter<iterations; ++iter)
	{
		auto offset = 0;

		cudaDeviceSynchronize();

		timer ts;
		ts.start();

		sparseGrid.setGPUInsertBuffer(gridSize, blockSize);
		insertValues2D<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), offset, offset);
		sparseGrid.template flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
		cudaDeviceSynchronize();

		ts.stop();

		float mElemS = numElements / (1e6 * ts.getwct());
		measures.add(mElemS);
	}

	double mean = 0;
	double deviation = 0;
	standard_deviation(measures,mean,deviation);

    report_sparsegrid_funcs.graphs.put(base + ".Minsert.mean",mean);
    report_sparsegrid_funcs.graphs.put(base +".Minsert.dev",deviation);

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
	std::cout << "Throughput:\n\t" << mean << "M/s" << "\n";
}


BOOST_AUTO_TEST_SUITE(performance)

BOOST_AUTO_TEST_SUITE(SparseGridGpu_test)

BOOST_AUTO_TEST_CASE(testInsert_gridScaling_2)
{
    std::string testURI = suiteURI + ".device.insert.dense.single.2D.2.gridScaling";
    unsigned int counter = 0;
    testInsertSingle<2, 128>(testURI, counter++);
    testInsertSingle<2, 256>(testURI, counter++);
    testInsertSingle<2, 512>(testURI, counter++);
    testInsertSingle<2, 1024>(testURI, counter++);
    testSet.insert(testURI);
}
BOOST_AUTO_TEST_CASE(testInsert_gridScaling_4)
{
    std::string testURI = suiteURI + ".device.insert.dense.single.2D.4.gridScaling";
    unsigned int counter = 0;
    testInsertSingle<4, 64>(testURI, counter++);
    testInsertSingle<4, 128>(testURI, counter++);
    testInsertSingle<4, 256>(testURI, counter++);
    testInsertSingle<4, 512>(testURI, counter++);
    testSet.insert(testURI);
}
BOOST_AUTO_TEST_CASE(testInsert_gridScaling_8)
{
    std::string testURI = suiteURI + ".device.insert.dense.single.2D.8.gridScaling";
    unsigned int counter = 0;
    testInsertSingle<8, 32>(testURI, counter++);
    testInsertSingle<8, 64>(testURI, counter++);
    testInsertSingle<8, 128>(testURI, counter++);
    testInsertSingle<8, 256>(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testInsert_blockScaling)
{
    std::string testURI = suiteURI + ".device.insert.dense.single.2D.blockScaling";
    unsigned int counter = 0;
    testInsertSingle<2, 1024>(testURI, counter++);
    testInsertSingle<4, 512>(testURI, counter++);
    testInsertSingle<8, 256>(testURI, counter++);

    testSet.insert(testURI);
}



BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()



