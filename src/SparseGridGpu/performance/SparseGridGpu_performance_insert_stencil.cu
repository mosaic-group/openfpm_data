/*
 * SparseGridGpu_performance_insert_stencil.cu
 *
 *  Created on: Sep 10, 2019
 *      Author: i-bird
 */
#define SCAN_WITH_CUB
#define BOOST_TEST_DYN_LINK
#define OPENFPM_DATA_ENABLE_IO_MODULE
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
void testInsertStencil(std::string testURI, unsigned int i)
{
	auto testName = "Insert stencil";
	constexpr unsigned int dim = 2;
//	constexpr unsigned int blockEdgeSize = 8;
	constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
	typedef aggregate<float> AggregateT;
	typedef HeatStencil<dim,0,0> StencilT;

	unsigned int iterations = 10;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilInsertN");

    report_sparsegrid_funcs.graphs.put(base + ".dim",2);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*blockEdgeSize);

	dim3 gridSize(gridEdgeSize, gridEdgeSize);
	dim3 blockSize(blockEdgeSize, blockEdgeSize);
	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// Initialize the grid
	sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
	CUDA_LAUNCH_DIM3((insertConstantValue<0>),gridSize, blockSize,sparseGrid.toKernel(), 0);
	sparseGrid.template flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
	dim3 sourcePt(gridSize.x * blockEdgeSize / 2, gridSize.y * blockEdgeSize / 2, 0);
	insertOneValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), sourcePt, 100);
	sparseGrid.template flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

	unsigned long long numElements = gridEdgeSize*blockEdgeSize*gridEdgeSize*blockEdgeSize;

	for (unsigned int iter=0; iter<5; ++iter)
	{
		sparseGrid.template applyStencils<StencilT>(STENCIL_MODE_INSERT, 0.1);
		sparseGrid.template flush<smax_<0>>(ctx, flush_type::FLUSH_ON_DEVICE);
	}

    openfpm::vector<double> gElemSMeasures;
    openfpm::vector<double> gFlopsSMeasures;

    for (unsigned int iter=0; iter<iterations; ++iter)
	{
		timer ts;
		ts.start();

		cudaDeviceSynchronize();

		sparseGrid.template applyStencils<StencilT>(STENCIL_MODE_INSERT, 0.1);

		cudaDeviceSynchronize();

		ts.stop();

		float gElemS = numElements / (1e9 * ts.getwct());
		float gFlopsS = gElemS * StencilT::flops;

		gElemSMeasures.add(gElemS);
		gFlopsSMeasures.add(gFlopsS);
	}


	double elemMean=0, elemDeviation=0;
	standard_deviation(gElemSMeasures, elemMean, elemDeviation);
    report_sparsegrid_funcs.graphs.put(base + ".GElems.mean",elemMean);
    report_sparsegrid_funcs.graphs.put(base +".GElems.dev",elemDeviation);
    double flopsMean=0, flopsDeviation=0;
    standard_deviation(gFlopsSMeasures, flopsMean, flopsDeviation);
    report_sparsegrid_funcs.graphs.put(base + ".GFlops.mean",flopsMean);
    report_sparsegrid_funcs.graphs.put(base +".GFlops.dev",flopsDeviation);

	std::cout << "Test: " << testName << "\n";
	std::cout << "Block: " << blockEdgeSize << "x" << blockEdgeSize << "\n";
	std::cout << "Grid: " << gridEdgeSize*blockEdgeSize << "x" << gridEdgeSize*blockEdgeSize << "\n";
    double dataOccupancyMean, dataOccupancyDev;
    sparseGrid.deviceToHost();
    sparseGrid.measureBlockOccupancy(dataOccupancyMean, dataOccupancyDev);std::cout << "Data Occupancy: " << dataOccupancyMean << " dev:" << dataOccupancyDev << std::endl;
    report_sparsegrid_funcs.graphs.put(base + ".dataOccupancy.mean",dataOccupancyMean);
    report_sparsegrid_funcs.graphs.put(base +".dataOccupancy.dev",dataOccupancyDev);
    std::cout << "Iterations: " << iterations << "\n";
	std::cout << "Throughput:\n\t" << elemMean << " GElem/s dev: " << elemDeviation << " GElem/s" << std::endl
	            << "\t" << flopsMean << " GFlops/s dev: " << flopsDeviation << " GFlops/s" << std::endl;
}


BOOST_AUTO_TEST_SUITE(performance)

BOOST_AUTO_TEST_SUITE(SparseGridGpu_test)

BOOST_AUTO_TEST_CASE(testStencilHeatInsert_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencilInsert.dense.N.2D.gridScaling";
    unsigned int counter = 0;
    testInsertStencil<8, 64>(testURI, counter++);
    testInsertStencil<8, 128>(testURI, counter++);
    testInsertStencil<8, 256>(testURI, counter++);
    testInsertStencil<8, 512>(testURI, counter++);
    testInsertStencil<8, 1024>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatInsert_gridScaling_8)
{
    std::string testURI = suiteURI + ".device.stencilInsert.dense.N.2D.8.gridScaling";
    unsigned int counter = 0;
    testInsertStencil<8, 64>(testURI, counter++);
	testInsertStencil<8, 128>(testURI, counter++);
	testInsertStencil<8, 256>(testURI, counter++);
	testInsertStencil<8, 512>(testURI, counter++);
	testInsertStencil<8, 1024>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatInsert_gridScaling_16)
{
    std::string testURI = suiteURI + ".device.stencilInsert.dense.N.2D.16.gridScaling";
    unsigned int counter = 0;
    testInsertStencil<16, 32>(testURI, counter++);
    testInsertStencil<16, 64>(testURI, counter++);
    testInsertStencil<16, 128>(testURI, counter++);
    testInsertStencil<16, 256>(testURI, counter++);
    testInsertStencil<16, 512>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatInsert_blockScaling)
{
    std::string testURI = suiteURI + ".device.stencilInsert.dense.N.2D.blockScaling";
    unsigned int counter = 0;
    testInsertStencil<4, 1024>(testURI, counter++);
    testInsertStencil<8, 512>(testURI, counter++);
    testInsertStencil<16, 256>(testURI, counter++);
    testInsertStencil<32, 128>(testURI, counter++);

    testSet.insert(testURI);
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
