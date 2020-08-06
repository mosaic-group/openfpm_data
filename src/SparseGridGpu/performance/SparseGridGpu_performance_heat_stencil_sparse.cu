/*
 * SparseGridGpu_performance_heat_stencil_sparse.cu
 *
 *  Created on: Sep 10, 2019
 *      Author: i-bird
 */

#define SCAN_WITH_CUB
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

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
void testStencilHeatSparse_perf(unsigned int i, std::string base, float fillMultiplier=1, float voidMultiplier=1)
{
    auto testName = "In-place sparse stencil";
//    unsigned int gridEdgeSize = 128;
    constexpr unsigned int dim = SparseGridZ::dims;
//    const unsigned int blockEdgeSize = SparseGridZ::blockEdgeSize_;

    typedef HeatStencil<dim, 0, 1> Stencil01T;
    typedef HeatStencil<dim, 1, 0> Stencil10T;

//    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").stencil");

    report_sparsegrid_funcs.graphs.put(base + ".dim",2);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*blockEdgeSize);

    unsigned int iterations = 100;

    openfpm::vector<double> measures_gf;
    openfpm::vector<double> measures_tm;

    dim3 gridSize(gridEdgeSize, gridEdgeSize);
    dim3 blockSize(blockEdgeSize,blockEdgeSize);
    unsigned int spatialEdgeSize = 1000000;
    size_t sz[2] = {spatialEdgeSize, spatialEdgeSize};
    typename SparseGridZ::grid_info blockGeometry(sz);
    SparseGridZ sparseGrid(blockGeometry);
    mgpu::ofp_context_t ctx;
    sparseGrid.template setBackgroundValue<0>(0);

    ///// Insert sparse content, a set of concentric spheres /////
    float allMultiplier = fillMultiplier + voidMultiplier;
    const unsigned int numSpheres = gridEdgeSize / (2*allMultiplier);
//    const unsigned int numSpheres = 1;
    unsigned int centerPoint = spatialEdgeSize / 2;

    for (int i = 1; i <= numSpheres; ++i)
    {
        unsigned int rBig = allMultiplier*i * blockEdgeSize;
        unsigned int rSmall = (allMultiplier*i - fillMultiplier) * blockEdgeSize;
        // Sphere i-th
        grid_key_dx<dim, int> start1({centerPoint, centerPoint});
        sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
        CUDA_LAUNCH_DIM3((insertSphere<0>),
                         gridSize, dim3(blockEdgeSize * blockEdgeSize, 1, 1),
                         sparseGrid.toKernel(), start1, rBig, rSmall, 5);
        cudaDeviceSynchronize();
        sparseGrid.template flush<smax_<0 >>(ctx, flush_type::FLUSH_ON_DEVICE);
        cudaDeviceSynchronize();
    }
    ///// /////

    sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
    sparseGrid.tagBoundaries(ctx);

    sparseGrid.template deviceToHost<0>(); // NECESSARY as count takes place on Host!
    auto existingElements = sparseGrid.countExistingElements();
    auto boundaryElements = sparseGrid.countBoundaryElements();
    unsigned long long numElements = existingElements - boundaryElements;

    // Now apply some boundary conditions
    sparseGrid.template applyStencils<BoundaryStencilSetXRescaled<dim,0,0>>(STENCIL_MODE_INPLACE,
            centerPoint, centerPoint + 2*blockEdgeSize*gridEdgeSize,
            0.0, 10.0);

    iterations /= 2;
    for (unsigned int iter=0; iter<iterations; ++iter)
    {
        cudaDeviceSynchronize();

        timer ts;
        ts.start();

        sparseGrid.template applyStencils<Stencil01T>(STENCIL_MODE_INPLACE, 0.1);
        cudaDeviceSynchronize();
        sparseGrid.template applyStencils<Stencil10T>(STENCIL_MODE_INPLACE, 0.1);
        cudaDeviceSynchronize();

        ts.stop();

        measures_tm.add(ts.getwct());

        float gElemS = 2 * numElements / (1e9 * ts.getwct());
        float gFlopsS = gElemS * Stencil01T::flops;

        measures_gf.add(gFlopsS);
    }

    double mean_tm = 0;
    double deviation_tm = 0;
    standard_deviation(measures_tm,mean_tm,deviation_tm);

    double mean_gf = 0;
    double deviation_gf = 0;
    standard_deviation(measures_gf,mean_gf,deviation_gf);

    // All times above are in ms

    float gElemS =  2 * numElements / (1e9 * mean_tm);
    float gFlopsS = gElemS * Stencil01T::flops;
    std::cout << "Test: " << testName << std::endl;
    std::cout << "Block: " << blockEdgeSize << "x" << blockEdgeSize << std::endl;
    std::cout << "Grid: " << gridEdgeSize*blockEdgeSize << "x" << gridEdgeSize*blockEdgeSize << std::endl;
    double dataOccupancyMean, dataOccupancyDev;
    sparseGrid.deviceToHost();
    sparseGrid.measureBlockOccupancy(dataOccupancyMean, dataOccupancyDev);
    std::cout << "Data Occupancy: " << dataOccupancyMean << " dev:" << dataOccupancyDev << std::endl;
    report_sparsegrid_funcs.graphs.put(base + ".dataOccupancy.mean",dataOccupancyMean);
    report_sparsegrid_funcs.graphs.put(base +".dataOccupancy.dev",dataOccupancyDev);
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\tStencil: " << mean_gf << " dev:" << deviation_gf << " s" << std::endl;
    std::cout << "Throughput: " << std::endl << "\t " << gElemS << " GElem/s " << std::endl << "\t " << gFlopsS << " GFlops/s" << std::endl;

    report_sparsegrid_funcs.graphs.put(base + ".GFlops.mean",mean_gf);
    report_sparsegrid_funcs.graphs.put(base +".GFlops.dev",deviation_gf);
    report_sparsegrid_funcs.graphs.put(base + ".time.mean",mean_tm);
    report_sparsegrid_funcs.graphs.put(base +".time.dev",deviation_tm);
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilHeatSparse_perf(std::string testURI, unsigned int i,
        float fillMultiplier=1, float voidMultiplier=1, std::string occupancyStr="05")
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilNSparse"+occupancyStr);

    testStencilHeatSparse_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize, long int>>(i, base,
                    fillMultiplier, voidMultiplier);
    cudaDeviceSynchronize();
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilHeatSparseZ_perf(std::string testURI, unsigned int i,
                                       float fillMultiplier=1, float voidMultiplier=1, std::string occupancyStr="05")
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilNSparse"+occupancyStr);

    testStencilHeatSparse_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu_z<dim, AggregateT, blockEdgeSize, chunkSize, long int>>(i, base,
                                                                                fillMultiplier, voidMultiplier);
    cudaDeviceSynchronize();
}

BOOST_AUTO_TEST_SUITE(performance)

BOOST_AUTO_TEST_SUITE(SparseGridGpu_test)

BOOST_AUTO_TEST_CASE(testStencilHeatSparse05_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.05.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeatSparse_perf<blockEdgeSize, 128>(testURI, counter++, 1.45, 1, "05");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 256>(testURI, counter++, 1.45, 1, "05");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 512>(testURI, counter++, 1.45, 1, "05");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 1024>(testURI, counter++, 1.45, 1, "05");
////    launch_testStencilHeatSparse_perf<blockEdgeSize, 2048>(testURI, counter++, 1.45, 1, "05);

    testSet.insert(testURI);
}


BOOST_AUTO_TEST_CASE(testStencilHeatSparse08_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.08.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeatSparse_perf<blockEdgeSize, 128>(testURI, counter++, 2, 0.20, "08");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 256>(testURI, counter++, 2, 0.20, "08");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 512>(testURI, counter++, 2, 0.20, "08");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 1024>(testURI, counter++, 2, 0.20, "08");

    testSet.insert(testURI);
}


BOOST_AUTO_TEST_CASE(testStencilHeatSparse09_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.09.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeatSparse_perf<blockEdgeSize, 128>(testURI, counter++, 2.3, 0.07, "09");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 256>(testURI, counter++, 2.3, 0.07, "09");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 512>(testURI, counter++, 2.3, 0.07, "09");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 1024>(testURI, counter++, 2.3, 0.07, "09");

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatSparseZ05_gridScaling)
{
	std::string testURI = suiteURI + ".device.stencil.sparse.Z.2D.05.gridScaling";
	unsigned int counter = 0;
	constexpr unsigned int blockEdgeSize = 8;
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 128>(testURI, counter++, 1.45, 1, "05");
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 256>(testURI, counter++, 1.45, 1, "05");
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 512>(testURI, counter++, 1.45, 1, "05");
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 1024>(testURI, counter++, 1.45, 1, "05");

	testSet.insert(testURI);
}


BOOST_AUTO_TEST_CASE(testStencilHeatSparseZ08_gridScaling)
{
	std::string testURI = suiteURI + ".device.stencil.sparse.Z.2D.08.gridScaling";
	unsigned int counter = 0;
	constexpr unsigned int blockEdgeSize = 8;
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 128>(testURI, counter++, 2, 0.20, "08");
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 256>(testURI, counter++, 2, 0.20, "08");
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 512>(testURI, counter++, 2, 0.20, "08");
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 1024>(testURI, counter++, 2, 0.20, "08");

	testSet.insert(testURI);
}


BOOST_AUTO_TEST_CASE(testStencilHeatSparseZ09_gridScaling)
{
	std::string testURI = suiteURI + ".device.stencil.sparse.Z.2D.09.gridScaling";
	unsigned int counter = 0;
	constexpr unsigned int blockEdgeSize = 8;
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 128>(testURI, counter++, 2.3, 0.07, "09");
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 256>(testURI, counter++, 2.3, 0.07, "09");
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 512>(testURI, counter++, 2.3, 0.07, "09");
	launch_testStencilHeatSparseZ_perf<blockEdgeSize, 1024>(testURI, counter++, 2.3, 0.07, "09");

	testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatSparse05_32Block_2048Grid_Case)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.05.32_2048";
    unsigned int counter = 0;
    launch_testStencilHeatSparse_perf<32, 2048/32>(testURI, counter++, 1.45, 1, "05");

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
