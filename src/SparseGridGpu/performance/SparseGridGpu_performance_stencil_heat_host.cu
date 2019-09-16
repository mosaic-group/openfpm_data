/*
 * SparseGridGpu_performance_host.cu
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

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
void testStencilHeatSparseHost_perf(unsigned int i, std::string base, float fillMultiplier=1, float voidMultiplier=1)
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
    sparseGrid.tagBoundaries();

    sparseGrid.template deviceToHost<0>(); // NECESSARY as count takes place on Host!
    auto existingElements = sparseGrid.countExistingElements();
    auto boundaryElements = sparseGrid.countBoundaryElements();
    unsigned long long numElements = existingElements - boundaryElements;

    // Now apply some boundary conditions
    sparseGrid.template applyBoundaryStencils<BoundaryStencilSetXRescaled<dim,0,0>>(
            centerPoint, centerPoint + 2*blockEdgeSize*gridEdgeSize,
            0.0, 10.0);
    cudaDeviceSynchronize();

    sparseGrid.template deviceToHost<0>(); // NECESSARY as stencils are applied on Host!

    iterations /= 2;
    for (unsigned int iter=0; iter<iterations; ++iter)
    {
        cudaDeviceSynchronize();

        timer ts;
        ts.start();

        sparseGrid.template applyStencilsHost<Stencil01T>(STENCIL_MODE_INPLACE, 0.1);
        cudaDeviceSynchronize();
        sparseGrid.template applyStencilsHost<Stencil10T>(STENCIL_MODE_INPLACE, 0.1);
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
    sparseGrid.measureBlockOccupancy(dataOccupancyMean, dataOccupancyDev);std::cout << "Data Occupancy: " << dataOccupancyMean << " dev:" << dataOccupancyDev << std::endl;
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


template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
void testStencilHeatHost_perf(unsigned int i, std::string base)
{
    // todo: Make sure to reimplement the host stencil application function to pre-load to a block of memory both content and ghost
    // this way we can avoid binary searches...
    auto testName = "In-place stencil HOST";
    typedef HeatStencil<SparseGridZ::dims,0,1> Stencil01T;
    typedef HeatStencil<SparseGridZ::dims,1,0> Stencil10T;


    constexpr unsigned int dim = 2;

//    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").stencil");

    report_sparsegrid_funcs.graphs.put(base + ".dim",dim);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*SparseGridZ::blockEdgeSize_);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*SparseGridZ::blockEdgeSize_);

//    unsigned int iterations = 100;
    unsigned int iterations = 10;
//    unsigned int iterations = 2;
//    unsigned int iterations = 1; // Debug

    openfpm::vector<double> measures_gf;
    openfpm::vector<double> measures_tm;

    dim3 gridSize(gridEdgeSize, gridEdgeSize);
    dim3 blockSize(SparseGridZ::blockEdgeSize_,SparseGridZ::blockEdgeSize_);
    typename SparseGridZ::grid_info blockGeometry(gridSize);
    SparseGridZ sparseGrid(blockGeometry);
    mgpu::ofp_context_t ctx;
    sparseGrid.template setBackgroundValue<0>(0);

    unsigned long long numElements = gridEdgeSize*SparseGridZ::blockEdgeSize_*gridEdgeSize*SparseGridZ::blockEdgeSize_;

    // Initialize the grid
    sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
    CUDA_LAUNCH_DIM3((insertConstantValue<0>),gridSize, blockSize,sparseGrid.toKernel(), 0);
    sparseGrid.template flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
    dim3 sourcePt(gridSize.x * SparseGridZ::blockEdgeSize_ / 2, gridSize.y * SparseGridZ::blockEdgeSize_ / 2, 0);
    insertOneValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), sourcePt, 100);
    sparseGrid.template flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
    cudaDeviceSynchronize();

    sparseGrid.template deviceToHost<0,1>();

    iterations /= 2;
    for (unsigned int iter=0; iter<iterations; ++iter)
    {
        cudaDeviceSynchronize();

        timer ts;
        ts.start();

        sparseGrid.template applyStencilsHost<Stencil01T>(STENCIL_MODE_INPLACE, 0.1);
        cudaDeviceSynchronize();
        sparseGrid.template applyStencilsHost<Stencil10T>(STENCIL_MODE_INPLACE, 0.1);
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

    float gElemS = 2 * numElements / (1e9 * mean_tm);
    float gFlopsS = gElemS * Stencil01T::flops;

    std::cout << "Test: " << testName << std::endl;
    std::cout << "Host: " << SparseGridZ::blockEdgeSize_ << "x" << SparseGridZ::blockEdgeSize_ << std::endl;
    std::cout << "Grid: " << gridEdgeSize*SparseGridZ::blockEdgeSize_ << "x" << gridEdgeSize*SparseGridZ::blockEdgeSize_ << std::endl;
    double dataOccupancyMean=0, dataOccupancyDev=0;
    sparseGrid.deviceToHost();
    sparseGrid.measureBlockOccupancy(dataOccupancyMean, dataOccupancyDev);std::cout << "Data Occupancy: " << dataOccupancyMean << " dev:" << dataOccupancyDev << std::endl;
    report_sparsegrid_funcs.graphs.put(base + ".dataOccupancy.mean",dataOccupancyMean);
    report_sparsegrid_funcs.graphs.put(base +".dataOccupancy.dev",dataOccupancyDev);
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\tStencil: " << mean_gf << " dev:" << deviation_gf << " s" << std::endl;
    std::cout << "Throughput: " << std::endl << "\t " << gElemS << " GElem/s " << std::endl
                << "\t " << gFlopsS << " GFlops/s" << std::endl;

    report_sparsegrid_funcs.graphs.put(base + ".GFlops.mean",mean_gf);
    report_sparsegrid_funcs.graphs.put(base +".GFlops.dev",deviation_gf);
    report_sparsegrid_funcs.graphs.put(base + ".time.mean",mean_tm);
    report_sparsegrid_funcs.graphs.put(base +".time.dev",deviation_tm);
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilHeatSparseHost_perf(std::string testURI, unsigned int i,
        float fillMultiplier=1, float voidMultiplier=1, std::string occupancyStr="05")
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilNSparseHost"+occupancyStr);

    testStencilHeatSparseHost_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize, long int>>(i, base,
                    fillMultiplier, voidMultiplier);
    cudaDeviceSynchronize();
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilHeatHost_perf(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilN_Host");

    testStencilHeatHost_perf<blockEdgeSize, gridEdgeSize,
        SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize>>(i, base);
}

BOOST_AUTO_TEST_SUITE(performance)

BOOST_AUTO_TEST_SUITE(SparseGridGpu_test)

BOOST_AUTO_TEST_CASE(testStencilHeatHost_gridScaling)
{
   std::string testURI = suiteURI + ".host.stencil.dense.N.2D.gridScaling";
   unsigned int counter = 0;
   launch_testStencilHeatHost_perf<8, 128>(testURI, counter++);
   launch_testStencilHeatHost_perf<8, 256>(testURI, counter++);
   launch_testStencilHeatHost_perf<8, 512>(testURI, counter++);
//   launch_testStencilHeatHost_perf<8, 1024>(testURI, counter++);

   testSet.insert(testURI);
}
BOOST_AUTO_TEST_CASE(testStencilHeatHost_blockScaling)
{
   std::string testURI = suiteURI + ".host.stencil.dense.N.2D.blockScaling";
   unsigned int counter = 0;
   launch_testStencilHeatHost_perf<4, 2048>(testURI, counter++);
   launch_testStencilHeatHost_perf<8, 1024>(testURI, counter++);
   launch_testStencilHeatHost_perf<16, 512>(testURI, counter++);
//   launch_testStencilHeatHost_perf<32, 256>(testURI, counter++);

   testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatSparseHost_gridScaling)
{
   std::string testURI = suiteURI + ".host.stencil.sparse.N.2D.05.gridScaling";
   unsigned int counter = 0;
   constexpr unsigned int blockSize = 32;
   constexpr unsigned int __referenceBlockSize = 8;
   constexpr unsigned int __bsFactor = blockSize / __referenceBlockSize;
   launch_testStencilHeatSparseHost_perf<blockSize, 128/__bsFactor>(testURI, counter++, 1.45, 1, "05");
   launch_testStencilHeatSparseHost_perf<blockSize, 256/__bsFactor>(testURI, counter++, 1.45, 1, "05");
   launch_testStencilHeatSparseHost_perf<blockSize, 512/__bsFactor>(testURI, counter++, 1.45, 1, "05");
   launch_testStencilHeatSparseHost_perf<blockSize, 1024/__bsFactor>(testURI, counter++, 1.45, 1, "05");

   testSet.insert(testURI);
}
BOOST_AUTO_TEST_CASE(testStencilHeatSparseHost_blockScaling)
{
   std::string testURI = suiteURI + ".host.stencil.sparse.N.2D.05.blockScaling";
   unsigned int counter = 0;
   launch_testStencilHeatSparseHost_perf<4, 2048>(testURI, counter++, 1.45, 1, "05");
   launch_testStencilHeatSparseHost_perf<8, 1024>(testURI, counter++, 1.45, 1, "05");
   launch_testStencilHeatSparseHost_perf<16, 512>(testURI, counter++, 1.45, 1, "05");
   launch_testStencilHeatSparseHost_perf<32, 256>(testURI, counter++, 1.45, 1, "05");

   testSet.insert(testURI);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

