//
// Created by tommaso on 4/07/19.
//

#define SCAN_WITH_CUB
#define BOOST_TEST_DYN_LINK
#define DISABLE_MPI_WRITTERS

//#define SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE __launch_bounds__(512)
#define SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE_NO_SHARED __launch_bounds__(BLOCK_SIZE_STENCIL,12)

#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/SparseGridGpu.hpp"
#include "cuda_macro.h"
#include "util/stat/common_statistics.hpp"
#include "Plot/GoogleChart.hpp"
#include "util/performance/performance_util.hpp"
#include "SparseGridGpu/tests/utils/SparseGridGpu_testKernels.cuh"
#include <set>
#include "performancePlots.hpp"
#include "SparseGridGpu/tests/utils/SparseGridGpu_util_test.cuh"

extern char * test_dir;

// Property tree

report_sparse_grid_tests report_sparsegrid_funcs;
std::string suiteURI = "performance.SparseGridGpu";
std::set<std::string> testSet;


BOOST_AUTO_TEST_SUITE(performance)

BOOST_AUTO_TEST_SUITE(SparseGridGpu_test)


template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
void testStencilHeatGet_perf(unsigned int i, std::string base)
{
    auto testName = "In-place GET stencil";
    typedef HeatStencilGet<SparseGridZ::dims,0,1> Stencil01T;
    typedef HeatStencilGet<SparseGridZ::dims,1,0> Stencil10T;

    // typedef HeatStencilGet<SparseGridZ::dims,0,0> Stencil01T;
    // typedef HeatStencilGet<SparseGridZ::dims,0,0> Stencil10T;

    report_sparsegrid_funcs.graphs.put(base + ".dim",2);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*SparseGridZ::blockEdgeSize_);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*SparseGridZ::blockEdgeSize_);

    unsigned int iterations = 100;

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

    float gElemS = 2 * numElements / (1e9 * mean_tm);
    float gFlopsS = gElemS * Stencil01T::flops;
    std::cout << "Test: " << testName << std::endl;
    std::cout << "Block: " << SparseGridZ::blockEdgeSize_ << "x" << SparseGridZ::blockEdgeSize_ << std::endl;
    std::cout << "Grid: " << gridEdgeSize*SparseGridZ::blockEdgeSize_ << "x" << gridEdgeSize*SparseGridZ::blockEdgeSize_ << std::endl;
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
template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilHeatGet_perf(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    // typedef aggregate<float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilN");

    testStencilHeatGet_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize>>(i, base);
    cudaDeviceSynchronize();
}
template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilHeatGetZ_perf(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    // typedef aggregate<float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilZ");

    testStencilHeatGet_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu_z<dim, AggregateT, blockEdgeSize, chunkSize>>(i, base);
    cudaDeviceSynchronize();
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
void testStencilSkeleton_perf(unsigned int i, std::string base)
{
    auto testName = "In-place stencil";
    typedef SkeletonStencil<SparseGridZ::dims,0,1> Stencil01T;
    typedef SkeletonStencil<SparseGridZ::dims,1,0> Stencil10T;

    // typedef SkeletonStencil<SparseGridZ::dims,0,0> Stencil01T;
    // typedef SkeletonStencil<SparseGridZ::dims,0,0> Stencil10T;

    report_sparsegrid_funcs.graphs.put(base + ".dim",2);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*SparseGridZ::blockEdgeSize_);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*SparseGridZ::blockEdgeSize_);

    unsigned int iterations = 100;

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

    float gElemS = 2 * numElements / (1e9 * mean_tm);
    float gFlopsS = gElemS * Stencil01T::flops;
    std::cout << "Test: " << testName << std::endl;
    std::cout << "Block: " << SparseGridZ::blockEdgeSize_ << "x" << SparseGridZ::blockEdgeSize_ << std::endl;
    std::cout << "Grid: " << gridEdgeSize*SparseGridZ::blockEdgeSize_ << "x" << gridEdgeSize*SparseGridZ::blockEdgeSize_ << std::endl;
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

void launch_testConv3x3x3_perf_z_morton(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 3;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<8,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","Conv3x3x3");

    testConv3x3x3_perf<SparseGridGpu_z<dim, AggregateT, 8, chunkSize,long int>>("Convolution 3x3x3 Z-morton");
}

void launch_testConv3x3x3_perf(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 3;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<8,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","Conv3x3x3");

    testConv3x3x3_perf<SparseGridGpu<dim, AggregateT, 8, chunkSize,long int>>("Convolution 3x3x3 ");
}

void launch_testConv3x3x3_perf_no_shared_z_morton(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 3;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<8,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","Conv3x3x3");

    testConv3x3x3_no_shared_perf<SparseGridGpu_z<dim, AggregateT, 8, 512, long int>>("Convolution 3x3x3_noshared z-morton");
}

void launch_testConv3x3x3_perf_no_shared(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 3;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<8,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","Conv3x3x3");

    testConv3x3x3_no_shared_perf<SparseGridGpu<dim, AggregateT, 8, 512, long int>>("Convolution 3x3x3_noshared");
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilSkeleton_perf(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    // typedef aggregate<float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilN");

    testStencilSkeleton_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize>>(i, base);
    cudaDeviceSynchronize();
}


template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilSkeletonZ_perf(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    // typedef aggregate<float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilZ");

    testStencilSkeleton_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu_z<dim, AggregateT, blockEdgeSize, chunkSize>>(i, base);
    cudaDeviceSynchronize();
}

BOOST_AUTO_TEST_CASE(testConv3x3x3_noshared)
{
    std::string testURI = suiteURI + ".device.conv3x3x3_no_shared.sparse.N.3D.gridScaling";
    unsigned int counter = 0;
    launch_testConv3x3x3_perf_no_shared(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testConv3x3x3_noshared_z_morton)
{
    std::string testURI = suiteURI + ".device.conv3x3x3_no_shared.sparse.N.3D.gridScaling";
    unsigned int counter = 0;
    launch_testConv3x3x3_perf_no_shared_z_morton(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testConv3x3x3)
{
    std::string testURI = suiteURI + ".device.conv3x3x3.sparse.N.3D.gridScaling";
    unsigned int counter = 0;
    launch_testConv3x3x3_perf(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testConv3x3x3_zmorton)
{

    std::string testURI = suiteURI + ".device.conv3x3x3_zmorton.sparse.N.3D.gridScaling";
    unsigned int counter = 0;
    launch_testConv3x3x3_perf_z_morton(testURI, counter++);
    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilSkeleton_gridScaling)
{
    std::string testURI = suiteURI + ".device.skeleton.dense.N.2D.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilSkeleton_perf<blockEdgeSize, 128>(testURI, counter++);
    launch_testStencilSkeleton_perf<blockEdgeSize, 256>(testURI, counter++);
    launch_testStencilSkeleton_perf<blockEdgeSize, 512>(testURI, counter++);
    launch_testStencilSkeleton_perf<blockEdgeSize, 1024>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatGet_gridScaling_2)
{
   std::string testURI = suiteURI + ".device.stencilGet.dense.N.2D.2.gridScaling";
   unsigned int counter = 0;
   launch_testStencilHeatGet_perf<2, 512>(testURI, counter++);
   launch_testStencilHeatGet_perf<2, 1024>(testURI, counter++);
   launch_testStencilHeatGet_perf<2, 2048>(testURI, counter++);
   // launch_testStencilHeatGet_perf<2, 4096>(testURI, counter++); // test

   testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatGet_gridScaling_4)
{
   std::string testURI = suiteURI + ".device.stencilGet.dense.N.2D.4.gridScaling";
   unsigned int counter = 0;
   launch_testStencilHeatGet_perf<4, 256>(testURI, counter++);
   launch_testStencilHeatGet_perf<4, 512>(testURI, counter++);
   launch_testStencilHeatGet_perf<4, 1024>(testURI, counter++);
   launch_testStencilHeatGet_perf<4, 2048>(testURI, counter++);

   testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatGet_gridScaling_8)
{
   std::string testURI = suiteURI + ".device.stencilGet.dense.N.2D.8.gridScaling";
   unsigned int counter = 0;
   launch_testStencilHeatGet_perf<8, 128>(testURI, counter++);
   launch_testStencilHeatGet_perf<8, 256>(testURI, counter++);
   launch_testStencilHeatGet_perf<8, 512>(testURI, counter++);
   launch_testStencilHeatGet_perf<8, 1024>(testURI, counter++);

   testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatGet_gridScaling_16)
{
   std::string testURI = suiteURI + ".device.stencilGet.dense.N.2D.16.gridScaling";
   unsigned int counter = 0;
   launch_testStencilHeatGet_perf<16, 64>(testURI, counter++);
   launch_testStencilHeatGet_perf<16, 128>(testURI, counter++);
   launch_testStencilHeatGet_perf<16, 256>(testURI, counter++);
   launch_testStencilHeatGet_perf<16, 512>(testURI, counter++);

   testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatGet_gridScaling_32)
{
   std::string testURI = suiteURI + ".device.stencilGet.dense.N.2D.32.gridScaling";
   unsigned int counter = 0;
   launch_testStencilHeatGet_perf<32, 32>(testURI, counter++);
   launch_testStencilHeatGet_perf<32, 64>(testURI, counter++);
   launch_testStencilHeatGet_perf<32, 128>(testURI, counter++);
   launch_testStencilHeatGet_perf<32, 256>(testURI, counter++); // test

   testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatGet_blockScaling)
{
    std::string testURI = suiteURI + ".device.stencilGet.dense.N.2D.blockScaling";
    unsigned int counter = 0;
    // Note - blockEdgeSize == 2 doesn't work
    launch_testStencilHeatGet_perf<4, 2048>(testURI, counter++);
    launch_testStencilHeatGet_perf<8, 1024>(testURI, counter++);
    launch_testStencilHeatGet_perf<16, 512>(testURI, counter++);
    launch_testStencilHeatGet_perf<32, 256>(testURI, counter++);

    testSet.insert(testURI);
}


BOOST_AUTO_TEST_CASE(write_teport)
{
    write_test_report(report_sparsegrid_funcs, testSet);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
