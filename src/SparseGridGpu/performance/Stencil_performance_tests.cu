//
// Created by tommaso on 4/07/19.
//

#define BOOST_TEST_DYN_LINK
#define OPENFPM_DATA_ENABLE_IO_MODULE
#define DISABLE_MPI_WRITTERS

#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/SparseGridGpu.hpp"
#include "cuda_macro.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "util/stat/common_statistics.hpp"
#include "Plot/GoogleChart.hpp"
#include "util/performance/performance_util.hpp"
#include "SparseGridGpu/tests/utils/SparseGridGpu_testKernels.cuh"
#include <set>
#include "performancePlots.hpp"

extern char * test_dir;

// Property tree

report_sparse_grid_tests report_sparsegrid_funcs;
std::string suiteURI = "performance.SparseGridGpu";
std::set<std::string> testSet;

struct Fixture
{
    Fixture()
    {
        BOOST_TEST_MESSAGE( "Setup fixture" );
    }

    ~Fixture()
    {
        BOOST_TEST_MESSAGE( "Teardown fixture" );
        write_test_report(report_sparsegrid_funcs, testSet);
    }
};


template<unsigned int p, typename SparseGridType>
__global__ void insertValues2D(SparseGridType sparseGrid, const int offsetX=0, const int offsetY=0)
{
    sparseGrid.init();

    const auto bDimX = blockDim.x;
    const auto bDimY = blockDim.y;
    const auto bIdX = blockIdx.x;
    const auto bIdY = blockIdx.y;
    const auto tIdX = threadIdx.x;
    const auto tIdY = threadIdx.y;
    int x = bIdX * bDimX + tIdX + offsetX;
    int y = bIdY * bDimY + tIdY + offsetY;
    grid_key_dx<SparseGridType::d, int> coord({x, y});

    sparseGrid.template insert<p>(coord) = x*x*y*y; // some function...

    __syncthreads();

    sparseGrid.flush_block_insert();

    // Compiler avoid warning
    x++;
    y++;
}

template<unsigned int p, unsigned int chunksPerBlock, unsigned int blockEdgeSize, typename SparseGridType>
__global__ void insertValues2DBlocked(SparseGridType sparseGrid, const int sOffsetX=0, const int sOffsetY=0)
{
    constexpr unsigned int pMask = SparseGridType::pMask;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p> BlockT;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, pMask> MaskBlockT;

    sparseGrid.init();

    __shared__ BlockT *blocks[chunksPerBlock];
    __shared__ MaskBlockT *masks[chunksPerBlock];

    int posX = blockIdx.x * blockDim.x + threadIdx.x + sOffsetX;
    int posY = blockIdx.y * blockDim.y + threadIdx.y + sOffsetY;
    const unsigned int offsetX = posX % blockEdgeSize;
    const unsigned int offsetY = posY % blockEdgeSize;

    const unsigned int blockDimX = blockDim.x / blockEdgeSize;
    const unsigned int blockOffsetX = threadIdx.x / blockEdgeSize;
    const unsigned int blockOffsetY = threadIdx.y / blockEdgeSize;

    const unsigned int dataBlockNum = blockOffsetY*blockDimX + blockOffsetX;
    const unsigned int offset = offsetY * blockEdgeSize + offsetX;

//    if (offset == 0) // Just one thread per data block
//    {
        grid_key_dx<SparseGridType::d, int> blockCoord({posX / blockEdgeSize, posY / blockEdgeSize});
        auto encap = sparseGrid.insertBlockNew(sparseGrid.getBlockLinId(blockCoord));
        blocks[dataBlockNum] = &(encap.template get<p>());
        masks[dataBlockNum] = &(encap.template get<pMask>());
//    }

    __syncthreads();

    blocks[dataBlockNum]->block[offset] = posX*posX * posY*posY;
    BlockMapGpu_ker<>::setExist(masks[dataBlockNum]->block[offset]);

    __syncthreads();

    sparseGrid.flush_block_insert();
}



template<unsigned int p, unsigned int chunksPerBlock=1, typename SparseGridType, typename ScalarT>
__global__ void insertConstantValue(SparseGridType sparseGrid, ScalarT value)
{
    constexpr unsigned int pMask = SparseGridType::pMask;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p> BlockT;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, pMask> MaskBlockT;

    sparseGrid.init();

    __shared__ BlockT *blocks[chunksPerBlock];
    __shared__ MaskBlockT *masks[chunksPerBlock];

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    grid_key_dx<SparseGridType::d, size_t> coord({x, y, z});

    auto pos = sparseGrid.getLinId(coord);
    unsigned int dataBlockId = pos / BlockT::size;
    unsigned int offset = pos % BlockT::size;
    unsigned int dataBlockNum = dataBlockId % chunksPerBlock;

//    if (offset == 0) // Just one thread per data block
//    {
        auto encap = sparseGrid.insertBlockNew(dataBlockId);
        blocks[dataBlockNum] = &(encap.template get<p>());
        masks[dataBlockNum] = &(encap.template get<pMask>());
//    }

    __syncthreads();
    blocks[dataBlockNum]->block[offset] = value;
    BlockMapGpu_ker<>::setExist(masks[dataBlockNum]->block[offset]);

    __syncthreads();

    sparseGrid.flush_block_insert();

    // Compiler avoid warning
    x++;
    y++;
    z++;
}


template<unsigned int p, typename SparseGridType, typename ValueT>
__global__ void insertOneValue(SparseGridType sparseGrid, dim3 pt, ValueT value)
{
    sparseGrid.init();

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    dim3 thCoord(x, y, z);
    if (thCoord.x == pt.x && thCoord.y == pt.y && thCoord.z == pt.z)
    {
        grid_key_dx<SparseGridType::d, size_t> coord({x, y, z});
        sparseGrid.template insert<p>(coord) = value;
    }
    __syncthreads();

    sparseGrid.flush_block_insert();

    // Compiler avoid warning
    y++;
    z++;
}

template<unsigned int p, typename SparseGridType, typename VectorOutType>
__global__ void copyBlocksToOutput(SparseGridType sparseGrid, VectorOutType output)
{
    const auto bDimX = blockDim.x;
    const auto bDimY = blockDim.y;
    const auto bDimZ = blockDim.z;
    const auto bIdX = blockIdx.x;
    const auto bIdY = blockIdx.y;
    const auto bIdZ = blockIdx.z;
    const auto tIdX = threadIdx.x;
    const auto tIdY = threadIdx.y;
    const auto tIdZ = threadIdx.z;
    int x = bIdX * bDimX + tIdX;
    int y = bIdY * bDimY + tIdY;
    int z = bIdZ * bDimZ + tIdZ;
    grid_key_dx<SparseGridType::d, size_t> coord({x, y, z});

    size_t pos = sparseGrid.getLinId(coord);

    auto value = sparseGrid.template get<p>(coord);

    output.template get<p>(pos) = value;

    // Compiler avoid warning
    x++;
    y++;
    z++;
}

template<unsigned int dim, unsigned int p_src, unsigned int p_dst>
struct HeatStencil
{
	typedef NNStar stencil_type;

    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int flops = 3 + 2*dim;

    static constexpr unsigned int supportRadius = 1;

    /*! \brief Stencil function
     *
     * \param sparseGrid This is the sparse grid data-structure
     * \param dataBlockId The id of the block
     * \param offset index in local coordinate of the point where we are working
	 * \param dataBlockLoad dataBlock from where we read
	 * \param dataBlockStore dataBlock from where we write
	 * \param isActive the point is active if exist and is not padding
	 * \param dt delta t
     *
     *
     */
    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            const unsigned int dataBlockId,
            const openfpm::sparse_index<unsigned int> dataBlockIdPos,
            const unsigned int offset,
            const grid_key_dx<dim, int> & pointCoord,
            const DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore,
            bool isActive,
            float dt)
    {
        typedef typename SparseGridT::AggregateBlockType AggregateT;
        typedef ScalarTypeOf<AggregateT, p_src> ScalarT;

        constexpr unsigned int enlargedBlockSize = IntPow<
                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

        __shared__ ScalarT enlargedBlock[enlargedBlockSize];

        sparseGrid.loadGhostBlock<p_src>(dataBlockLoad, dataBlockIdPos, enlargedBlock);

        __syncthreads();

        if (isActive)
        {
            const auto coord = sparseGrid.getCoordInEnlargedBlock(offset);
            const auto linId = sparseGrid.getLinIdInEnlargedBlock(offset);
            ScalarT cur = enlargedBlock[linId];
            ScalarT laplacian = -2.0 * dim * cur; // The central part of the stencil

            for (int d = 0; d < dim; ++d)
            {
                auto nPlusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, 1);
                auto nMinusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, -1);
                ScalarT neighbourPlus = enlargedBlock[nPlusId];
                ScalarT neighbourMinus = enlargedBlock[nMinusId];
                laplacian += neighbourMinus + neighbourPlus;
            }
            enlargedBlock[linId] = cur + dt * laplacian;
        }
    }

    /*! \brief Stencil Host function
    *
    * \param sparseGrid This is the sparse grid data-structure
    * \param dataBlockId The id of the block
    * \param offset index in local coordinate of the point where we are working
    * \param dataBlockLoad dataBlock from where we read
    * \param dataBlockStore dataBlock from where we write
    * \param isActive the point is active if exist and is not padding
    * \param dt delta t
    *
    *
    */
    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __host__ void stencilHost(
            SparseGridT & sparseGrid,
            const unsigned int dataBlockId,
            const openfpm::sparse_index<unsigned int> dataBlockIdPos,
            const unsigned int offset,
            const grid_key_dx<dim, int> & pointCoord,
            const DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore,
            bool isActive,
            float dt)
    {
        constexpr unsigned int blockEdgeSize = SparseGridT::getBlockEdgeSize();

        if (isActive)
        {
            auto cur = dataBlockLoad.template get<p_src>()[offset];
            auto laplacian = -2.0 * dim * cur; // The central part of the stencil

            auto neighbourCoord = pointCoord;
            auto counter = offset;
            unsigned int dimStride = 1;
            for (int d = 0; d < dim; ++d)
            {
                const auto localOffset = counter % blockEdgeSize;

                if (localOffset == 0) // This means we are at the lower boundary for this dimension
                {
                    neighbourCoord.set_d(d, neighbourCoord.get(d) - 1);
                    laplacian += sparseGrid.template get<p_src>(neighbourCoord);
                    neighbourCoord.set_d(d, neighbourCoord.get(d) + 1);
                }
                else
                {
                    laplacian += dataBlockLoad.template get<p_src>()[offset - dimStride];
                }
                if (localOffset == blockEdgeSize - 1) // This means we are at the lower boundary for this dimension
                {
                neighbourCoord.set_d(d, neighbourCoord.get(d) + 1);
                laplacian += sparseGrid.template get<p_src>(neighbourCoord);
                neighbourCoord.set_d(d, neighbourCoord.get(d) - 1);
                }
                else
                {
                    laplacian += dataBlockLoad.template get<p_src>()[offset + dimStride];
                }
                //
                counter /= blockEdgeSize;
                dimStride *= blockEdgeSize;
            }
            dataBlockStore.template get<p_dst>()[offset] = cur + dt * laplacian;
        }
    }

    template <typename SparseGridT, typename CtxT>
    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
    {
        sparseGrid.template flush <sRight_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
void testStencilHeat_perf(unsigned int i, std::string base)
{
    auto testName = "In-place stencil";
    typedef HeatStencil<SparseGridZ::dims,0,1> StencilT;

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

	for (unsigned int iter=0; iter<iterations; ++iter)
	{
		cudaDeviceSynchronize();

		timer ts;
		ts.start();

		sparseGrid.template applyStencils<StencilT>(STENCIL_MODE_INPLACE, 0.1);

		cudaDeviceSynchronize();
		ts.stop();

		measures_tm.add(ts.getwct());

	    float gElemS = numElements / (1e9 * ts.getwct());
	    float gFlopsS = gElemS * StencilT::flops;

		measures_gf.add(gFlopsS);
	}

	double mean_tm = 0;
	double deviation_tm = 0;
	standard_deviation(measures_tm,mean_tm,deviation_tm);

	double mean_gf = 0;
	double deviation_gf = 0;
	standard_deviation(measures_gf,mean_gf,deviation_gf);

    // All times above are in ms

    float gElemS = numElements / (1e9 * mean_tm);
    float gFlopsS = gElemS * StencilT::flops;
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
void launch_testStencilHeat_perf(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilN");

    testStencilHeat_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize>>(i, base);
    cudaDeviceSynchronize();
}
template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilHeatZ_perf(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 2;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilZ");

    testStencilHeat_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu_z<dim, AggregateT, blockEdgeSize, chunkSize>>(i, base);
    cudaDeviceSynchronize();
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
void testStencilHeat3D_perf(unsigned int i, std::string base)
{
    auto testName = "In-place 3D stencil";
//    unsigned int gridEdgeSize = 128;
//    unsigned int gridEdgeSize = 64;
    typedef HeatStencil<SparseGridZ::dims,0,1> StencilT;

    report_sparsegrid_funcs.graphs.put(base + ".dim",3);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*SparseGridZ::blockEdgeSize_);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*SparseGridZ::blockEdgeSize_);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.z",gridEdgeSize*SparseGridZ::blockEdgeSize_);

    unsigned int iterations = 100;

    openfpm::vector<double> measures_gf;
    openfpm::vector<double> measures_tm;

    dim3 gridSize(gridEdgeSize, gridEdgeSize, gridEdgeSize);
    dim3 blockSize(SparseGridZ::blockEdgeSize_, SparseGridZ::blockEdgeSize_, SparseGridZ::blockEdgeSize_);

    typename SparseGridZ::grid_info blockGeometry(gridSize);
    SparseGridZ sparseGrid(blockGeometry);
    mgpu::ofp_context_t ctx;
    sparseGrid.template setBackgroundValue<0>(0);

    unsigned long long numElements = gridEdgeSize*SparseGridZ::blockEdgeSize_
            *gridEdgeSize*SparseGridZ::blockEdgeSize_
            *gridEdgeSize*SparseGridZ::blockEdgeSize_;

    // Initialize the grid
    sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
    CUDA_LAUNCH_DIM3((insertConstantValue<0>),gridSize, blockSize,sparseGrid.toKernel(), 0);
    sparseGrid.template flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
    dim3 sourcePt(gridSize.x * SparseGridZ::blockEdgeSize_ / 2,
            gridSize.y * SparseGridZ::blockEdgeSize_ / 2,
            gridSize.z * SparseGridZ::blockEdgeSize_ / 2);
    insertOneValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), sourcePt, 100);
    sparseGrid.template flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

    for (unsigned int iter=0; iter<iterations; ++iter)
    {
        cudaDeviceSynchronize();

        timer ts;
        ts.start();

        sparseGrid.template applyStencils<StencilT>(STENCIL_MODE_INPLACE, 0.1);

        cudaDeviceSynchronize();
        ts.stop();

        measures_tm.add(ts.getwct());

        float gElemS = numElements / (1e9 * ts.getwct());
        float gFlopsS = gElemS * StencilT::flops;

        measures_gf.add(gFlopsS);
    }

    double mean_tm = 0;
    double deviation_tm = 0;
    standard_deviation(measures_tm,mean_tm,deviation_tm);

    double mean_gf = 0;
    double deviation_gf = 0;
    standard_deviation(measures_gf,mean_gf,deviation_gf);

    // All times above are in ms

    float gElemS = numElements / (1e9 * mean_tm);
    float gFlopsS = gElemS * StencilT::flops;
    std::cout << "Test: " << testName << std::endl;
    std::cout << "Block: " << SparseGridZ::blockEdgeSize_
              << "x" << SparseGridZ::blockEdgeSize_
              << "x" << SparseGridZ::blockEdgeSize_
              << std::endl;
    std::cout << "Grid: " << gridEdgeSize*SparseGridZ::blockEdgeSize_
        << "x" << gridEdgeSize*SparseGridZ::blockEdgeSize_
        << "x" << gridEdgeSize*SparseGridZ::blockEdgeSize_
        << std::endl;
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
void launch_testStencilHeat3D_perf(std::string testURI, unsigned int i)
{
    constexpr unsigned int dim = 3;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilN3D");

    testStencilHeat3D_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize>>(i, base);
    cudaDeviceSynchronize();
}

//template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
//void testStencilHeatSparse_perf(unsigned int i, std::string base)
//{
//    auto testName = "In-place sparse stencil";
////    unsigned int gridEdgeSize = 128;
//    constexpr unsigned int dim = SparseGridZ::dims;
////    const unsigned int blockEdgeSize = SparseGridZ::blockEdgeSize_;
//
//    typedef HeatStencil<dim, 0, 1> Stencil01T;
//    typedef HeatStencil<dim, 1, 0> Stencil10T;
//
////    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").stencil");
//
//    report_sparsegrid_funcs.graphs.put(base + ".dim",2);
//    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
//    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*blockEdgeSize);
//    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*blockEdgeSize);
//
//    unsigned int iterations = 100;
//
//    openfpm::vector<double> measures_gf;
//    openfpm::vector<double> measures_tm;
//
//    dim3 gridSize(gridEdgeSize, gridEdgeSize);
//    dim3 blockSize(blockEdgeSize,blockEdgeSize);
//    size_t sz[2] = {1000000,1000000};
//    typename SparseGridZ::grid_info blockGeometry(sz);
//    SparseGridZ sparseGrid(blockGeometry);
//    mgpu::ofp_context_t ctx;
//    sparseGrid.template setBackgroundValue<0>(0);
//
//    ///// Insert sparse content, a set of 3 hollow spheres /////
//    constexpr unsigned int rBig = gridEdgeSize * blockEdgeSize / 2;
//    constexpr unsigned int rSmall = rBig/2;
//    constexpr unsigned int rBig2 = rBig;
//    constexpr unsigned int rSmall2 = rBig2 - (rBig2/16);
//    constexpr unsigned int rBig3 = rBig/8;
//    constexpr unsigned int rSmall3 = rBig3 - (rBig3/10);
//    // Sphere 1
//    grid_key_dx<2,int> start1({500000,500000});
//    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
//    CUDA_LAUNCH_DIM3((insertSphere<0>),
//                     gridSize, dim3(blockEdgeSize*blockEdgeSize,1,1),
//                     sparseGrid.toKernel(), start1, rBig, rSmall, 1);
//    cudaDeviceSynchronize();
//    sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//    cudaDeviceSynchronize();
//
//    // Sphere 2
//    grid_key_dx<2,int> start2({500000+rBig,500000+rBig});
//    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
//    CUDA_LAUNCH_DIM3((insertSphere<0>),
//                     gridSize, dim3(blockEdgeSize*blockEdgeSize,1,1),
//                     sparseGrid.toKernel(), start2, rBig2, rSmall2, 1);
//    cudaDeviceSynchronize();
//    sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//    cudaDeviceSynchronize();
//
//    // Sphere 3
//    grid_key_dx<2,int> start3({500000+rBig,500000});
//    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
//    CUDA_LAUNCH_DIM3((insertSphere<0>),
//                     gridSize, dim3(blockEdgeSize*blockEdgeSize,1,1),
//                     sparseGrid.toKernel(), start3, rBig3, rSmall3, 1);
//    cudaDeviceSynchronize();
//    sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//    cudaDeviceSynchronize();
//    ///// /////
//
//    sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
//    sparseGrid.tagBoundaries();
//
//    sparseGrid.template deviceToHost<0>(); // NECESSARY as count takes place on Host!
//    auto existingElements = sparseGrid.countExistingElements();
//    auto boundaryElements = sparseGrid.countBoundaryElements();
//    unsigned long long numElements = existingElements - boundaryElements;
//
//    // Now apply some boundary conditions
//    sparseGrid.template applyBoundaryStencils<BoundaryStencilSetXRescaled<dim,0,0>>(
//            500000, 500000+(2*rBig),
//            0.0, 10.0);
//    cudaDeviceSynchronize();
//
//    for (unsigned int iter=0; iter<iterations; ++iter)
//    {
//        cudaDeviceSynchronize();
//
//        timer ts;
//        ts.start();
//
//        sparseGrid.template applyStencils<Stencil01T>(STENCIL_MODE_INPLACE, 0.1);
//        cudaDeviceSynchronize();
//        sparseGrid.template applyStencils<Stencil10T>(STENCIL_MODE_INPLACE, 0.1);
//        cudaDeviceSynchronize();
//
//        ts.stop();
//
//        measures_tm.add(ts.getwct());
//
//        float gElemS = numElements / (1e9 * ts.getwct());
//        float gFlopsS = gElemS * Stencil01T::flops;
//
//        measures_gf.add(gFlopsS);
//    }
//
//    double mean_tm = 0;
//    double deviation_tm = 0;
//    standard_deviation(measures_tm,mean_tm,deviation_tm);
//
//    double mean_gf = 0;
//    double deviation_gf = 0;
//    standard_deviation(measures_gf,mean_gf,deviation_gf);
//
//    // All times above are in ms
//
//    float gElemS = numElements / (1e9 * mean_tm);
//    float gFlopsS = gElemS * Stencil01T::flops;
//    std::cout << "Test: " << testName << std::endl;
//    std::cout << "Block: " << blockEdgeSize << "x" << blockEdgeSize << std::endl;
//    std::cout << "Grid: " << gridEdgeSize*blockEdgeSize << "x" << gridEdgeSize*blockEdgeSize << std::endl;
//    double dataOccupancyMean, dataOccupancyDev;
//    sparseGrid.deviceToHost();
//    sparseGrid.measureBlockOccupancy(dataOccupancyMean, dataOccupancyDev);std::cout << "Data Occupancy: " << dataOccupancyMean << " dev:" << dataOccupancyDev << std::endl;
//    report_sparsegrid_funcs.graphs.put(base + ".dataOccupancy.mean",dataOccupancyMean);
//    report_sparsegrid_funcs.graphs.put(base +".dataOccupancy.dev",dataOccupancyDev);
//    std::cout << "Iterations: " << iterations << std::endl;
//    std::cout << "\tStencil: " << mean_gf << " dev:" << deviation_gf << " s" << std::endl;
//    std::cout << "Throughput: " << std::endl << "\t " << gElemS << " GElem/s " << std::endl << "\t " << gFlopsS << " GFlops/s" << std::endl;
//
//    report_sparsegrid_funcs.graphs.put(base + ".GFlops.mean",mean_gf);
//    report_sparsegrid_funcs.graphs.put(base +".GFlops.dev",deviation_gf);
//    report_sparsegrid_funcs.graphs.put(base + ".time.mean",mean_tm);
//    report_sparsegrid_funcs.graphs.put(base +".time.dev",deviation_tm);
//
////    // DEBUG
////    sparseGrid.template deviceToHost<0,1>();
////    sparseGrid.write("SparseGridGPU_testStencilHeatSparse_perf_DEBUG.vtk");
//}
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

        float gElemS = numElements / (1e9 * ts.getwct());
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

    float gElemS = numElements / (1e9 * mean_tm);
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

//    // DEBUG
//    sparseGrid.template deviceToHost<0,1>();
//    sparseGrid.write("SparseGridGPU_testStencilHeatSparse_perf_DEBUG.vtk");
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

//template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
//void testStencilHeat3DSparse_perf(unsigned int i, std::string base)
//{
//    auto testName = "In-place 3D sparse stencil";
////    unsigned int gridEdgeSize = 32;
//    constexpr unsigned int dim = SparseGridZ::dims;
////    const unsigned int blockEdgeSize = SparseGridZ::blockEdgeSize_;
//
//    typedef HeatStencil<dim, 0, 1> Stencil01T;
//    typedef HeatStencil<dim, 1, 0> Stencil10T;
//
////    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").stencil");
//
//    report_sparsegrid_funcs.graphs.put(base + ".dim",dim);
//    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
//    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x", gridEdgeSize * blockEdgeSize);
//    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y", gridEdgeSize * blockEdgeSize);
//    report_sparsegrid_funcs.graphs.put(base + ".gridSize.z", gridEdgeSize * blockEdgeSize);
//
//    unsigned int iterations = 100;
//
//    openfpm::vector<double> measures_gf;
//    openfpm::vector<double> measures_tm;
//
//    dim3 gridSize(gridEdgeSize, gridEdgeSize, gridEdgeSize);
//    dim3 blockSize(blockEdgeSize, blockEdgeSize, blockEdgeSize);
//    unsigned int spatialEdgeSize = 10000;
//    size_t sz[3] = {spatialEdgeSize, spatialEdgeSize, spatialEdgeSize};
//    typename SparseGridZ::grid_info blockGeometry(sz);
//    SparseGridZ sparseGrid(blockGeometry);
//    mgpu::ofp_context_t ctx;
//    sparseGrid.template setBackgroundValue<0>(0);
//
//    ///// Insert sparse content, a set of 3 hollow spheres /////
//    constexpr unsigned int rBig = gridEdgeSize * blockEdgeSize / 2;
//    constexpr unsigned int rSmall = rBig/2;
//    constexpr unsigned int rBig2 = rBig;
//    constexpr unsigned int rSmall2 = rBig2 - (rBig2/3);
//    constexpr unsigned int rBig3 = rBig/4;
//    constexpr unsigned int rSmall3 = rBig3 - (rBig3/4);
//    // Sphere 1
//    unsigned int centerPoint = spatialEdgeSize/2;
//    grid_key_dx<dim,int> start1({centerPoint, centerPoint, centerPoint});
//    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
//    CUDA_LAUNCH_DIM3((insertSphere3D<0>),
//                     gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
//                     sparseGrid.toKernel(), start1, rBig, rSmall, 1);
//    cudaDeviceSynchronize();
//    sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//    cudaDeviceSynchronize();
//
//    sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
//    sparseGrid.tagBoundaries();
//
//    // Sphere 2
//    grid_key_dx<dim,int> start2({centerPoint - rBig, centerPoint - rBig, centerPoint - rBig});
//    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
//    CUDA_LAUNCH_DIM3((insertSphere3D<0>),
//                     gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
//                     sparseGrid.toKernel(), start2, rBig2, rSmall2, 1);
//    cudaDeviceSynchronize();
//    sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//    cudaDeviceSynchronize();
//
//    // Sphere 3
//    grid_key_dx<dim,int> start3({centerPoint + rBig, centerPoint - rBig, centerPoint - rBig});
//    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
//    CUDA_LAUNCH_DIM3((insertSphere3D<0>),
//                     gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
//                     sparseGrid.toKernel(), start3, rBig3, rSmall3, 1);
//    cudaDeviceSynchronize();
//    sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//    cudaDeviceSynchronize();
//    ///// /////
//
//    sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
//    sparseGrid.tagBoundaries();
//
//    sparseGrid.template deviceToHost<0>(); // NECESSARY as count takes place on Host!
//    auto existingElements = sparseGrid.countExistingElements();
//    auto boundaryElements = sparseGrid.countBoundaryElements();
//    unsigned long long numElements = existingElements - boundaryElements;
//
//    // Now apply some boundary conditions
//    sparseGrid.template applyBoundaryStencils<BoundaryStencilSetXRescaled<dim,0,0>>(
//            centerPoint - rBig, centerPoint + (2 * rBig),
//            0.0, 10.0);
//    cudaDeviceSynchronize();
//
//    for (unsigned int iter=0; iter<iterations; ++iter)
//    {
//        cudaDeviceSynchronize();
//
//        timer ts;
//        ts.start();
//
//        sparseGrid.template applyStencils<Stencil01T>(STENCIL_MODE_INPLACE, 0.1);
//        cudaDeviceSynchronize();
//        sparseGrid.template applyStencils<Stencil10T>(STENCIL_MODE_INPLACE, 0.1);
//        cudaDeviceSynchronize();
//
//        ts.stop();
//
//        measures_tm.add(ts.getwct());
//
//        float gElemS = numElements / (1e9 * ts.getwct());
//        float gFlopsS = gElemS * Stencil01T::flops;
//
//        measures_gf.add(gFlopsS);
//    }
//
//    double mean_tm = 0;
//    double deviation_tm = 0;
//    standard_deviation(measures_tm,mean_tm,deviation_tm);
//
//    double mean_gf = 0;
//    double deviation_gf = 0;
//    standard_deviation(measures_gf,mean_gf,deviation_gf);
//
//    // All times above are in ms
//
//    float gElemS = numElements / (1e9 * mean_tm);
//    float gFlopsS = gElemS * Stencil01T::flops;
//    std::cout << "Test: " << testName << std::endl;
//    std::cout << "Block: " << blockEdgeSize << "x" << blockEdgeSize << "x" << blockEdgeSize << std::endl;
//    std::cout << "Grid: " << gridEdgeSize * blockEdgeSize
//              << "x" << gridEdgeSize * blockEdgeSize
//              << "x" << gridEdgeSize * blockEdgeSize
//              << std::endl;
//    double dataOccupancyMean, dataOccupancyDev;
//    sparseGrid.deviceToHost();
//    sparseGrid.measureBlockOccupancy(dataOccupancyMean, dataOccupancyDev);std::cout << "Data Occupancy: " << dataOccupancyMean << " dev:" << dataOccupancyDev << std::endl;
//    report_sparsegrid_funcs.graphs.put(base + ".dataOccupancy.mean",dataOccupancyMean);
//    report_sparsegrid_funcs.graphs.put(base +".dataOccupancy.dev",dataOccupancyDev);
//    std::cout << "Iterations: " << iterations << std::endl;
//    std::cout << "\tStencil: " << mean_gf << " dev:" << deviation_gf << " s" << std::endl;
//    std::cout << "Throughput: " << std::endl << "\t " << gElemS << " GElem/s " << std::endl << "\t " << gFlopsS << " GFlops/s" << std::endl;
//
//    report_sparsegrid_funcs.graphs.put(base + ".GFlops.mean",mean_gf);
//    report_sparsegrid_funcs.graphs.put(base +".GFlops.dev",deviation_gf);
//    report_sparsegrid_funcs.graphs.put(base + ".time.mean",mean_tm);
//    report_sparsegrid_funcs.graphs.put(base +".time.dev",deviation_tm);
//
////    // DEBUG
////    sparseGrid.template deviceToHost<0,1>();
////    sparseGrid.write("SparseGridGPU_testStencilHeat3DSparse_perf_DEBUG.vtk");
//}
template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
void testStencilHeat3DSparse_perf(unsigned int i, std::string base, float fillMultiplier=1, float voidMultiplier=1)
{
    auto testName = "In-place 3D sparse stencil";
//    unsigned int gridEdgeSize = 32;
    constexpr unsigned int dim = SparseGridZ::dims;
//    const unsigned int blockEdgeSize = SparseGridZ::blockEdgeSize_;

    typedef HeatStencil<dim, 0, 1> Stencil01T;
    typedef HeatStencil<dim, 1, 0> Stencil10T;

//    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").stencil");

    report_sparsegrid_funcs.graphs.put(base + ".dim",dim);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x", gridEdgeSize * blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y", gridEdgeSize * blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.z", gridEdgeSize * blockEdgeSize);

    unsigned int iterations = 100;

    openfpm::vector<double> measures_gf;
    openfpm::vector<double> measures_tm;

    dim3 gridSize(gridEdgeSize, gridEdgeSize, gridEdgeSize);
    dim3 blockSize(blockEdgeSize, blockEdgeSize, blockEdgeSize);
    unsigned int spatialEdgeSize = 10000;
    size_t sz[3] = {spatialEdgeSize, spatialEdgeSize, spatialEdgeSize};
    typename SparseGridZ::grid_info blockGeometry(sz);
    SparseGridZ sparseGrid(blockGeometry);
    mgpu::ofp_context_t ctx;
    sparseGrid.template setBackgroundValue<0>(0);

    ///// Insert sparse content, a set of concentric spheres /////
    float allMultiplier = fillMultiplier + voidMultiplier;
    const unsigned int numSpheres = gridEdgeSize / (2*allMultiplier);
    unsigned int centerPoint = spatialEdgeSize / 2;

    for (int i = 1; i <= numSpheres; ++i)
    {
        unsigned int rBig = allMultiplier*i * blockEdgeSize;
        unsigned int rSmall = (allMultiplier*i - fillMultiplier) * blockEdgeSize;
        // Sphere i-th
        grid_key_dx<dim, int> start1({centerPoint, centerPoint, centerPoint});
        sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
        CUDA_LAUNCH_DIM3((insertSphere3D<0>),
                         gridSize, dim3(blockEdgeSize * blockEdgeSize * blockEdgeSize, 1, 1),
                         sparseGrid.toKernel(), start1, rBig, rSmall, 1);
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

        float gElemS = numElements / (1e9 * ts.getwct());
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

    float gElemS = numElements / (1e9 * mean_tm);
    float gFlopsS = gElemS * Stencil01T::flops;
    std::cout << "Test: " << testName << std::endl;
    std::cout << "Block: " << blockEdgeSize << "x" << blockEdgeSize << "x" << blockEdgeSize << std::endl;
    std::cout << "Grid: " << gridEdgeSize * blockEdgeSize
              << "x" << gridEdgeSize * blockEdgeSize
              << "x" << gridEdgeSize * blockEdgeSize
              << std::endl;
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

//    // DEBUG
//    sparseGrid.template deviceToHost<0,1>();
//    sparseGrid.write("SparseGridGPU_testStencilHeat3DSparse_perf_DEBUG.vtk");
}
template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void launch_testStencilHeat3DSparse_perf(std::string testURI, unsigned int i,
        float fillMultiplier=1, float voidMultiplier=1)
{
    constexpr unsigned int dim = 3;
    typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","StencilN3DSparse");

    testStencilHeat3DSparse_perf<blockEdgeSize, gridEdgeSize,
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize, long int>>(i, base,
                    fillMultiplier, voidMultiplier);
    cudaDeviceSynchronize();
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize, typename SparseGridZ>
void testStencilHeatHost_perf(unsigned int i, std::string base)
{
    // todo: Make sure to reimplement the host stencil application function to pre-load to a block of memory both content and ghost
    // this way we can avoid binary searches...
    auto testName = "In-place stencil HOST";
    typedef HeatStencil<SparseGridZ::dims,0,1> StencilT;

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

    sparseGrid.template deviceToHost<0>();

    for (unsigned int iter=0; iter<iterations; ++iter)
    {
        cudaDeviceSynchronize();

        timer ts;
        ts.start();

        sparseGrid.template applyStencilsHost<StencilT>(STENCIL_MODE_INPLACE, 0.1);

        cudaDeviceSynchronize();
        ts.stop();

        measures_tm.add(ts.getwct());

        float gElemS = numElements / (1e9 * ts.getwct());
        float gFlopsS = gElemS * StencilT::flops;

        measures_gf.add(gFlopsS);
    }

    double mean_tm = 0;
    double deviation_tm = 0;
    standard_deviation(measures_tm,mean_tm,deviation_tm);

    double mean_gf = 0;
    double deviation_gf = 0;
    standard_deviation(measures_gf,mean_gf,deviation_gf);

    // All times above are in ms

    float gElemS = numElements / (1e9 * mean_tm);
    float gFlopsS = gElemS * StencilT::flops;

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

BOOST_AUTO_TEST_SUITE(performance, *boost::unit_test::fixture<Fixture>())

BOOST_AUTO_TEST_SUITE(SparseGridGpu_test)

//BOOST_AUTO_TEST_CASE(testStencilHeatHost_gridScaling)
//{
//    std::string testURI = suiteURI + ".host.stencil.dense.N.2D.gridScaling";
//    unsigned int counter = 0;
//    launch_testStencilHeatHost_perf<8, 128>(testURI, counter++);
//    launch_testStencilHeatHost_perf<8, 256>(testURI, counter++);
//    launch_testStencilHeatHost_perf<8, 512>(testURI, counter++);
//    launch_testStencilHeatHost_perf<8, 1024>(testURI, counter++);
////    launch_testStencilHeatHost_perf<8, 2048>(testURI, counter++);

//    testSet.insert(testURI);
//}
//
//BOOST_AUTO_TEST_CASE(testStencilHeatHost_blockScaling)
//{
//    std::string testURI = suiteURI + ".host.stencil.dense.N.2D.blockScaling";
//    unsigned int counter = 0;
//    launch_testStencilHeatHost_perf<4, 2048>(testURI, counter++);
//    launch_testStencilHeatHost_perf<8, 1024>(testURI, counter++);
//    launch_testStencilHeatHost_perf<16, 512>(testURI, counter++);
//    launch_testStencilHeatHost_perf<32, 256>(testURI, counter++);

//    testSet.insert(testURI);
//}

BOOST_AUTO_TEST_CASE(testStencilHeat_gridScaling)
//BOOST_AUTO_TEST_CASE(testStencilHeatSparse10_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.dense.N.2D.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeat_perf<blockEdgeSize, 128>(testURI, counter++);
    launch_testStencilHeat_perf<blockEdgeSize, 256>(testURI, counter++);
    launch_testStencilHeat_perf<blockEdgeSize, 512>(testURI, counter++);
    launch_testStencilHeat_perf<blockEdgeSize, 1024>(testURI, counter++);
//    launch_testStencilHeat_perf<blockEdgeSize, 2048>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeat_blockScaling)
{
    std::string testURI = suiteURI + ".device.stencil.dense.N.2D.blockScaling";
    unsigned int counter = 0;
    // Note - blockEdgeSize == 2 doesn't work
    launch_testStencilHeat_perf<4, 2048>(testURI, counter++);
    launch_testStencilHeat_perf<8, 1024>(testURI, counter++);
    launch_testStencilHeat_perf<16, 512>(testURI, counter++);
    launch_testStencilHeat_perf<32, 256>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatZ_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.dense.Z.2D.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeatZ_perf<blockEdgeSize, 128>(testURI, counter++);
    launch_testStencilHeatZ_perf<blockEdgeSize, 256>(testURI, counter++);
    launch_testStencilHeatZ_perf<blockEdgeSize, 512>(testURI, counter++);
    launch_testStencilHeatZ_perf<blockEdgeSize, 1024>(testURI, counter++);
//    launch_testStencilHeatZ_perf<blockEdgeSize, 2048>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatZ_blockScaling)
{
    std::string testURI = suiteURI + ".device.stencil.dense.Z.2D.blockScaling";
    unsigned int counter = 0;
    // Note - blockEdgeSize == 2 doesn't work
    launch_testStencilHeatZ_perf<4, 2048>(testURI, counter++);
    launch_testStencilHeatZ_perf<8, 1024>(testURI, counter++);
    launch_testStencilHeatZ_perf<16, 512>(testURI, counter++);
    launch_testStencilHeatZ_perf<32, 256>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeat3D_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.dense.N.3D.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeat3D_perf<blockEdgeSize, 8>(testURI, counter++);
    launch_testStencilHeat3D_perf<blockEdgeSize, 16>(testURI, counter++);
    launch_testStencilHeat3D_perf<blockEdgeSize, 32>(testURI, counter++);
    launch_testStencilHeat3D_perf<blockEdgeSize, 64>(testURI, counter++);
//    launch_testStencilHeat3D_perf<blockEdgeSize, 128>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeat3D_blockScaling)
{
    std::string testURI = suiteURI + ".device.stencil.dense.N.3D.blockScaling";
    unsigned int counter = 0;
    launch_testStencilHeat3D_perf<2, 128>(testURI, counter++);
    launch_testStencilHeat3D_perf<4, 64>(testURI, counter++);
    launch_testStencilHeat3D_perf<8, 32>(testURI, counter++);
//    launch_testStencilHeat3D_perf<16, 16>(testURI, counter++); // Too big, it doesn't work

    testSet.insert(testURI);
}

//BOOST_AUTO_TEST_CASE(testStencilHeatZ3D)
//{
//    constexpr unsigned int dim = 3;
//    constexpr unsigned int blockEdgeSize = 4;
//
//    typedef aggregate<float,float> AggregateT;
//    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
//
//    report_sparsegrid_funcs.graphs.put("performance.SparseGridGpu(1).stencil.test.name","StencilZ3D");
//
//    testStencilHeat3D_perf<SparseGridGpu_z<dim, AggregateT, blockEdgeSize, chunkSize>>(1);
//}

BOOST_AUTO_TEST_CASE(testStencilHeatSparse05_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.05.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeatSparse_perf<blockEdgeSize, 128>(testURI, counter++, 1, 1, "05");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 256>(testURI, counter++, 1, 1, "05");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 512>(testURI, counter++, 1, 1, "05");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 1024>(testURI, counter++, 1, 1, "05");
////    launch_testStencilHeatSparse_perf<blockEdgeSize, 2048>(testURI, counter++, 1, 1, "05);

    testSet.insert(testURI);
}
BOOST_AUTO_TEST_CASE(testStencilHeatSparse05_blockScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.05.blockScaling";
    unsigned int counter = 0;
    // Note - blockEdgeSize == 2 doesn't work
    launch_testStencilHeatSparse_perf<2, 2048>(testURI, counter++, 1, 1, "05");
    launch_testStencilHeatSparse_perf<4, 1024>(testURI, counter++, 1, 1, "05");
    launch_testStencilHeatSparse_perf<8, 512>(testURI, counter++, 1, 1, "05");
    launch_testStencilHeatSparse_perf<16, 256>(testURI, counter++, 1, 1, "05");
    launch_testStencilHeatSparse_perf<32, 128>(testURI, counter++, 1, 1, "05");

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatSparse08_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.08.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeatSparse_perf<blockEdgeSize, 128>(testURI, counter++, 1, 0.25, "08");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 256>(testURI, counter++, 1, 0.25, "08");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 512>(testURI, counter++, 1, 0.25, "08");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 1024>(testURI, counter++, 1, 0.25, "08");
////    launch_testStencilHeatSparse_perf<blockEdgeSize, 2048>(testURI, counter++, 1, 0.25, "08");

    testSet.insert(testURI);
}
BOOST_AUTO_TEST_CASE(testStencilHeatSparse08_blockScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.08.blockScaling";
    unsigned int counter = 0;
    launch_testStencilHeatSparse_perf<2, 2048>(testURI, counter++, 1, 0.25, "08");
    launch_testStencilHeatSparse_perf<4, 1024>(testURI, counter++, 1, 0.25, "08");
    launch_testStencilHeatSparse_perf<8, 512>(testURI, counter++, 1, 0.25, "08");
    launch_testStencilHeatSparse_perf<16, 256>(testURI, counter++, 1, 0.25, "08");
    launch_testStencilHeatSparse_perf<32, 128>(testURI, counter++, 1, 0.25, "08");

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeatSparse09_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.09.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeatSparse_perf<blockEdgeSize, 128>(testURI, counter++, 1, 0.1, "09");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 256>(testURI, counter++, 1, 0.1, "09");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 512>(testURI, counter++, 1, 0.1, "09");
    launch_testStencilHeatSparse_perf<blockEdgeSize, 1024>(testURI, counter++, 1, 0.1, "09");
////    launch_testStencilHeatSparse_perf<blockEdgeSize, 2048>(testURI, counter++, 1, 0.1, "09");

    testSet.insert(testURI);
}
BOOST_AUTO_TEST_CASE(testStencilHeatSparse09_blockScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.2D.09.blockScaling";
    unsigned int counter = 0;
    // Note - blockEdgeSize == 2 doesn't work
    launch_testStencilHeatSparse_perf<2, 2048>(testURI, counter++, 1, 0.1, "09");
    launch_testStencilHeatSparse_perf<4, 1024>(testURI, counter++, 1, 0.1, "09");
    launch_testStencilHeatSparse_perf<8, 512>(testURI, counter++, 1, 0.1, "09");
    launch_testStencilHeatSparse_perf<16, 256>(testURI, counter++, 1, 0.1, "09");
    launch_testStencilHeatSparse_perf<32, 128>(testURI, counter++, 1, 0.1, "09");

    testSet.insert(testURI);
}

        BOOST_AUTO_TEST_CASE(testStencilHeatSparseZ05_gridScaling)
        {
            std::string testURI = suiteURI + ".device.stencil.sparse.Z.2D.05.gridScaling";
            unsigned int counter = 0;
            constexpr unsigned int blockEdgeSize = 8;
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 128>(testURI, counter++, 1, 1, "05");
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 256>(testURI, counter++, 1, 1, "05");
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 512>(testURI, counter++, 1, 1, "05");
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 1024>(testURI, counter++, 1, 1, "05");
////    launch_testStencilHeatSparseZ_perf<blockEdgeSize, 2048>(testURI, counter++, 1, 1, "05);

            testSet.insert(testURI);
        }
        BOOST_AUTO_TEST_CASE(testStencilHeatSparseZ05_blockScaling)
        {
            std::string testURI = suiteURI + ".device.stencil.sparse.Z.2D.05.blockScaling";
            unsigned int counter = 0;
            // Note - blockEdgeSize == 2 doesn't work
            launch_testStencilHeatSparseZ_perf<2, 2048>(testURI, counter++, 1, 1, "05");
            launch_testStencilHeatSparseZ_perf<4, 1024>(testURI, counter++, 1, 1, "05");
            launch_testStencilHeatSparseZ_perf<8, 512>(testURI, counter++, 1, 1, "05");
            launch_testStencilHeatSparseZ_perf<16, 256>(testURI, counter++, 1, 1, "05");
            launch_testStencilHeatSparseZ_perf<32, 128>(testURI, counter++, 1, 1, "05");

            testSet.insert(testURI);
        }

        BOOST_AUTO_TEST_CASE(testStencilHeatSparseZ08_gridScaling)
        {
            std::string testURI = suiteURI + ".device.stencil.sparse.Z.2D.08.gridScaling";
            unsigned int counter = 0;
            constexpr unsigned int blockEdgeSize = 8;
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 128>(testURI, counter++, 1, 0.25, "08");
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 256>(testURI, counter++, 1, 0.25, "08");
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 512>(testURI, counter++, 1, 0.25, "08");
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 1024>(testURI, counter++, 1, 0.25, "08");
////    launch_testStencilHeatSparseZ_perf<blockEdgeSize, 2048>(testURI, counter++, 1, 0.25, "08");

            testSet.insert(testURI);
        }
        BOOST_AUTO_TEST_CASE(testStencilHeatSparseZ08_blockScaling)
        {
            std::string testURI = suiteURI + ".device.stencil.sparse.Z.2D.08.blockScaling";
            unsigned int counter = 0;
            launch_testStencilHeatSparseZ_perf<2, 2048>(testURI, counter++, 1, 0.25, "08");
            launch_testStencilHeatSparseZ_perf<4, 1024>(testURI, counter++, 1, 0.25, "08");
            launch_testStencilHeatSparseZ_perf<8, 512>(testURI, counter++, 1, 0.25, "08");
            launch_testStencilHeatSparseZ_perf<16, 256>(testURI, counter++, 1, 0.25, "08");
            launch_testStencilHeatSparseZ_perf<32, 128>(testURI, counter++, 1, 0.25, "08");

            testSet.insert(testURI);
        }

        BOOST_AUTO_TEST_CASE(testStencilHeatSparseZ09_gridScaling)
        {
            std::string testURI = suiteURI + ".device.stencil.sparse.Z.2D.09.gridScaling";
            unsigned int counter = 0;
            constexpr unsigned int blockEdgeSize = 8;
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 128>(testURI, counter++, 1, 0.1, "09");
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 256>(testURI, counter++, 1, 0.1, "09");
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 512>(testURI, counter++, 1, 0.1, "09");
            launch_testStencilHeatSparseZ_perf<blockEdgeSize, 1024>(testURI, counter++, 1, 0.1, "09");
////    launch_testStencilHeatSparseZ_perf<blockEdgeSize, 2048>(testURI, counter++, 1, 0.1, "09");

            testSet.insert(testURI);
        }
        BOOST_AUTO_TEST_CASE(testStencilHeatSparseZ09_blockScaling)
        {
            std::string testURI = suiteURI + ".device.stencil.sparse.Z.2D.09.blockScaling";
            unsigned int counter = 0;
            // Note - blockEdgeSize == 2 doesn't work
            launch_testStencilHeatSparseZ_perf<2, 2048>(testURI, counter++, 1, 0.1, "09");
            launch_testStencilHeatSparseZ_perf<4, 1024>(testURI, counter++, 1, 0.1, "09");
            launch_testStencilHeatSparseZ_perf<8, 512>(testURI, counter++, 1, 0.1, "09");
            launch_testStencilHeatSparseZ_perf<16, 256>(testURI, counter++, 1, 0.1, "09");
            launch_testStencilHeatSparseZ_perf<32, 128>(testURI, counter++, 1, 0.1, "09");

            testSet.insert(testURI);
        }

BOOST_AUTO_TEST_CASE(testStencilHeat3DSparse_gridScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.3D.gridScaling";
    unsigned int counter = 0;
    constexpr unsigned int blockEdgeSize = 8;
    launch_testStencilHeat3DSparse_perf<blockEdgeSize, 8>(testURI, counter++, 1, 1);
    launch_testStencilHeat3DSparse_perf<blockEdgeSize, 16>(testURI, counter++, 1, 1);
    launch_testStencilHeat3DSparse_perf<blockEdgeSize, 32>(testURI, counter++, 1, 1);
    launch_testStencilHeat3DSparse_perf<blockEdgeSize, 64>(testURI, counter++, 1, 1);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testStencilHeat3DSparse_blockScaling)
{
    std::string testURI = suiteURI + ".device.stencil.sparse.N.3D.blockScaling";
    unsigned int counter = 0;
    launch_testStencilHeat3DSparse_perf<2, 128>(testURI, counter++, 1, 1);
    launch_testStencilHeat3DSparse_perf<4, 64>(testURI, counter++, 1, 1);
    launch_testStencilHeat3DSparse_perf<8, 32>(testURI, counter++, 1, 1);
//    launch_testStencilHeat3DSparse_perf<16, 16>(testURI, counter++, 1, 1); // Too big, it doesn't work

    testSet.insert(testURI);
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void testInsertStencil(std::string testURI, unsigned int i)
{
	auto testName = "Insert stencil";
	constexpr unsigned int dim = 2;
//	constexpr unsigned int blockEdgeSize = 8;
	constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
	typedef aggregate<float,float> AggregateT;
	typedef HeatStencil<dim,0,1> StencilT;

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
		cudaDeviceSynchronize();

		timer ts;
		ts.start();

		sparseGrid.template applyStencils<StencilT>(STENCIL_MODE_INSERT, 0.1);
		sparseGrid.template flush<smax_<0>>(ctx, flush_type::FLUSH_ON_DEVICE);

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

BOOST_AUTO_TEST_CASE(testInsert_gridScaling)
{
    std::string testURI = suiteURI + ".device.insert.dense.single.2D.gridScaling";
    unsigned int counter = 0;
    testInsertSingle<8, 64>(testURI, counter++);
    testInsertSingle<8, 128>(testURI, counter++);
    testInsertSingle<8, 256>(testURI, counter++);
//    testInsertSingle<8, 512>(testURI, counter++);
//    testInsertSingle<8, 1024>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testInsert_blockScaling)
{
    std::string testURI = suiteURI + ".device.insert.dense.single.2D.blockScaling";
    unsigned int counter = 0;
    testInsertSingle<2, 1024>(testURI, counter++);
    testInsertSingle<4, 512>(testURI, counter++);
    testInsertSingle<8, 256>(testURI, counter++);
//    testInsertSingle<16, 128>(testURI, counter++);
//    testInsertSingle<32, 64>(testURI, counter++);

    testSet.insert(testURI);
}

template<unsigned int blockEdgeSize, unsigned int gridEdgeSize>
void test_insert_block(std::string testURI, unsigned int i)
{
	auto testName = "Insert (one chunk per block)";
	constexpr unsigned int dim = 2;
//	constexpr unsigned int blockEdgeSize = 8;
	constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
	typedef aggregate<float> AggregateT;

//	std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").insert");
    std::string base(testURI + "(" + std::to_string(i) + ")");
    report_sparsegrid_funcs.graphs.put(base + ".test.name","InsertBlock");

	report_sparsegrid_funcs.graphs.put(base + ".name","Block insert");
    report_sparsegrid_funcs.graphs.put(base + ".dim",dim);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*blockEdgeSize);

	unsigned int iterations = 10;

	openfpm::vector<double> measures;

	unsigned long long numElements = gridEdgeSize*blockEdgeSize*gridEdgeSize*blockEdgeSize;
	dim3 gridSize(gridEdgeSize, gridEdgeSize);
	dim3 blockSize(blockEdgeSize, blockEdgeSize);
	dim3 blockSizeBlockedInsert(1, 1);
	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// Warmup
	for (unsigned int iter=0; iter<5; ++iter)
	{
		auto offset = 0;
		sparseGrid.setGPUInsertBuffer(gridSize, blockSizeBlockedInsert);
		insertValues2DBlocked<0, 1, blockEdgeSize> << < gridSize, blockSize >> >
				(sparseGrid.toKernel(), offset, offset);
		sparseGrid.template flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
	}


	cudaDeviceSynchronize();


	for (unsigned int iter=0; iter<iterations; ++iter)
	{
		auto offset = 0;

		cudaDeviceSynchronize();

		timer ts;
		ts.start();

		sparseGrid.setGPUInsertBuffer(gridSize, blockSizeBlockedInsert);
		insertValues2DBlocked<0, 1, blockEdgeSize> << < gridSize, blockSize >> >
				(sparseGrid.toKernel(), offset, offset);
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
	std::cout << "\tInsert: " << mean << " dev: " << deviation << " s" << std::endl;
	std::cout << "Throughput:\n\t" << mean << " MElem/s\n";
}

BOOST_AUTO_TEST_CASE(testInsertBlocked_gridScaling)
{
    std::string testURI = suiteURI + ".device.insert.dense.block.2D.gridScaling";
    unsigned int counter = 0;
    test_insert_block<8,64>(testURI, counter++);
    test_insert_block<8,128>(testURI, counter++);
    test_insert_block<8,256>(testURI, counter++);
    test_insert_block<8,512>(testURI, counter++);
    test_insert_block<8,1024>(testURI, counter++);
//    test_insert_block<8,2048>(testURI, counter++); // Out of memory

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(testInsertBlocked_blockScaling)
{
    std::string testURI = suiteURI + ".device.insert.dense.block.2D.blockScaling";
    unsigned int counter = 0;
    test_insert_block<2,2048>(testURI, counter++);
    test_insert_block<4,1024>(testURI, counter++);
    test_insert_block<8,512>(testURI, counter++);
    test_insert_block<16,256>(testURI, counter++);
    test_insert_block<32,128>(testURI, counter++);

    testSet.insert(testURI);
}

BOOST_AUTO_TEST_CASE(write_teport)
{
    write_test_report(report_sparsegrid_funcs, testSet);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
