//
// Created by tommaso on 4/07/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/SparseGridGpu.hpp"
#include "cuda_macro.h"


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
    const unsigned int dataBlockIdX = posX / blockEdgeSize;
    const unsigned int dataBlockIdY = posY / blockEdgeSize;
    const unsigned int offsetX = posX % blockEdgeSize;
    const unsigned int offsetY = posY % blockEdgeSize;

    const unsigned int blockDimX = blockDim.x / blockEdgeSize;
    const unsigned int blockOffsetX = threadIdx.x / blockEdgeSize;
    const unsigned int blockOffsetY = threadIdx.y / blockEdgeSize;

    const unsigned int dataBlockNum = blockOffsetY*blockDimX + blockOffsetX;
    const unsigned int offset = offsetY * blockEdgeSize + offsetX;

    if (offset == 0) // Just one thread per data block
    {
        grid_key_dx<SparseGridType::d, int> blockCoord({dataBlockIdX, dataBlockIdY});
        auto encap = sparseGrid.insertBlock(sparseGrid.getBlockLinId(blockCoord));
        blocks[dataBlockNum] = &(encap.template get<p>());
        masks[dataBlockNum] = &(encap.template get<pMask>());
    }

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

    if (offset == 0) // Just one thread per data block
    {
        auto encap = sparseGrid.insertBlock(dataBlockId);
        blocks[dataBlockNum] = &(encap.template get<p>());
        masks[dataBlockNum] = &(encap.template get<pMask>());
    }

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

//    printf("copyBlocksToOutput: bDim=(%d,%d), bId=(%d,%d), tId=(%d,%d) : "
//           "pos=%ld, coord={%d,%d}, value=%d\n",
//           bDimX, bDimY,
//           bIdX, bIdY,
//           tIdX, tIdY,
//           pos,
//           x, y,
//           static_cast<int>(value)); //debug

    output.template get<p>(pos) = value;

    // Compiler avoid warning
    x++;
    y++;
    z++;
}

template<unsigned int dim, unsigned int p>
struct EmptyStencil
{
    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int flops = 0;

    static constexpr unsigned int supportRadius = 0;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            const unsigned int dataBlockId,
            const int * neighboursPositions,
            unsigned int offset,
            grid_key_dx<dim, int> & pointCoord,
            DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore,
            bool applyStencilHere,
            float dt)
    {
        // Nothing to do here
    }

    template <typename SparseGridT, typename CtxT>
    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
    {
        sparseGrid.template flush <smin_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

template<unsigned int dim, unsigned int p>
struct SkeletonStencil
{
    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int flops = 0;

    static constexpr unsigned int supportRadius = 1;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            const unsigned int dataBlockId,
            const int * neighboursPositions,
            unsigned int offset,
            grid_key_dx<dim, int> & pointCoord,
            DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore,
            bool applyStencilHere,
            float dt)
    {
        typedef typename SparseGridT::AggregateBlockType AggregateT;
        typedef ScalarTypeOf<AggregateT, p> ScalarT;

        constexpr unsigned int enlargedBlockSize = IntPow<
                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

        __shared__ ScalarT enlargedBlock[enlargedBlockSize];
        sparseGrid.loadBlock<p>(dataBlockLoad, enlargedBlock);
        sparseGrid.loadGhost<p>(dataBlockId, neighboursPositions, enlargedBlock);
        __syncthreads();
        sparseGrid.storeBlock<p>(dataBlockStore, enlargedBlock);
    }

    template <typename SparseGridT, typename CtxT>
    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
    {
        sparseGrid.template flush <sRight_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

template<unsigned int dim, unsigned int p_src, unsigned int p_dst>
struct HeatStencil
{
    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int flops = 3 + 2*dim;

    static constexpr unsigned int supportRadius = 1;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            const unsigned int dataBlockId,
            const int * neighboursPositions,
            unsigned int offset,
            grid_key_dx<dim, int> & pointCoord,
            const DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore,
            bool applyStencilHere,
            float dt)
    {
        typedef typename SparseGridT::AggregateBlockType AggregateT;
        typedef ScalarTypeOf<AggregateT, p_src> ScalarT;

        constexpr unsigned int enlargedBlockSize = IntPow<
                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

        __shared__ ScalarT enlargedBlock[enlargedBlockSize];
//        sparseGrid.loadGhost<p_src>(dataBlockId, neighboursPositions, enlargedBlock);
//        sparseGrid.loadBlock<p_src>(dataBlockLoad, enlargedBlock);

        sparseGrid.loadGhostBlock<p_src>(dataBlockLoad,dataBlockId, neighboursPositions, enlargedBlock);

        __syncthreads();

        if (applyStencilHere)
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

        __syncthreads();
        sparseGrid.storeBlock<p_dst>(dataBlockStore, enlargedBlock);
    }

    template <typename SparseGridT, typename CtxT>
    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
    {
        sparseGrid.template flush <sRight_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

template<typename SparseGridZ>
void testStencilHeat_perf()
{
    auto testName = "In-place stencil";
    unsigned int gridEdgeSize = 1024;
//        constexpr unsigned int blockEdgeSize = 16;
//        unsigned int gridEdgeSize = 256;
//        unsigned int gridEdgeSize = 8;
//        typedef EmptyStencil<dim, 0> StencilT;
//        typedef SkeletonStencil<dim, 0> StencilT;
    typedef HeatStencil<SparseGridZ::dims,0,1> StencilT;

    unsigned int iterations = 100;
    unsigned int repetitions = 5;

    float timeStencilAvg = 0.0;

    for (int rep=0; rep<repetitions; ++rep)
    {
        dim3 gridSize(gridEdgeSize, gridEdgeSize);
        dim3 blockSize(SparseGridZ::blockEdgeSize_,SparseGridZ::blockEdgeSize_);
        typename SparseGridZ::grid_info blockGeometry(gridSize);
        SparseGridZ sparseGrid(blockGeometry);
        mgpu::ofp_context_t ctx;
        sparseGrid.template setBackgroundValue<0>(0);

        // Initialize the grid
        sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
        insertConstantValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), 0);
        sparseGrid.template flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
        dim3 sourcePt(gridSize.x * SparseGridZ::blockEdgeSize_ / 2, gridSize.y * SparseGridZ::blockEdgeSize_ / 2, 0);
        insertOneValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), sourcePt, 100);
        sparseGrid.template flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

        cudaDeviceSynchronize();

        timer ts;
        ts.start();

        for (unsigned int iter=0; iter<iterations; ++iter)
        {sparseGrid.template applyStencils<StencilT>(STENCIL_MODE_INPLACE, 0.1);}

        cudaDeviceSynchronize();

        ts.stop();

        timeStencilAvg += ts.getwct();
    }

    timeStencilAvg /= repetitions;

    // All times above are in ms

    unsigned long long numElements = gridEdgeSize*SparseGridZ::blockEdgeSize_*gridEdgeSize*SparseGridZ::blockEdgeSize_;
    float gElemS = numElements * iterations / (1e9 * timeStencilAvg);
    float gFlopsS = gElemS * StencilT::flops;
    float stencilSingleTimingMillis = timeStencilAvg/iterations;
    std::cout << "Test: " << testName << std::endl;
    std::cout << "Grid: " << gridEdgeSize*SparseGridZ::blockEdgeSize_ << "x" << gridEdgeSize*SparseGridZ::blockEdgeSize_ << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "Timing (avg on " << repetitions << " repetitions):" << std::endl
    		  << "\tStencil: " << timeStencilAvg << " ms" << std::endl;
    std::cout << "Stencil details:" << std::endl << "\tSingle application timing: " << stencilSingleTimingMillis << " ms" << std::endl;
    std::cout << "Throughput: " << std::endl << "\t " << gElemS << " GElem/s " << std::endl << "\t " << gFlopsS << " GFlops/s" << std::endl;
}

BOOST_AUTO_TEST_SUITE(SparseGridGpu_Stencil_Performance_tests)

    BOOST_AUTO_TEST_CASE(testStencilHeat)
    {
    	constexpr unsigned int dim = 2;
    	constexpr unsigned int blockEdgeSize = 8;

    	typedef aggregate<float,float> AggregateT;
        constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

		testStencilHeat_perf<SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize>>();
    }

BOOST_AUTO_TEST_CASE(testStencilHeatZ)
{
	constexpr unsigned int dim = 2;
	constexpr unsigned int blockEdgeSize = 8;

	typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

	testStencilHeat_perf<SparseGridGpu_z<dim, AggregateT, blockEdgeSize, chunkSize>>();
}

    BOOST_AUTO_TEST_CASE(testStencilHeatInsert)
    {
        auto testName = "Insert stencil";
        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
        unsigned int gridEdgeSize = 512;
//        constexpr unsigned int blockEdgeSize = 16;
//        unsigned int gridEdgeSize = 256;
//        unsigned int gridEdgeSize = 8;
        constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
        typedef aggregate<float,float> AggregateT;
//        typedef EmptyStencil<dim, 0> StencilT;
//        typedef SkeletonStencil<dim, 0> StencilT; //todo: Rimetti HeatStencil!
        typedef HeatStencil<dim,0,1> StencilT;

        unsigned int iterations = 10;
        unsigned int repetitions = 5;
//        unsigned int iterations = 5;
//        unsigned int repetitions = 1;

        float timeInitAvg;
        float timeStencilAvg;
        float timeTotalAvg;

        for (int rep=0; rep<repetitions; ++rep)
        {

            cudaEvent_t start, afterInit, stop;
            float timeInit;
            float timeStencil;
            float timeTotal;

            CUDA_SAFE_CALL(cudaEventCreate(&start));
            CUDA_SAFE_CALL(cudaEventCreate(&afterInit));
            CUDA_SAFE_CALL(cudaEventCreate(&stop));

            CUDA_SAFE_CALL(cudaEventRecord(start, 0));

            dim3 gridSize(gridEdgeSize, gridEdgeSize);
            dim3 blockSize(blockEdgeSize, blockEdgeSize);
            grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize> sparseGrid(blockGeometry);
            mgpu::ofp_context_t ctx;
            sparseGrid.template setBackgroundValue<0>(0);

            // Initialize the grid
            sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
            insertConstantValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), 0);
            sparseGrid.flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

            sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
            dim3 sourcePt(gridSize.x * blockEdgeSize / 2, gridSize.y * blockEdgeSize / 2, 0);
            insertOneValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), sourcePt, 100);
            sparseGrid.flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

            sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

            CUDA_SAFE_CALL(cudaEventRecord(afterInit, 0));
            CUDA_SAFE_CALL(cudaEventSynchronize(afterInit));

            for (unsigned int iter=0; iter<iterations; ++iter)
            {
                sparseGrid.applyStencils<StencilT>(STENCIL_MODE_INSERT, 0.1);
                cudaDeviceSynchronize();
            }

            CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
            CUDA_SAFE_CALL(cudaEventSynchronize(stop));
            CUDA_SAFE_CALL(cudaEventElapsedTime(&timeInit, start, afterInit));
            CUDA_SAFE_CALL(cudaEventElapsedTime(&timeStencil, afterInit, stop));
            CUDA_SAFE_CALL(cudaEventElapsedTime(&timeTotal, start, stop));

            timeInitAvg += timeInit;
            timeStencilAvg += timeStencil;
            timeTotalAvg += timeTotal;
        }

        timeInitAvg /= repetitions;
        timeStencilAvg /= repetitions;
        timeTotalAvg /= repetitions;

        // All times above are in ms

        unsigned long long numElements = gridEdgeSize*blockEdgeSize*gridEdgeSize*blockEdgeSize;
        float gElemS = numElements * iterations / (1e9 * timeStencilAvg/1000);
        float gFlopsS = gElemS * StencilT::flops;
        float stencilSingleTimingMillis = timeStencilAvg/iterations;
        printf("Test: %s\n", testName);
        printf("Grid: %ux%u\n", gridEdgeSize*blockEdgeSize, gridEdgeSize*blockEdgeSize);
        printf("Iterations: %u\n", iterations);
        printf("Timing (avg on %u repetitions):\n\tInit: %f ms\n\tStencil: %f ms\n\tTotal: %f ms\n",
               repetitions, timeInitAvg, timeStencilAvg, timeTotalAvg);
        printf("Stencil details:\n\tSingle application timing: %f ms\n", stencilSingleTimingMillis);
        printf("Throughput:\n\t%f GElem/s\n\t%f GFlops/s\n", gElemS, gFlopsS);

    }

    BOOST_AUTO_TEST_CASE(testInsert)
    {
        auto testName = "Insert (one chunk per element)";
        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
//        unsigned int gridEdgeSize = 512;
        unsigned int gridEdgeSize = 128;
//        constexpr unsigned int blockEdgeSize = 16;
//        unsigned int gridEdgeSize = 256;
//        unsigned int gridEdgeSize = 8;
        constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
        typedef aggregate<float> AggregateT;

        unsigned int iterations = 10;
        unsigned int repetitions = 5;
        bool prePopulateGrid = true;

        float timeInitAvg;
        float timeStencilAvg;
        float timeTotalAvg;

        for (int rep=0; rep<repetitions; ++rep)
        {

            cudaEvent_t start, afterInit, stop;
            float timeInit;
            float timeStencil;
            float timeTotal;

            CUDA_SAFE_CALL(cudaEventCreate(&start));
            CUDA_SAFE_CALL(cudaEventCreate(&afterInit));
            CUDA_SAFE_CALL(cudaEventCreate(&stop));

            CUDA_SAFE_CALL(cudaEventRecord(start, 0));

            dim3 gridSize(gridEdgeSize, gridEdgeSize);
            dim3 blockSize(blockEdgeSize, blockEdgeSize);
            grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize> sparseGrid(blockGeometry);
            mgpu::ofp_context_t ctx;
            sparseGrid.template setBackgroundValue<0>(0);

            // Initialize the grid
//            sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
//            insertConstantValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), 0);
//            sparseGrid.flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//
//            sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
//            dim3 sourcePt(gridSize.x * blockEdgeSize / 2, gridSize.y * blockEdgeSize / 2, 0);
//            insertOneValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), sourcePt, 100);
//            sparseGrid.flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//
//            sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

            if (prePopulateGrid)
            {
                // Pre-populate grid
                sparseGrid.setGPUInsertBuffer(gridSize, blockSize);
                insertValues2D<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), 0, 0);
                sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
                cudaDeviceSynchronize();
                ///
            }

            CUDA_SAFE_CALL(cudaEventRecord(afterInit, 0));
            CUDA_SAFE_CALL(cudaEventSynchronize(afterInit));

            for (unsigned int iter=0; iter<iterations; ++iter)
            {
//                auto offset = iter * 99999 % 32003;
                auto offset = 0;
                sparseGrid.setGPUInsertBuffer(gridSize, blockSize);
                insertValues2D<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), offset, offset);
                sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
                cudaDeviceSynchronize();
            }

            CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
            CUDA_SAFE_CALL(cudaEventSynchronize(stop));
            CUDA_SAFE_CALL(cudaEventElapsedTime(&timeInit, start, afterInit));
            CUDA_SAFE_CALL(cudaEventElapsedTime(&timeStencil, afterInit, stop));
            CUDA_SAFE_CALL(cudaEventElapsedTime(&timeTotal, start, stop));

            timeInitAvg += timeInit;
            timeStencilAvg += timeStencil;
            timeTotalAvg += timeTotal;
        }

        timeInitAvg /= repetitions;
        timeStencilAvg /= repetitions;
        timeTotalAvg /= repetitions;

        // All times above are in ms

        unsigned long long numElements = gridEdgeSize*blockEdgeSize*gridEdgeSize*blockEdgeSize;
        float mElemS = numElements * iterations / (1e6 * timeStencilAvg/1000);
//        float gFlopsS = gElemS * StencilT::flops;
        float stencilSingleTimingMillis = timeStencilAvg/iterations;
        printf("Test: %s\n", testName);
        printf("Grid: %ux%u\n", gridEdgeSize*blockEdgeSize, gridEdgeSize*blockEdgeSize);
        printf("Iterations: %u\n", iterations);
        printf("Timing (avg on %u repetitions):\n\tInit: %f ms\n\tStencil: %f ms\n\tTotal: %f ms\n",
               repetitions, timeInitAvg, timeStencilAvg, timeTotalAvg);
        printf("Stencil details:\n\tSingle application timing: %f ms\n", stencilSingleTimingMillis);
//        printf("Throughput:\n\t%f GElem/s\n\t%f GFlops/s\n", gElemS, gFlopsS);
        printf("Throughput:\n\t%f MElem/s\n", mElemS);

    }

    BOOST_AUTO_TEST_CASE(testInsertBlocked)
    {
        auto testName = "Insert (one chunk per block)";
        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
//        unsigned int gridEdgeSize = 512;
        unsigned int gridEdgeSize = 128;
//        constexpr unsigned int blockEdgeSize = 16;
//        unsigned int gridEdgeSize = 256;
//        unsigned int gridEdgeSize = 8;
        constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
        typedef aggregate<float> AggregateT;

        unsigned int iterations = 10;
        unsigned int repetitions = 5;
        bool prePopulateGrid = true;

        float timeInitAvg;
        float timeStencilAvg;
        float timeTotalAvg;

        for (int rep=0; rep<repetitions; ++rep)
        {

            cudaEvent_t start, afterInit, stop;
            float timeInit;
            float timeStencil;
            float timeTotal;

            CUDA_SAFE_CALL(cudaEventCreate(&start));
            CUDA_SAFE_CALL(cudaEventCreate(&afterInit));
            CUDA_SAFE_CALL(cudaEventCreate(&stop));

            CUDA_SAFE_CALL(cudaEventRecord(start, 0));

            dim3 gridSize(gridEdgeSize, gridEdgeSize);
            dim3 blockSize(blockEdgeSize, blockEdgeSize);
            dim3 blockSizeBlockedInsert(1, 1);
            grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
            SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize> sparseGrid(blockGeometry);
            mgpu::ofp_context_t ctx;
            sparseGrid.template setBackgroundValue<0>(0);

            // Initialize the grid
//            sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
//            insertConstantValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), 0);
//            sparseGrid.flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//
//            sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
//            dim3 sourcePt(gridSize.x * blockEdgeSize / 2, gridSize.y * blockEdgeSize / 2, 0);
//            insertOneValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), sourcePt, 100);
//            sparseGrid.flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
//
//            sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

            if (prePopulateGrid)
            {
                // Pre-populate grid
                sparseGrid.setGPUInsertBuffer(gridSize, blockSize);
                insertValues2D<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), 0, 0);
                sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
                cudaDeviceSynchronize();
                ///
            }

            CUDA_SAFE_CALL(cudaEventRecord(afterInit, 0));
            CUDA_SAFE_CALL(cudaEventSynchronize(afterInit));

            for (unsigned int iter=0; iter<iterations; ++iter)
            {
//                auto offset = iter * 99999 % 32003;
                auto offset = 0;
                sparseGrid.setGPUInsertBuffer(gridSize, blockSizeBlockedInsert);
                insertValues2DBlocked<0, 1, blockEdgeSize> << < gridSize, blockSize >> >
                        (sparseGrid.toKernel(), offset, offset);
                sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
                cudaDeviceSynchronize();
            }

            CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
            CUDA_SAFE_CALL(cudaEventSynchronize(stop));
            CUDA_SAFE_CALL(cudaEventElapsedTime(&timeInit, start, afterInit));
            CUDA_SAFE_CALL(cudaEventElapsedTime(&timeStencil, afterInit, stop));
            CUDA_SAFE_CALL(cudaEventElapsedTime(&timeTotal, start, stop));

            timeInitAvg += timeInit;
            timeStencilAvg += timeStencil;
            timeTotalAvg += timeTotal;
        }

        timeInitAvg /= repetitions;
        timeStencilAvg /= repetitions;
        timeTotalAvg /= repetitions;

        // All times above are in ms

        unsigned long long numElements = gridEdgeSize*blockEdgeSize*gridEdgeSize*blockEdgeSize;
        float mElemS = numElements * iterations / (1e6 * timeStencilAvg/1000);
//        float gFlopsS = gElemS * StencilT::flops;
        float stencilSingleTimingMillis = timeStencilAvg/iterations;
        printf("Test: %s\n", testName);
        printf("Grid: %ux%u\n", gridEdgeSize*blockEdgeSize, gridEdgeSize*blockEdgeSize);
        printf("Iterations: %u\n", iterations);
        printf("Timing (avg on %u repetitions):\n\tInit: %f ms\n\tStencil: %f ms\n\tTotal: %f ms\n",
               repetitions, timeInitAvg, timeStencilAvg, timeTotalAvg);
        printf("Stencil details:\n\tSingle application timing: %f ms\n", stencilSingleTimingMillis);
//        printf("Throughput:\n\t%f GElem/s\n\t%f GFlops/s\n", gElemS, gFlopsS);
        printf("Throughput:\n\t%f MElem/s\n", mElemS);

    }

BOOST_AUTO_TEST_SUITE_END()
