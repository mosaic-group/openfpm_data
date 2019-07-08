//
// Created by tommaso on 4/07/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/SparseGridGpu.hpp"

// Error handling from https://stackoverflow.com/a/13245319
static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

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
struct HeatStencilOrig
{
    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int flops = 3 + 2*dim;

    static constexpr unsigned int supportRadius = 1;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            grid_key_dx<dim> & dataBlockCoord,
            unsigned int offset,
            grid_key_dx<dim> & pointCoord,
            DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore,
            float dt)
    {
        typedef typename SparseGridT::AggregateBlockType AggregateT;
        typedef ScalarTypeOf<AggregateT, p> ScalarT;

        constexpr unsigned int enlargedBlockSize = IntPow<
                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

        __shared__ ScalarT enlargedBlock[enlargedBlockSize];
        sparseGrid.loadBlock<p>(dataBlockLoad, enlargedBlock);
        sparseGrid.loadGhost<p>(dataBlockCoord, enlargedBlock);
        __syncthreads();

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

        __syncthreads();
        sparseGrid.storeBlock<p>(dataBlockStore, enlargedBlock);
    }

    template <typename SparseGridT, typename CtxT>
    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
    {
        sparseGrid.template flush <smin_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

template<unsigned int dim, unsigned int p>
struct HeatStencil
{
    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int flops = 3 + 2*dim;

    static constexpr unsigned int supportRadius = 1;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            grid_key_dx<dim> & dataBlockCoord,
            unsigned int offset,
            grid_key_dx<dim> & pointCoord,
            DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore,
            float dt, unsigned int maxIter=1000)
    {
        typedef typename SparseGridT::AggregateBlockType AggregateT;
        typedef ScalarTypeOf<AggregateT, p> ScalarT;
        typedef BlockTypeOf<AggregateT, p> BlockT;
        constexpr unsigned int blockSize = BlockT::size;

        constexpr unsigned int enlargedBlockSize = IntPow<
                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

        const auto coord = sparseGrid.getCoordInEnlargedBlock(offset);
        const auto linId = sparseGrid.getLinIdInEnlargedBlock(offset);
        char boundaryDirection[dim];
        bool isBoundary = sparseGrid.getIfBoundaryElementInEnlargedBlock(coord, boundaryDirection);

        unsigned int nPlusId[dim], nMinusId[dim];
        for (int d=0; d<dim; ++d)
        {
            nPlusId[d] = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, 1);
            nMinusId[d] = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, -1);
        }

        __shared__ ScalarT enlargedBlock[enlargedBlockSize];

        ScalarT * nPlus[dim];
        ScalarT * nMinus[dim];
        for (int d=0; d<dim; ++d)
        {
            const auto boundaryDir = boundaryDirection[d];
            const auto nCoord = sparseGrid.getNeighbour(pointCoord, d, boundaryDir);
            const auto nOffset = sparseGrid.getLinId(nCoord) % blockSize;
            nPlus[d] = &(enlargedBlock[nPlusId[d]]);
            nMinus[d] = &(enlargedBlock[nMinusId[d]]);
            if (boundaryDir==1)
            {
//                nPlus[d] = &(sparseGrid.getBlock(nCoord).template get<p>()[nOffset]);
//                nPlus[d] = &(sparseGrid.template get<p>(nCoord));
                nPlus[d] = sparseGrid.getBlock(nCoord).template get<p>().block + nOffset;
            }
            else if (boundaryDir==-1)
            {
//                nMinus[d] = &(sparseGrid.getBlock(nCoord).template get<p>()[nOffset]);
//                nMinus[d] = &(sparseGrid.template get<p>(nCoord));
                nMinus[d] = sparseGrid.getBlock(nCoord).template get<p>().block + nOffset;
            }
        }

        sparseGrid.loadBlock<p>(dataBlockLoad, enlargedBlock);
        __syncthreads();
//            sparseGrid.loadGhost<p>(dataBlockCoord, enlargedBlock);
//            __syncthreads();

        ScalarT cur = enlargedBlock[linId];
        ScalarT laplacian = -2.0 * dim * cur; // The central part of the stencil
        for (int d = 0; d < dim; ++d)
        {
//                const auto boundary = boundaryDirection[d];
//                ScalarT neighbourPlus = enlargedBlock[nPlusId[d]];
//                ScalarT neighbourMinus = enlargedBlock[nMinusId[d]];
//                if (boundary == 1)
//                {
//                    neighbourPlus = *(nPlus[d]);
//                }
//                else if (boundary == -1)
//                {
//                    neighbourMinus = *(nMinus[d]);
//                }
            laplacian += *(nMinus[d]) + *(nPlus[d]);
//                laplacian += neighbourMinus + neighbourPlus;

        }
        enlargedBlock[linId] = cur + dt * laplacian;

        __syncthreads();
        sparseGrid.storeBlock<p>(dataBlockStore, enlargedBlock);
        __syncthreads();
    }

    template <typename SparseGridT, typename CtxT>
    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
    {
        sparseGrid.template flush <smin_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

BOOST_AUTO_TEST_SUITE(SparseGridGpu_Stencil_Performance_tests)
    BOOST_AUTO_TEST_CASE(testStencilHeatPerf)
    {
        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
        unsigned int gridEdgeSize = 256;
        constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
        typedef aggregate<float> AggregateT;
//        typedef HeatStencil<dim, 0> StencilT;
        typedef HeatStencilOrig<dim, 0> StencilT;

        unsigned int iterations = 5;
        unsigned int repetitions = 1;

        float timeInitAvg;
        float timeStencilAvg;
        float timeTotalAvg;

        for (int rep=0; rep<repetitions; ++rep)
        {

            cudaEvent_t start, afterInit, stop;
            float timeInit;
            float timeStencil;
            float timeTotal;

            HANDLE_ERROR(cudaEventCreate(&start));
            HANDLE_ERROR(cudaEventCreate(&afterInit));
            HANDLE_ERROR(cudaEventCreate(&stop));

            HANDLE_ERROR(cudaEventRecord(start, 0));

            dim3 gridSize(gridEdgeSize, gridEdgeSize);
            dim3 blockSize(blockEdgeSize, blockEdgeSize);
            BlockGeometry<dim, blockEdgeSize> blockGeometry(gridSize);
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

            HANDLE_ERROR(cudaEventRecord(afterInit, 0));
            HANDLE_ERROR(cudaEventSynchronize(afterInit));

            for (unsigned int iter=0; iter<iterations; ++iter)
            {
                sparseGrid.applyStencils<StencilT>(STENCIL_MODE_INPLACE, 0.1);
                cudaDeviceSynchronize();
            }

            HANDLE_ERROR(cudaEventRecord(stop, 0));
            HANDLE_ERROR(cudaEventSynchronize(stop));
            HANDLE_ERROR(cudaEventElapsedTime(&timeInit, start, afterInit));
            HANDLE_ERROR(cudaEventElapsedTime(&timeStencil, afterInit, stop));
            HANDLE_ERROR(cudaEventElapsedTime(&timeTotal, start, stop));

            timeInitAvg += timeInit;
            timeStencilAvg += timeStencil;
            timeTotalAvg += timeTotal;
        }

        timeInitAvg /= repetitions;
        timeStencilAvg /= repetitions;
        timeTotalAvg /= repetitions;

        unsigned int numElements = gridEdgeSize*blockEdgeSize*gridEdgeSize*blockEdgeSize;
        float gElemS = numElements * iterations / (1e9 * timeStencilAvg/1000);
        float gFlopsS = gElemS * StencilT::flops;
        printf("Grid: %ux%u\n", gridEdgeSize*blockEdgeSize, gridEdgeSize*blockEdgeSize);
        printf("Iterations: %u\n", iterations);
        printf("Timing (avg on %u repetitions):\n\tInit: %f ms\n\tStencil: %f ms\n\tTotal: %f ms\n",
                repetitions, timeInitAvg, timeStencilAvg, timeTotalAvg);
        printf("Throughput:\n\t%f GElem/s\n\t%f GFlops/s\n", gElemS, gFlopsS);

    }
BOOST_AUTO_TEST_SUITE_END()
