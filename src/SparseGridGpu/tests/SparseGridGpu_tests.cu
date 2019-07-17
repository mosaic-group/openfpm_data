//
// Created by tommaso on 10/06/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/SparseGridGpu.hpp"

template<unsigned int p, typename SparseGridType>
__global__ void insertValues(SparseGridType sparseGrid)
{
    sparseGrid.init();

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
//    printf("insertValues: bDim=(%d,%d), bId=(%d,%d), tId=(%d,%d) : "
//           "pos=%ld, coord={%d,%d}, value=%d\n",
//           bDimX, bDimY,
//           bIdX, bIdY,
//           tIdX, tIdY,
//           pos,
//           x, y,
//           x); //debug

    sparseGrid.template insert<p>(coord) = x;

    __syncthreads();

    sparseGrid.flush_block_insert();

    // Compiler avoid warning
    x++;
    y++;
    z++;
}

template<unsigned int p, typename SparseGridType, typename ScalarT>
__global__ void insertConstantValue(SparseGridType sparseGrid, ScalarT value)
{
    sparseGrid.init();

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

    sparseGrid.template insert<p>(coord) = value;

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

template<unsigned int p, typename SparseGridType, typename VectorOutType>
__global__ void copyToOutputIfPadding(SparseGridType sparseGrid, VectorOutType output)
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

//    auto value = sparseGrid.template get<p>(coord);
//    auto mask = sparseGrid.template get<1>(coord);
//    if (value == 1)
//    {
//        printf("copyBlocksToOutput: bDim=(%d,%d), bId=(%d,%d), tId=(%d,%d) : "
//               "pos=%ld, coord={%d,%d}, value=%d, mask=%u\n",
//               bDimX, bDimY,
//               bIdX, bIdY,
//               tIdX, tIdY,
//               pos,
//               x, y,
//               static_cast<int>(value),
//               static_cast<unsigned char>(mask)); //debug
//    }
//
//    if (sparseGrid.isPadding(coord))
//    {
//        printf("OUTPUT : Element isPadding! pos=%u\n",
//               pos); //debug
//    }

    output.template get<p>(pos) = sparseGrid.isPadding(coord) ? 1 : 0;

    // Compiler avoid warning
    x++;
    y++;
    z++;
}

template<unsigned int p, typename SparseGridType>
__global__ void insertBoundaryValuesHeat(SparseGridType sparseGrid)
{
    sparseGrid.init();

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

    float value = 0;
    if (x == 0)
    {
        value = 0;
    }
    else if (x == bDimX * gridDim.x - 1)
    {
        value = 10;
    }

    if (y == 0 || y == bDimY * gridDim.y - 1)
    {
        value = 10.0 * x / (bDimX * gridDim.x - 1);
    }

    sparseGrid.template insert<p>(coord) = value;

    __syncthreads();

    sparseGrid.flush_block_insert();

    // Compiler avoid warning
    x++;
    y++;
    z++;
}

template<unsigned int dim, unsigned int p>
struct LaplacianStencil
{
    // This is an example of a laplacian stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int supportRadius = 1;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            grid_key_dx<dim, int> & dataBlockCoord,
            unsigned int offset,
            grid_key_dx<dim, int> & pointCoord,
            DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore)
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
        ScalarT res = -2.0*dim*cur; // The central part of the stencil
        for (int d=0; d<dim; ++d)
        {
            auto nPlusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, 1);
            auto nMinusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, -1);
            ScalarT neighbourPlus = enlargedBlock[nPlusId];
            ScalarT neighbourMinus = enlargedBlock[nMinusId];
            res += neighbourMinus + neighbourPlus;
        }
        enlargedBlock[linId] = res;

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
//        sparseGrid.loadGhost<p>(dataBlockId, nullptr, enlargedBlock);
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
        sparseGrid.storeBlock<p>(dataBlockStore, enlargedBlock);
    }

    template <typename SparseGridT, typename CtxT>
    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
    {
        sparseGrid.template flush <smax_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

template<unsigned int dim, unsigned int p>
struct HeatStencil2
{
    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int flops = 3 + 2*dim;

    static constexpr unsigned int supportRadius = 1;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            grid_key_dx<dim, int> & dataBlockCoord,
            unsigned int offset,
            grid_key_dx<dim, int> & pointCoord,
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
//                nPlus[d] = sparseGrid.getBlock(nCoord).template getPtr<p>()->block + nOffset;
                nPlus[d] = sparseGrid.getBlock(nCoord).template get<p>().block + nOffset;
            }
            else if (boundaryDir==-1)
            {
//                nMinus[d] = &(sparseGrid.getBlock(nCoord).template get<p>()[nOffset]);
//                nMinus[d] = &(sparseGrid.template get<p>(nCoord));
//                nMinus[d] = sparseGrid.getBlock(nCoord).template getPtr<p>()->block + nOffset;
                nMinus[d] = sparseGrid.getBlock(nCoord).template get<p>().block + nOffset;
            }
        }

        sparseGrid.loadBlock<p>(dataBlockLoad, enlargedBlock);
        __syncthreads();
        for (unsigned int iter=0; iter<maxIter; ++iter)
        {
//            sparseGrid.loadGhost<p>(dataBlockCoord, enlargedBlock);
//            __syncthreads();

//todo: capisci come mai questa load non va mentre con la load ghost si!

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
//                laplacian += neighbourMinus + neighbourPlus;
                laplacian += *(nMinus[d]) + *(nPlus[d]);
            }
            enlargedBlock[linId] = cur + dt * laplacian;

            __syncthreads();
//            sparseGrid.storeBlock<p>(dataBlockLoad, enlargedBlock);
            if (isBoundary)
            {
                dataBlockLoad.template get<p>()[offset] = enlargedBlock[linId];
            }
            __syncthreads();
        }
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

//todo: Here write some sort of stencil kernel to test the client interface for stencil application

BOOST_AUTO_TEST_SUITE(SparseGridGpu_tests)

    BOOST_AUTO_TEST_CASE(testInsert)
    {

        std::cout << std::endl; //debug empty line

        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
        constexpr unsigned int dataBlockSize = blockEdgeSize * blockEdgeSize;
        typedef aggregate<DataBlock<float, dataBlockSize>> AggregateSGT;
        typedef aggregate<float> AggregateOutT;

        dim3 gridSize(2, 2);
        dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize);

        BlockGeometry<dim, blockEdgeSize> blockGeometry(gridSize);
        SparseGridGpu<dim, AggregateOutT, blockEdgeSize> sparseGrid(blockGeometry);

        sparseGrid.template setBackgroundValue<0>(666);
//        const unsigned int gridSizeLin = 4;
//        const unsigned int bufferPoolSize = 1024;
//        const unsigned int blockSizeInsert = 128;
//        const unsigned int gridSizeRead = gridSize + 1;
//        const unsigned int blockSizeRead = 128;

//        sparseGrid.setGPUInsertBuffer(gridSizeLin, bufferPoolSize);
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);

        insertValues<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel());

        mgpu::ofp_context_t ctx;
        sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        // Get output
        openfpm::vector_gpu<AggregateOutT> output;
        output.resize(4 * 64);

        copyBlocksToOutput<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
//            auto expectedValue = (i < gridSize * blockSizeInsert) ? i : 666;
            auto coord = sparseGrid.getCoord(i);
            auto expectedValue = coord.get(0);

            std::cout << "sparseGrid(" << coord.get(0) << "," << coord.get(1) << ") = "
                      << sparseGrid.template get<0>(coord)
                      << " == "
                      << expectedValue
                      << " == "
                      << output.template get<0>(i) << " = output(" << i << ")"
                      << std::endl;
            match &= output.template get<0>(i) == sparseGrid.template get<0>(coord);
            match &= output.template get<0>(i) == expectedValue;
        }

        BOOST_REQUIRE_EQUAL(match, true);
    }

    BOOST_AUTO_TEST_CASE(testInsert3D)
    {
        constexpr unsigned int dim = 3;
        constexpr unsigned int blockEdgeSize = 4;
        constexpr unsigned int dataBlockSize = blockEdgeSize * blockEdgeSize;
        typedef aggregate<DataBlock<float, dataBlockSize>> AggregateSGT;
        typedef aggregate<float> AggregateOutT;

        dim3 gridSize(2, 2, 2);
        dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize, blockEdgeSize);

        BlockGeometry<dim, blockEdgeSize> blockGeometry(gridSize);
        SparseGridGpu<dim, AggregateOutT, blockEdgeSize> sparseGrid(blockGeometry);

        sparseGrid.template setBackgroundValue<0>(666);

        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);

        insertValues<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel());

        mgpu::ofp_context_t ctx;
        sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        // Get output
        openfpm::vector_gpu<AggregateOutT> output;
        output.resize(sparseGrid.dim3SizeToInt(gridSize) * 64);

        copyBlocksToOutput<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
//            auto expectedValue = (i < gridSize * blockSizeInsert) ? i : 666;
            auto coord = sparseGrid.getCoord(i);
            auto expectedValue = coord.get(0);

//            std::cout << "sparseGrid(" << coord.get(0) << "," << coord.get(1) << "," << coord.get(2) << ") = "
//                      << sparseGrid.template get<0>(coord)
//                      << " == "
//                      << expectedValue
//                      << " == "
//                      << output.template get<0>(i) << " = output(" << i << ")"
//                      << std::endl;
            match &= output.template get<0>(i) == sparseGrid.template get<0>(coord);
            match &= output.template get<0>(i) == expectedValue;
        }

        BOOST_REQUIRE_EQUAL(match, true);
    }

    BOOST_AUTO_TEST_CASE(testTagBoundaries)
    {

        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
        typedef aggregate<float> AggregateT;

        dim3 gridSize(2, 2);
        dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize);

        BlockGeometry<dim, blockEdgeSize> blockGeometry(gridSize);
        SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(blockGeometry);

        sparseGrid.template setBackgroundValue<0>(666);
        mgpu::ofp_context_t ctx;

        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        dim3 pt1(0, 0, 0);
        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt1, 1);
        dim3 pt2(6, 6, 0);
        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt2, 1);
        dim3 pt3(7, 6, 0);
        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt3, 1);
        sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        dim3 pt4(8, 6, 0);
        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt4, 1);
        sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
        /////////
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        for (int y = 9; y <= 11; y++)
        {
            dim3 pt1(6, y, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt1, 1);
            dim3 pt2(7, y, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt2, 1);
        }
        sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        for (int y = 9; y <= 11; y++)
        {
            dim3 pt1(8, y, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt1, 1);
            dim3 pt2(9, y, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt2, 1);
        }
        sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

//        sparseGrid.hostToDevice(); //just sync masks
        sparseGrid.deviceToHost(); //just sync masks
//        sparseGrid.deviceToHost<0>();

        sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
        // Now tag the boundaries
        sparseGrid.tagBoundaries();

        // Get output
        openfpm::vector_gpu<AggregateT> output;
        output.resize(4 * 64);

        copyToOutputIfPadding<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
            auto coord = sparseGrid.getCoord(i);
            auto expectedValue =
                     ( i == 0
                    || i == 54
                    || i == 55
                    || i == 112
                    || i == 142 || i == 143 || i == 200 || i == 201 // (6,9), (7,9), (8,9), (9,9)
                    || i == 150 || i == 209 // (6,10), (9,10)
                    || i == 158 || i == 159 || i == 216 || i == 217 // (6,11), (7,11), (8,11), (9,11)
                     ) ? 1 : 0;

            std::cout
                    << "sparseGrid(" << coord.get(0) << "," << coord.get(1) << "," << coord.get(2) << ") = "
                    << sparseGrid.template get<0>(coord) << " | "
                    << expectedValue
                    << " == "
                    << output.template get<0>(i) << " = output(" << i << ")"
                    << std::endl;
            match &= output.template get<0>(i) == expectedValue;
        }

        BOOST_REQUIRE_EQUAL(match, true);
    }

        BOOST_AUTO_TEST_CASE(testTagBoundaries2)
    {
        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
        typedef aggregate<float> AggregateT;

        dim3 gridSize(2, 2);
        dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize);

        BlockGeometry<dim, blockEdgeSize> blockGeometry(gridSize);
        SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(blockGeometry);

        sparseGrid.template setBackgroundValue<0>(666);
        mgpu::ofp_context_t ctx;

        ///////
        {
            sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
            dim3 ptd1(6, 6, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd1, 1);
            dim3 ptd2(6, 7, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd2, 1);
            dim3 ptd3(7, 6, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd3, 1);
            dim3 ptd4(7, 7, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd4, 1);
            sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
        }
        {
            sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
            dim3 ptd1(8, 6, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd1, 1);
            dim3 ptd2(9, 6, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd2, 1);
            dim3 ptd3(8, 7, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd3, 1);
            dim3 ptd4(9, 7, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd4, 1);
            sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
        }
        {
            sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
            dim3 ptd1(6, 8, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd1, 1);
            dim3 ptd2(7, 8, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd2, 1);
            dim3 ptd3(6, 9, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd3, 1);
            dim3 ptd4(7, 9, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd4, 1);
            sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
        }
        {
            sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
            dim3 ptd1(8, 8, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd1, 1);
            dim3 ptd2(8, 9, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd2, 1);
            dim3 ptd3(9, 8, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd3, 1);
            dim3 ptd4(9, 9, 0);
            insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), ptd4, 1);
            sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
        }
        ///////

//        sparseGrid.hostToDevice(); //just sync masks
        sparseGrid.deviceToHost(); //just sync masks
//        sparseGrid.deviceToHost<0>();

        sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
        // Now tag the boundaries
        sparseGrid.tagBoundaries();

        // Get output
        openfpm::vector_gpu<AggregateT> output;
        output.resize(4 * 64);

        copyToOutputIfPadding<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
            auto coord = sparseGrid.getCoord(i);
            auto expectedValue =
                     (
                             i == 54 || i == 55 || i == 62 // (6,6), (7,6), (6,7)
                          || i == 134 || i == 142 || i == 143 // (6,8), (6,9), (7,9)
                          || i == 112 || i == 113 || i == 121 // (8,6), (9,6), (9,7)
                          || i == 200 || i == 193 || i == 201 // (8,9), (9,8), (9,9)
                     ) ? 1 : 0;

            std::cout
                    << "sparseGrid(" << coord.get(0) << "," << coord.get(1) << "," << coord.get(2) << ") = "
                    << sparseGrid.template get<0>(coord) << " | "
                    << expectedValue
                    << " == "
                    << output.template get<0>(i) << " = output(" << i << ")"
                    << std::endl;
            match &= output.template get<0>(i) == expectedValue;
        }

        BOOST_REQUIRE_EQUAL(match, true);
    }

    BOOST_AUTO_TEST_CASE(testStencilHeat)
    {
        printf("\n");

        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
        typedef aggregate<float> AggregateT;

        dim3 gridSize(2, 2);
        dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize);

        BlockGeometry<dim, blockEdgeSize> blockGeometry(gridSize);
        SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(blockGeometry);
        mgpu::ofp_context_t ctx;
        sparseGrid.template setBackgroundValue<0>(0);

        // Insert values on the grid
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        insertConstantValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), 0);
        sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        insertBoundaryValuesHeat<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel());
        sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

//        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
//        dim3 pt2(4, 4, 0);
//        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt2, 100);
//        sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

//        // Now tag the boundaries
        sparseGrid.tagBoundaries();

        // Now apply the laplacian operator
        const unsigned int maxIter = 1000;
//        const unsigned int maxIter = 10;
        for (unsigned int iter=0; iter<maxIter; ++iter)
        {
            sparseGrid.applyStencils<HeatStencil<dim, 0>>(STENCIL_MODE_INPLACE, 0.1);
            cudaDeviceSynchronize();
        }

        // Get output
        openfpm::vector_gpu<AggregateT> output;
        output.resize(4 * 64);

        copyBlocksToOutput<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
            auto coord = sparseGrid.getCoord(i);
            float expectedValue = 10.0 * coord.get(0) / (gridSize.x * blockEdgeSize - 1);

            std::cout
                    << "sparseGrid(" << coord.get(0) << "," << coord.get(1) << "," << coord.get(2) << ") = "
                    << sparseGrid.template get<0>(coord) << " | "
                    << expectedValue
                    << " == "
                    << output.template get<0>(i) << " = output(" << i << ")"
                    << std::endl;
            match &= fabs(output.template get<0>(i) - expectedValue) < 1e-2;

        }

        BOOST_REQUIRE_EQUAL(match, true);
//        BOOST_REQUIRE_CLOSE(output.template get<0>(255), 3.20309591e-05, 1e-6);
    }

    BOOST_AUTO_TEST_CASE(testStencilHeatInsert)
    {
        printf("\n");

        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
        typedef aggregate<float> AggregateT;

        dim3 gridSize(2, 2);
        dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize);

        BlockGeometry<dim, blockEdgeSize> blockGeometry(gridSize);
        SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(blockGeometry);
        mgpu::ofp_context_t ctx;
        sparseGrid.template setBackgroundValue<0>(0);

        // Insert values on the grid
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        insertConstantValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), 0);
        sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        insertBoundaryValuesHeat<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel());
        sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

//        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
//        dim3 pt2(4, 4, 0);
//        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt2, 100);
//        sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

        sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

//        // Now tag the boundaries
        sparseGrid.tagBoundaries();

        // Now apply the laplacian operator
        const unsigned int maxIter = 1000;
//        const unsigned int maxIter = 10;
        for (unsigned int iter=0; iter<maxIter; ++iter)
        {
            sparseGrid.applyStencils<HeatStencil<dim, 0>>(STENCIL_MODE_INSERT, 0.1);
//            sparseGrid.applyStencils<HeatStencil<dim, 0>>(STENCIL_MODE_INPLACE, 0.1);
            cudaDeviceSynchronize();
        }

        // Get output
        openfpm::vector_gpu<AggregateT> output;
        output.resize(4 * 64);

        copyBlocksToOutput<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
            grid_key_dx<dim, int> coord = sparseGrid.getCoord(i);
            float expectedValue = 10.0 * coord.get(0) / (gridSize.x * blockEdgeSize - 1);

            unsigned int check = sparseGrid.getLinId(coord);

            std::cout
                    << "invLinId=" << check << ", "
//                    << "sparseGrid(" << coord.get(0) << "," << coord.get(1) << "," << coord.get(2) << ") = "
                    << "sparseGrid(" << coord.get(0) << "," << coord.get(1) << ") = "
                    << sparseGrid.template get<0>(coord) << " | "
                    << expectedValue
                    << " == "
                    << output.template get<0>(i) << " = output(" << i << ")"
                    << std::endl;
            match &= fabs(output.template get<0>(i) - expectedValue) < 1e-2;

        }

        BOOST_REQUIRE_EQUAL(match, true);
//        BOOST_REQUIRE_CLOSE(output.template get<0>(255), 3.20309591e-05, 1e-6);
    }

BOOST_AUTO_TEST_SUITE_END()
