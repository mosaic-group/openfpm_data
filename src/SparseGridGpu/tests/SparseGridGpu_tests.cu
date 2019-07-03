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

template<unsigned int dim, unsigned int p>
struct LaplacianStencil
{
    // This is an example of a laplacian stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int supportRadius = 1;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            grid_key_dx<dim> & dataBlockCoord,
            unsigned int offset,
            grid_key_dx<dim> & pointCoord,
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
        for (unsigned int iter=0; iter<1000; ++iter)
        {
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
            __syncthreads();
        }
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
        sparseGrid.initializeGPUInsertBuffer();

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
        sparseGrid.initializeGPUInsertBuffer();

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
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        sparseGrid.initializeGPUInsertBuffer();

        // Insert only one point in the SparseGrid
//        grid_key_dx<dim, int> pt({6,6});
//        sparseGrid.template insert<0>(pt) = 1;
        dim3 pt1(0, 0, 0);
        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt1, 1);
        dim3 pt2(6, 6, 0);
        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt2, 1);
        dim3 pt3(7, 6, 0);
        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt3, 1);
        mgpu::ofp_context_t ctx;
        sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

//        sparseGrid.hostToDevice(); //just sync masks
        sparseGrid.deviceToHost(); //just sync masks
//        sparseGrid.deviceToHost<0>();

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
            auto expectedValue = (i == 0 || i == 54 || i == 55) ? 1 : 0;

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
        constexpr unsigned int dim = 2;
        constexpr unsigned int blockEdgeSize = 8;
        typedef aggregate<float> AggregateT;

        dim3 gridSize(2, 2);
        dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize);

        BlockGeometry<dim, blockEdgeSize> blockGeometry(gridSize);
        SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(blockGeometry);

        sparseGrid.template setBackgroundValue<0>(0);

        // Insert values on the grid
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        sparseGrid.initializeGPUInsertBuffer();

        insertConstantValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), 0);

        mgpu::ofp_context_t ctx;
        sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
        //
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        sparseGrid.initializeGPUInsertBuffer();

        dim3 pt2(4, 4, 0);
        insertOneValue<0> << < gridSize, blockSizeInsert >> > (sparseGrid.toKernel(), pt2, 100);

//        mgpu::ofp_context_t ctx;
        sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

//        // Now tag the boundaries
//        sparseGrid.tagBoundaries();

        // Now apply the laplacian operator
        sparseGrid.applyStencils<HeatStencil<dim, 0>>(STENCIL_MODE_INPLACE, 0.1);

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
            auto expectedValue = 0;

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
    }

BOOST_AUTO_TEST_SUITE_END()
