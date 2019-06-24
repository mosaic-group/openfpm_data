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
}

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

        dim3 gridSize(2, 2, 2);
        dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize, blockEdgeSize);

        BlockGeometry<dim, blockEdgeSize> blockGeometry(gridSize);
        SparseGridGpu<dim, blockEdgeSize, AggregateSGT> sparseGrid(blockGeometry);

        sparseGrid.template setBackground<0>(666);
//        const unsigned int gridSizeLin = 4;
//        const unsigned int bufferPoolSize = 1024;
//        const unsigned int blockSizeInsert = 128;
//        const unsigned int gridSizeRead = gridSize + 1;
//        const unsigned int blockSizeRead = 128;

//        sparseGrid.setGPUInsertBuffer(gridSizeLin, bufferPoolSize);
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        sparseGrid.initializeGPUInsertBuffer<0>();

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
        SparseGridGpu<dim, blockEdgeSize, AggregateSGT> sparseGrid(blockGeometry);

        sparseGrid.template setBackground<0>(666);

        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
        sparseGrid.initializeGPUInsertBuffer<0>();

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

BOOST_AUTO_TEST_SUITE_END()
