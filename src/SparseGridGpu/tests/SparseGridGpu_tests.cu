//
// Created by tommaso on 14/05/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/SparseGridGpu.hpp"
#include "SparseGridGpu/SparseGridGpu_ker.cuh"
#include "SparseGridGpu/SparseGridGpu_kernels.cuh"
#include "SparseGridGpu/DataBlock.cuh"

#include <limits>

template<unsigned int p, typename SparseGridType, typename VectorOutType>
__global__ void copyBlocksToOutput(SparseGridType sparseGrid, VectorOutType output)
{
    int pos = blockIdx.x * blockDim.x + threadIdx.x;
    output.template get<p>(pos) = sparseGrid.template get<p>(pos);
}

template<unsigned int p, typename SparseGridType>
__global__ void insertValues(SparseGridType sparseGrid)
{
    sparseGrid.init();

    int pos = blockIdx.x * blockDim.x + threadIdx.x;

    sparseGrid.template insert<p>(pos) = pos;

    __syncthreads();

    sparseGrid.flush_block_insert();
}

template<unsigned int p, unsigned int chunksPerBlock, typename SparseGridType>
__global__ void insertValuesBlocked(SparseGridType sparseGrid)
{
    constexpr unsigned int pMask = SparseGridType::pMask;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p> BlockT;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, pMask> MaskBlockT;

    sparseGrid.init();

    __shared__ BlockT *blocks[chunksPerBlock];
    __shared__ MaskBlockT *masks[chunksPerBlock];

    int pos = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int dataBlockId = pos / BlockT::size;
    unsigned int offset = pos % BlockT::size;

    unsigned int dataBlockNum = dataBlockId % chunksPerBlock;

    if (offset == 0) // Just one thread per data block
    {
        auto & encap = sparseGrid.insertBlock(dataBlockId);
        blocks[dataBlockNum] = &(encap.template get<p>());
        masks[dataBlockNum] = &(encap.template get<pMask>());
    }

    __syncthreads();

//    (*(blocks[dataBlockNum]))[offset] = pos;
    blocks[dataBlockNum]->block[offset] = pos;
    BlockT::setElement(masks[dataBlockNum]->block[offset]);

    sparseGrid.flush_block_insert();
}

BOOST_AUTO_TEST_SUITE(SparseGridGpu_tests)

    BOOST_AUTO_TEST_CASE(testBackground)
    {
        typedef aggregate<DataBlock<float, 64>> AggregateSGT;
        typedef aggregate<float> AggregateOutT;
        SparseGridGpu<AggregateSGT> sparseGrid;
        sparseGrid.template setBackground<0>(666);

        const unsigned int gridSize = 10;
        const unsigned int blockSize = 128;

        // Get output
        openfpm::vector_gpu<AggregateOutT> output;
        output.resize(gridSize * blockSize);
        copyBlocksToOutput<0> <<< gridSize, blockSize >>> (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
//            std::cout << "output(" << i << ") = " << output.template get<0>(i) << std::endl;
            match &= output.template get<0>(i) == sparseGrid.template get<0>(i);
        }

        BOOST_REQUIRE_EQUAL(match, true);
    }

    BOOST_AUTO_TEST_CASE(testInsert)
    {
        typedef aggregate<DataBlock<float, 64>> AggregateT;
        typedef aggregate<float> AggregateOutT;
        SparseGridGpu<AggregateT, 128> sparseGrid;
        sparseGrid.template setBackground<0>(666);

        const unsigned int gridSize = 10;
        const unsigned int blockSizeInit = 128; // Should be multiple of BlockT::size
        const unsigned int blockSizeInsert = 128;
        const unsigned int gridSizeRead = gridSize + 1;
        const unsigned int blockSizeRead = 128;

        // Prealloc insert buffer
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInit);

        // Initialize the insert buffer
        sparseGrid.initializeGPUInsertBuffer<0>();

        // Insert values
        insertValues<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel());
//        insertValuesBlocked<0, 2> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel());

        // Flush inserts
        mgpu::ofp_context_t ctx;
        sparseGrid.flush<smax_<0>>(ctx, flush_type::FLUSH_ON_DEVICE);

        // Get output
        openfpm::vector_gpu<AggregateOutT> output;
        output.resize(gridSizeRead * blockSizeRead);

        copyBlocksToOutput<0> <<< gridSizeRead, blockSizeRead >>> (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
            auto expectedValue = (i < gridSize * blockSizeInsert) ? i : 666;
            std::cout << "sparseGrid(" << i << ") = " << sparseGrid.template get<0>(i)
                    << " == "
                    << expectedValue
                    << " == "
                    << output.template get<0>(i) << " = output(" << i << ")"
                            << std::endl;
            match &= output.template get<0>(i) == sparseGrid.template get<0>(i);
            match &= output.template get<0>(i) == expectedValue;
        }

        BOOST_REQUIRE_EQUAL(match, true);
    }

//todo: Write a test that actually uses sparseGridGpu, to test the whole workflow

BOOST_AUTO_TEST_SUITE_END()

