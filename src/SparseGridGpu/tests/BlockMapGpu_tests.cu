//
// Created by tommaso on 14/05/19.
//

#define BOOST_TEST_DYN_LINK

#include <hip/hip_runtime.h>
#include "config.h"
#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/BlockMapGpu.hpp"
#include "SparseGridGpu/BlockMapGpu_ker.cuh"
#include "SparseGridGpu/BlockMapGpu_kernels.cuh"
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


    int pos = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int dataBlockId = pos / BlockT::size;
    unsigned int offset = pos % BlockT::size;

    auto encap = sparseGrid.template insertBlock<chunksPerBlock>(dataBlockId,BlockT::size);


    encap.template get<p>()[offset] = pos;
    BlockMapGpu_ker<>::setExist(encap.template get<pMask>()[offset]);

    __syncthreads();

    sparseGrid.flush_block_insert();
}

template<unsigned int p, typename SparseGridType>
__global__ void insertValuesHalfBlock(SparseGridType sparseGrid)
{
    sparseGrid.init();

    int pos = blockIdx.x * blockDim.x + threadIdx.x;

    constexpr unsigned int dataChunkSize = BlockTypeOf<typename SparseGridType::AggregateType, p>::size;
    if (threadIdx.x % dataChunkSize < dataChunkSize/ 2)
    {
        sparseGrid.template insert<p>(pos) = pos;
    }

    __syncthreads();

    sparseGrid.flush_block_insert();
}

BOOST_AUTO_TEST_SUITE(BlockMapGpu_tests)

BOOST_AUTO_TEST_CASE(testBitwiseOps)
{
	BOOST_REQUIRE(BlockMapGpu_ker<>::getBit(1,0));
	BOOST_REQUIRE(!BlockMapGpu_ker<>::getBit(2,0));
	BOOST_REQUIRE(BlockMapGpu_ker<>::getBit(2,1));
	BOOST_REQUIRE(BlockMapGpu_ker<>::getBit(3,0));
	BOOST_REQUIRE(BlockMapGpu_ker<>::getBit(3,1));
	unsigned int m = 0;
	BOOST_REQUIRE(!BlockMapGpu_ker<>::getBit(m,0));
	BlockMapGpu_ker<>::setBit(m, 0);
	BOOST_REQUIRE(BlockMapGpu_ker<>::getBit(m,0));
	BlockMapGpu_ker<>::unsetBit(m, 0);
	BOOST_REQUIRE(!BlockMapGpu_ker<>::getBit(m,0));
}

BOOST_AUTO_TEST_CASE(testBackground)
{
	typedef aggregate<DataBlock<float, 64>> AggregateSGT;
	typedef aggregate<float> AggregateOutT;
	BlockMapGpu<AggregateSGT> sparseGrid;
	sparseGrid.template setBackgroundValue<0>(666);

	const unsigned int gridSize = 10;
	const unsigned int blockSize = 128;

	// Get output
	openfpm::vector_gpu<AggregateOutT> output;
	output.resize(gridSize * blockSize);
	hipLaunchKernelGGL(HIP_KERNEL_NAME(copyBlocksToOutput<0>), dim3(gridSize), dim3(blockSize), 0, 0, sparseGrid.toKernel(), output.toKernel());

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
	BlockMapGpu<AggregateT, 128> blockMap;
	blockMap.template setBackgroundValue<0>(666);

	const unsigned int gridSize = 3;
	const unsigned int bufferPoolSize = 128; // Should be multiple of BlockT::size
	const unsigned int blockSizeInsert = 128;
	const unsigned int gridSizeRead = gridSize + 1;
	const unsigned int blockSizeRead = 128;

	// Prealloc insert buffer
	blockMap.setGPUInsertBuffer(gridSize, bufferPoolSize);

	// Insert values
	hipLaunchKernelGGL(HIP_KERNEL_NAME(insertValues<0>), dim3(gridSize), dim3(blockSizeInsert), 0, 0, blockMap.toKernel());

	// Flush inserts
	mgpu::ofp_context_t ctx;
	blockMap.flush<smax_<0>>(ctx, flush_type::FLUSH_ON_DEVICE);

	// Get output
	openfpm::vector_gpu<AggregateOutT> output;
	output.resize(gridSizeRead * blockSizeRead);

	hipLaunchKernelGGL(HIP_KERNEL_NAME(copyBlocksToOutput<0>), dim3(gridSizeRead), dim3(blockSizeRead), 0, 0, blockMap.toKernel(), output.toKernel());

	output.template deviceToHost<0>();
	blockMap.template deviceToHost<0>();

	// Compare
	bool match = true;
	for (size_t i = 0; i < output.size(); i++)
	{
		auto expectedValue = (i < gridSize * blockSizeInsert) ? i : 666;
		std::cout << "blockMap(" << i << ") = " << blockMap.template get<0>(i)
				<< " == "
				<< expectedValue
				<< " == "
				<< output.template get<0>(i) << " = output(" << i << ")"
						<< std::endl;
		match &= output.template get<0>(i) == blockMap.template get<0>(i);
		match &= output.template get<0>(i) == expectedValue;
	}

	BOOST_REQUIRE_EQUAL(match, true);
}

BOOST_AUTO_TEST_CASE(testInsert_halfBlock) //todo
{
	typedef aggregate<DataBlock<float, 64>> AggregateT;
	typedef aggregate<float> AggregateOutT;
	BlockMapGpu<AggregateT, 128> blockMap;
	blockMap.template setBackgroundValue<0>(666);

	const unsigned int gridSize = 3;
	const unsigned int bufferPoolSize = 128; // Should be multiple of BlockT::size
	const unsigned int blockSizeInsert = 128;
	const unsigned int gridSizeRead = gridSize + 1;
	const unsigned int blockSizeRead = 128;

	// Prealloc insert buffer
	blockMap.setGPUInsertBuffer(gridSize, bufferPoolSize);

	// Insert values
	hipLaunchKernelGGL(HIP_KERNEL_NAME(insertValuesHalfBlock<0>), dim3(gridSize), dim3(blockSizeInsert), 0, 0, blockMap.toKernel());

	// Flush inserts
	mgpu::ofp_context_t ctx;
	blockMap.flush<smax_<0>>(ctx, flush_type::FLUSH_ON_DEVICE);

	// Get output
	openfpm::vector_gpu<AggregateOutT> output;
	output.resize(gridSizeRead * blockSizeRead);

	hipLaunchKernelGGL(HIP_KERNEL_NAME(copyBlocksToOutput<0>), dim3(gridSizeRead), dim3(blockSizeRead), 0, 0, blockMap.toKernel(), output.toKernel());

	output.template deviceToHost<0>();
	blockMap.template deviceToHost<0>();

	// Compare
	bool match = true;
	for (size_t i = 0; i < output.size(); i++)
	{
		auto expectedValue = (i < gridSize * blockSizeInsert) ? i : 666;
		constexpr unsigned int dataChunkSize = BlockTypeOf<AggregateT, 0>::size;
		int offset = i % dataChunkSize;
		if (! (offset < dataChunkSize / 2))
		{
			expectedValue = 666; // Just the first half of each block was inserted
		}
		std::cout << "blockMap(" << i << ") = " << blockMap.template get<0>(i)
				  << " == "
				  << expectedValue
				  << " == "
				  << output.template get<0>(i) << " = output(" << i << ")"
				  << std::endl;
		match &= output.template get<0>(i) == blockMap.template get<0>(i);
		match &= output.template get<0>(i) == expectedValue;
	}

	BOOST_REQUIRE_EQUAL(match, true);
}

BOOST_AUTO_TEST_CASE(testInsert_blocked)
{
	typedef aggregate<DataBlock<float, 64>> AggregateT;
	typedef aggregate<float> AggregateOutT;
	BlockMapGpu<AggregateT, 128> sparseGrid;
	sparseGrid.template setBackgroundValue<0>(666);

	const unsigned int gridSize = 3;
	const unsigned int bufferPoolSize = 4; // Should be multiple of BlockT::size
	const unsigned int blockSizeInsert = 128;
	const unsigned int gridSizeRead = gridSize + 1;
	const unsigned int blockSizeRead = 128;

////////// DEBUG
//        const unsigned int gridSize = 2;
//        const unsigned int bufferPoolSize = 128; // Should be multiple of BlockT::size
//        const unsigned int blockSizeInsert = 64;
//        const unsigned int gridSizeRead = gridSize + 1;
//        const unsigned int blockSizeRead = 64;
//////////

	// Prealloc insert buffer
	sparseGrid.setGPUInsertBuffer(gridSize, bufferPoolSize);

	// Insert values
	hipLaunchKernelGGL(HIP_KERNEL_NAME(insertValuesBlocked<0, 2>), dim3(gridSize), dim3(blockSizeInsert), 0, 0, sparseGrid.toKernel());

	// Flush inserts
	mgpu::ofp_context_t ctx;
	sparseGrid.flush<smax_<0>>(ctx, flush_type::FLUSH_ON_DEVICE);

	// Get output
	openfpm::vector_gpu<AggregateOutT> output;
	output.resize(gridSizeRead * blockSizeRead);

	hipLaunchKernelGGL(HIP_KERNEL_NAME(copyBlocksToOutput<0>), dim3(gridSizeRead), dim3(blockSizeRead), 0, 0, sparseGrid.toKernel(), output.toKernel());

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

BOOST_AUTO_TEST_SUITE_END()

