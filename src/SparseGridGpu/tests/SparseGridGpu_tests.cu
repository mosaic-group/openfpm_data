//
// Created by tommaso on 10/06/19.
//

#define BOOST_TEST_DYN_LINK

#define DISABLE_MPI_WRITTERS

#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/SparseGridGpu.hpp"
#include "SparseGridGpu/tests/utils/SparseGridGpu_testKernels.cuh"
#include "SparseGridGpu/tests/utils/SparseGridGpu_util_test.cuh"

template<unsigned int p1 , unsigned int p2, unsigned int chunksPerBlock=1, typename SparseGridType, typename ScalarT>
__global__ void insertConstantValue2(SparseGridType sparseGrid, ScalarT value)
{
    constexpr unsigned int pMask = SparseGridType::pMask;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p1> BlockT;

    sparseGrid.init();

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    grid_key_dx<SparseGridType::d, size_t> coord({x, y, z});

    auto pos = sparseGrid.getLinId(coord);
    unsigned int dataBlockId = pos / BlockT::size;
    unsigned int offset = pos % BlockT::size;

    auto encap = sparseGrid.insertBlock(dataBlockId);
    encap.template get<p1>()[offset] = value;
    encap.template get<p2>()[offset] = value;
    BlockMapGpu_ker<>::setExist(encap.template get<pMask>()[offset]);

    __syncthreads();

    sparseGrid.flush_block_insert();

    // Compiler avoid warning
    x++;
    y++;
    z++;
}

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

    sparseGrid.template insert<p>(coord) = x;

    __syncthreads();

    sparseGrid.flush_block_insert();

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

	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateOutT, blockEdgeSize> sparseGrid(blockGeometry);

	sparseGrid.template setBackgroundValue<0>(666);
	sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);

	insertValues<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel());

	mgpu::ofp_context_t ctx;
	sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	// Get output
	openfpm::vector_gpu<AggregateOutT> output;
	output.resize(4 * 64);

	copyBlocksToOutput<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), output.toKernel());

	output.template deviceToHost<0>();
	sparseGrid.template deviceToHost<0>();

	// Compare
	bool match = true;
	for (size_t i = 0; i < output.size(); i++)
	{
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

	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateOutT, blockEdgeSize> sparseGrid(blockGeometry);

	sparseGrid.template setBackgroundValue<0>(666);

	sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);

	insertValues<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel());

	mgpu::ofp_context_t ctx;
	sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	// Get output
	openfpm::vector_gpu<AggregateOutT> output;
	output.resize(sparseGrid.dim3SizeToInt(gridSize) * 64);

	copyBlocksToOutput<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), output.toKernel());

	output.template deviceToHost<0>();
	sparseGrid.template deviceToHost<0>();

	// Compare
	bool match = true;
	for (size_t i = 0; i < output.size(); i++)
	{
		auto coord = sparseGrid.getCoord(i);
		auto expectedValue = coord.get(0);

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

	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(blockGeometry);

	sparseGrid.template setBackgroundValue<0>(666);
	mgpu::ofp_context_t ctx;

	sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
	dim3 pt1(0, 0, 0);
	CUDA_LAUNCH_DIM3((insertOneValue<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), pt1, 1);
	dim3 pt2(6, 6, 0);
	CUDA_LAUNCH_DIM3((insertOneValue<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), pt2, 1);
	dim3 pt3(7, 6, 0);
	CUDA_LAUNCH_DIM3((insertOneValue<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), pt3, 1);
	sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
	dim3 pt4(8, 6, 0);
	CUDA_LAUNCH_DIM3((insertOneValue<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), pt4, 1);
	sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
	/////////
	sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
	for (int y = 9; y <= 11; y++)
	{
		dim3 pt1(6, y, 0);
		CUDA_LAUNCH_DIM3((insertOneValue<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), pt1, 1);
		dim3 pt2(7, y, 0);
		CUDA_LAUNCH_DIM3((insertOneValue<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), pt2, 1);
	}
	sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
	for (int y = 9; y <= 11; y++)
	{
		dim3 pt1(8, y, 0);
		CUDA_LAUNCH_DIM3((insertOneValue<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), pt1, 1);
		dim3 pt2(9, y, 0);
		CUDA_LAUNCH_DIM3((insertOneValue<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), pt2, 1);
	}
	sparseGrid.template flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

//        sparseGrid.hostToDevice(); //just sync masks
	sparseGrid.deviceToHost(); //just sync masks
	sparseGrid.deviceToHost<0>();

	sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
	// Now tag the boundaries
	sparseGrid.tagBoundaries(ctx);

	// Get output
	openfpm::vector_gpu<AggregateT> output;
	output.resize(4 * 64);

	CUDA_LAUNCH_DIM3((copyToOutputIfPadding<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), output.toKernel());

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
				<< "sparseGrid(" << coord.get(0) << "," << coord.get(1) << ") = "
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

	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(blockGeometry);

	sparseGrid.template setBackgroundValue<0>(666);
	mgpu::ofp_context_t ctx;

	///////
	{
		sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
		dim3 ptd1(6, 6, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd1, 1);
		dim3 ptd2(6, 7, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd2, 1);
		dim3 ptd3(7, 6, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd3, 1);
		dim3 ptd4(7, 7, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd4, 1);
		sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
	}
	{
		sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
		dim3 ptd1(8, 6, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd1, 1);
		dim3 ptd2(9, 6, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd2, 1);
		dim3 ptd3(8, 7, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd3, 1);
		dim3 ptd4(9, 7, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd4, 1);
		sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
	}
	{
		sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
		dim3 ptd1(6, 8, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd1, 1);
		dim3 ptd2(7, 8, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd2, 1);
		dim3 ptd3(6, 9, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd3, 1);
		dim3 ptd4(7, 9, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd4, 1);
		sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
	}
	{
		sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
		dim3 ptd1(8, 8, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd1, 1);
		dim3 ptd2(8, 9, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd2, 1);
		dim3 ptd3(9, 8, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd3, 1);
		dim3 ptd4(9, 9, 0);
		insertOneValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), ptd4, 1);
		sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
	}
	///////

	sparseGrid.deviceToHost(); //just sync masks

	sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
	// Now tag the boundaries
	sparseGrid.tagBoundaries(ctx);

	// Get output
	openfpm::vector_gpu<AggregateT> output;
	output.resize(4 * 64);

	copyToOutputIfPadding<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), output.toKernel());

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
	typedef aggregate<float,float> AggregateT;

	dim3 gridSize(2, 2);
	dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize);

	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// Insert values on the grid
	sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
	insertConstantValue<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), 0);
	sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
	sparseGrid.tagBoundaries(ctx);

    cudaDeviceSynchronize();

    sparseGrid.template applyBoundaryStencils<BoundaryStencilSetXRescaled<dim,0,0>>(0.0 ,gridSize.x * blockEdgeSize, 0.0, 10.0);

	cudaDeviceSynchronize();

        // Now apply the laplacian operator
	const unsigned int maxIter = 1000;
//    const unsigned int maxIter = 100;
	for (unsigned int iter=0; iter<maxIter; ++iter)
	{
		sparseGrid.applyStencils<HeatStencil<dim, 0, 1>>(STENCIL_MODE_INPLACE, 0.1);
        cudaDeviceSynchronize();
        sparseGrid.applyStencils<HeatStencil<dim, 1, 0>>(STENCIL_MODE_INPLACE, 0.1);
        cudaDeviceSynchronize();
	}

	sparseGrid.deviceToHost<0,1>();

	// Get output
	openfpm::vector_gpu<AggregateT> output;
	output.resize(4 * 64);

	copyBlocksToOutput<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel(), output.toKernel());

	output.template deviceToHost<0>();
	sparseGrid.template deviceToHost<0>();

	// Compare
	bool match = true;
	for (size_t i = 0; i < output.size(); i++)
	{
		auto coord = sparseGrid.getCoord(i);
		float expectedValue = 10.0 * coord.get(0) / (gridSize.x * blockEdgeSize - 1);

		std::cout
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

BOOST_AUTO_TEST_CASE(testStencilHeatInsert)
{
	printf("\n");

	constexpr unsigned int dim = 2;
	constexpr unsigned int blockEdgeSize = 8;
	typedef aggregate<float,float> AggregateT;

	dim3 gridSize(2, 2);
	dim3 blockSizeInsert(blockEdgeSize, blockEdgeSize);

	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// Insert values on the grid
	sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
	CUDA_LAUNCH_DIM3((insertConstantValue2<0,1>),gridSize, blockSizeInsert,sparseGrid.toKernel(), 0);
	sparseGrid.flush < smax_< 0 >, smax_<1>> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInsert);
	CUDA_LAUNCH_DIM3((insertBoundaryValuesHeat<0>),gridSize, blockSizeInsert,sparseGrid.toKernel());
	sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

//        // Now tag the boundaries
	sparseGrid.tagBoundaries(ctx);

	sparseGrid.template deviceToHost<0,1>();

	// Now apply the laplacian operator
	const unsigned int maxIter = 1000;
//        const unsigned int maxIter = 10;
	for (unsigned int iter=0; iter<maxIter; ++iter)
	{
		sparseGrid.applyStencils<HeatStencil<dim, 0,0>>(STENCIL_MODE_INSERT, 0.1);
	}


	// Get output
	openfpm::vector_gpu<AggregateT> output;
	output.resize(4 * 64);

	CUDA_LAUNCH_DIM3((copyBlocksToOutput<0>),gridSize, blockSizeInsert,sparseGrid.toKernel(), output.toKernel());

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

template<typename sparsegrid_type>
__global__ void sparse_grid_get_test(sparsegrid_type sparseGrid, grid_key_dx<3> key, float * data)
{
	*data = sparseGrid.template get<0>(key);
}

BOOST_AUTO_TEST_CASE(testFlushInsert)
{
	printf("\n");

	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 4;
	typedef aggregate<float,float> AggregateT;


	size_t sz[] = {137,100,57};

	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64> sparseGrid(sz);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	sparseGrid.insertFlush<0>(grid_key_dx<3>({3,6,7})) = 2.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({13,16,17})) = 3.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({13,46,27})) = 4.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({36,63,11})) = 5.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({37,96,47})) = 6.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({130,56,37})) = 7.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({131,76,17})) = 8.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({36,86,27})) = 9.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({34,36,7})) = 10.0;

	////// Add points in the same blocks

	sparseGrid.insertFlush<0>(grid_key_dx<3>({4,6,7})) = 2.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({12,16,17})) = 3.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({12,46,27})) = 4.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({35,63,11})) = 5.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({38,96,47})) = 6.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({131,56,37})) = 7.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({132,76,17})) = 8.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({37,86,27})) = 9.0;
	sparseGrid.insertFlush<0>(grid_key_dx<3>({35,36,7})) = 10.0;

	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({3,6,7})),2.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({13,16,17})),3.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({13,46,27})),4.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({36,63,11})),5.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({37,96,47})),6.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({130,56,37})),7.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({131,76,17})),8.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({36,86,27})),9.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({34,36,7})),10.0);

	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({4,6,7})),2.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({12,16,17})),3.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({12,46,27})),4.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({35,63,11})),5.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({38,96,47})),6.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({131,56,37})),7.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({132,76,17})),8.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({37,86,27})),9.0);
	BOOST_REQUIRE_EQUAL(sparseGrid.get<0>(grid_key_dx<3>({35,36,7})),10.0);

	sparseGrid.template hostToDevice<0>();

	// Check on device I can get information

	CudaMemory mem;

	mem.allocate(sizeof(float));

	grid_key_dx<3> key({3,6,7});

	sparse_grid_get_test<<<1,1>>>(sparseGrid.toKernel(),key,(float *)mem.getDevicePointer());

	mem.deviceToHost();

	BOOST_REQUIRE_EQUAL(*(float *)mem.getPointer(),2.0);

	grid_key_dx<3> key2({131,76,17});

	sparse_grid_get_test<<<1,1>>>(sparseGrid.toKernel(),key2,(float *)mem.getDevicePointer());

	mem.deviceToHost();

	BOOST_REQUIRE_EQUAL(*(float *)mem.getPointer(),8.0);
}

struct conv_coeff
{
	float coeff[3][3][3];
};

template<unsigned int dim, unsigned int p_src, unsigned int p_dst>
struct Conv3x3x3
{
    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

	typedef NNFull<dim> stencil_type;

    static constexpr unsigned int supportRadius = 1;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            const unsigned int dataBlockId,
            openfpm::sparse_index<unsigned int> dataBlockIdPos,
            unsigned int offset,
            grid_key_dx<dim, int> & pointCoord,
            DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore,
            bool applyStencilHere,
            conv_coeff & cc)
    {
        typedef typename SparseGridT::AggregateBlockType AggregateT;
        typedef ScalarTypeOf<AggregateT, p_src> ScalarT;

        constexpr unsigned int enlargedBlockSize = IntPow<
                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

        __shared__ ScalarT enlargedBlock[enlargedBlockSize];

        sparseGrid.template loadGhostBlock<p_src>(dataBlockLoad,dataBlockIdPos,enlargedBlock);

        __syncthreads();

        if (applyStencilHere)
        {
            const auto coord = sparseGrid.getCoordInEnlargedBlock(offset);
            const auto linId = sparseGrid.getLinIdInEnlargedBlock(offset);
            ScalarT tot = 0.0;
            for (int i = 0; i < dim; ++i)
            {
                for (int j = 0; j < dim; ++j)
                {
                    for (int k = 0; k < dim; ++k)
                    {
                    	grid_key_dx<dim,int> key;

                    	key.set_d(0,i-1);
                    	key.set_d(1,j-1);
                    	key.set_d(2,k-1);

                    	auto nPlusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, key);
                    	tot += enlargedBlock[nPlusId] * cc.coeff[i][j][k];
                    }
                }
            }

            dataBlockStore.template get<p_dst>()[offset] = tot;
        }
    }

    template <typename SparseGridT, typename CtxT>
    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
    {
        sparseGrid.template flush <smax_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

template<unsigned int dim, unsigned int p_src, unsigned int p_dst>
struct Conv3x3x3_noshared
{
    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

	typedef NNFull<dim> stencil_type;

    static constexpr unsigned int supportRadius = 1;

    template<typename SparseGridT, typename DataBlockWrapperT>
    static inline __device__ void stencil(
            SparseGridT & sparseGrid,
            const unsigned int dataBlockId,
            openfpm::sparse_index<unsigned int> dataBlockIdPos,
            unsigned int offset,
            grid_key_dx<dim, int> & pointCoord,
            DataBlockWrapperT & dataBlockLoad,
            DataBlockWrapperT & dataBlockStore,
            bool applyStencilHere,
            conv_coeff & cc)
    {
        typedef typename SparseGridT::AggregateBlockType AggregateT;
        typedef ScalarTypeOf<AggregateT, p_src> ScalarT;

        __syncthreads();

        __shared__ block_offset<int> pos[BLOCK_SIZE_STENCIL];

        if (applyStencilHere)
        {
            ScalarT tot = 0.0;
            for (int i = 0; i < dim; ++i)
            {
                for (int j = 0; j < dim; ++j)
                {
                    for (int k = 0; k < dim; ++k)
                    {

                    	grid_key_dx<dim,int> key;

                    	key.set_d(0,k-1);
                    	key.set_d(1,j-1);
                    	key.set_d(2,i-1);

                    	pos[threadIdx.x] = sparseGrid.template getNNPoint<stencil_type>(dataBlockIdPos, offset, key);

                    	tot += sparseGrid.template get<p_src>(pos[threadIdx.x]) * cc.coeff[i][j][k];
                    }
                }
            }

            dataBlockStore.template get<p_dst>()[offset] = tot;
        }
    }

    template <typename SparseGridT, typename CtxT>
    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
    {
        sparseGrid.template flush <smax_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

template<typename SparseGridZ>
void test_convolution_3x3x3()
{
	size_t sz[] = {1000,1000,1000};

	SparseGridZ sparseGrid(sz);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// now create 3 3D sphere

    grid_key_dx<3,int> start({256,256,256});

    dim3 gridSize(32,32,32);

    // Insert values on the grid
    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
    CUDA_LAUNCH_DIM3((insertSphere3D_radius<0>),
            gridSize, dim3(SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_,1,1),
            sparseGrid.toKernel(), start,64, 56, 1);

    sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    sparseGrid.template findNeighbours<NNFull<3>>(); // Pre-compute the neighbours pos for each block!

    sparseGrid.template setNNType<NNFull<3>>();
    sparseGrid.template tagBoundaries<NNFull<3>>(ctx);

    conv_coeff cc;

    for (int i = 0 ; i < 3 ; i++)
    {
    	for (int j = 0 ; j < 3 ; j++)
    	{
    		for (int k = 0 ; k < 3 ; k++)
    		{
    			cc.coeff[k][j][i] = 1.0;
    		}
    	}
    }


    sparseGrid.template applyStencils<Conv3x3x3<3,0,1>>(STENCIL_MODE_INPLACE,cc);

    sparseGrid.template deviceToHost<0,1>();

	auto & bm = sparseGrid.private_get_blockMap();
	auto & dataVector = bm.getDataBuffer();

	bool match = true;

    for (size_t i = 0 ; i < dataVector.size() ; i++)
    {
        for (size_t j = 0 ; j < 64 ; j++)
        {
			if (dataVector.template get<2>(i)[j] == 1)
			{
				match &= dataVector.template get<0>(i)[j]*27 == dataVector.template get<1>(i)[j];
			}
        }
    }

    BOOST_REQUIRE_EQUAL(match,true);
}

template<typename SparseGridZ>
void test_convolution_3x3x3_no_shared()
{
	size_t sz[] = {1000,1000,1000};

	SparseGridZ sparseGrid(sz);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// now create 3 3D sphere

    grid_key_dx<3,int> start({256,256,256});

    dim3 gridSize(32,32,32);

    // Insert values on the grid
    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
    CUDA_LAUNCH_DIM3((insertSphere3D_radius<0>),
            gridSize, dim3(SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_,1,1),
            sparseGrid.toKernel(), start,64, 56, 1);

    sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    sparseGrid.template findNeighbours<NNFull<3>>(); // Pre-compute the neighbours pos for each block!

    sparseGrid.template setNNType<NNFull<SparseGridZ::dims>>();
    sparseGrid.template tagBoundaries<NNFull<3>>(ctx,tag_boundaries::CALCULATE_EXISTING_POINTS);

    conv_coeff cc;

    for (int i = 0 ; i < 3 ; i++)
    {
    	for (int j = 0 ; j < 3 ; j++)
    	{
    		for (int k = 0 ; k < 3 ; k++)
    		{
    			cc.coeff[k][j][i] = 1.0;
    		}
    	}
    }

    sparseGrid.template applyStencils<Conv3x3x3_noshared<SparseGridZ::dims,0,1>>(STENCIL_MODE_INPLACE_NO_SHARED,cc);

    sparseGrid.template deviceToHost<0,1>();

	auto & bm = sparseGrid.private_get_blockMap();
	auto & dataVector = bm.getDataBuffer();

	bool match = true;

    for (size_t i = 0 ; i < dataVector.size() ; i++)
    {
        for (size_t j = 0 ; j < 64 ; j++)
        {
			if (dataVector.template get<2>(i)[j] == 1)
			{
				match &= dataVector.template get<0>(i)[j]*27 == dataVector.template get<1>(i)[j];
			}
        }
    }

    BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE(test3x3x3convolution_no_shared)
{
	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 4;
	typedef aggregate<float,float> AggregateT;

	test_convolution_3x3x3_no_shared<SparseGridGpu<dim, AggregateT, blockEdgeSize, 64, long int>>();
}

BOOST_AUTO_TEST_CASE(test3x3x3convolution_no_shared_z_morton)
{
	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 4;
	typedef aggregate<float,float> AggregateT;

	test_convolution_3x3x3_no_shared<SparseGridGpu_z<dim, AggregateT, blockEdgeSize, 64, long int>>();
}

BOOST_AUTO_TEST_CASE(test3x3x3convolution)
{
	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 4;
	typedef aggregate<float,float> AggregateT;

	test_convolution_3x3x3<SparseGridGpu<dim, AggregateT, blockEdgeSize, 64, long int>>();
}

BOOST_AUTO_TEST_CASE(test3x3x3convolution_morton_z)
{
	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 4;
	typedef aggregate<float,float> AggregateT;

	test_convolution_3x3x3<SparseGridGpu_z<dim, AggregateT, blockEdgeSize, 64, long int>>();
}

BOOST_AUTO_TEST_CASE(test_sparse_grid_iterator_sub_host)
{
	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 4;
	typedef aggregate<float, float> AggregateT;

	size_t sz[3] = {768,768,768};
	dim3 gridSize(32,32,32);

	grid_smb<dim, blockEdgeSize> blockGeometry(sz);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64, long int> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	///// Insert sparse content, a set of 3 hollow spheres /////
	// Sphere 1
	grid_key_dx<3,int> start1({256,256,256});
	sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
	CUDA_LAUNCH_DIM3((insertSphere3D<0>),
					 gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
					 sparseGrid.toKernel(), start1, 32, 0, 1);

	sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.template deviceToHost<0>();

	bool match = true;

	int count = 0;

	grid_key_dx<3> start({303,303,303});
	grid_key_dx<3> stop({337,337,337});

	auto it = sparseGrid.getIterator(start,stop);

	while (it.isNext())
	{
		auto key = it.get();

		match &= sparseGrid.template get<0>(key) == 1.0;

		sparseGrid.template get<0>(key) = 5.0;

		count++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(match,true);
	BOOST_REQUIRE_EQUAL(count,42875);
}




BOOST_AUTO_TEST_CASE(test_sparse_grid_iterator_host)
{
	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 4;
	typedef aggregate<float, float> AggregateT;

	size_t sz[3] = {512,512,512};
	dim3 gridSize(32,32,32);

	grid_smb<dim, blockEdgeSize> blockGeometry(sz);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64, long int> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	///// Insert sparse content, a set of 3 hollow spheres /////
	// Sphere 1
	grid_key_dx<3,int> start1({256,256,256});
	sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
	CUDA_LAUNCH_DIM3((insertSphere3D<0>),
					 gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
					 sparseGrid.toKernel(), start1, 64, 32, 1);

	sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.template deviceToHost<0>();

	bool match = true;

	int count = 0;

	auto it = sparseGrid.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		match &= sparseGrid.template get<0>(key) == 1.0;
		//unsigned char bl = sparseGrid.template get<2>(key);

		count++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(sparseGrid.countExistingElements(),count);
	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE(test_pack_request)
{
	size_t sz[] = {1000,1000,1000};

	constexpr int blockEdgeSize = 4;
	constexpr int dim = 3;

	typedef SparseGridGpu<dim, aggregate<float>, blockEdgeSize, 64, long int> SparseGridZ;

	SparseGridZ sparseGrid(sz);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// now create a 3D sphere

    grid_key_dx<3,int> start({256,256,256});

    dim3 gridSize(32,32,32);

    // Insert values on the grid
    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
    CUDA_LAUNCH_DIM3((insertSphere3D_radius<0>),
            gridSize, dim3(SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_,1,1),
            sparseGrid.toKernel(), start,64, 56, 1);

    sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
    sparseGrid.template deviceToHost<0>();

    size_t cnt = sparseGrid.countExistingElements();

    size_t req = 0;
    sparseGrid.packRequest<0>(req,ctx);

    size_t tot = 8 +                // how many chunks
    		     sparseGrid.private_get_index_array().size()*16 + 8 +// store the scan + chunk indexes
    		     cnt*(sizeof(float) + 2); // how much data

    BOOST_REQUIRE_EQUAL(req,tot);
}

BOOST_AUTO_TEST_CASE(test_pack_request_with_iterator)
{
	size_t sz[] = {1000,1000,1000};

	constexpr int blockEdgeSize = 4;
	constexpr int dim = 3;

	typedef SparseGridGpu<dim, aggregate<float>, blockEdgeSize, 64, long int> SparseGridZ;

	SparseGridZ sparseGrid(sz);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// now create a 3D sphere

    grid_key_dx<3,int> start({256,256,256});

    dim3 gridSize(32,32,32);

    // Insert values on the grid
    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
    CUDA_LAUNCH_DIM3((insertSphere3D_radius<0>),
            gridSize, dim3(SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_,1,1),
            sparseGrid.toKernel(), start,64, 56, 1);

    sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    size_t req = 0;
    sparseGrid.packReset();

    {
    grid_key_dx<3> start1({0,0,0});
    grid_key_dx<3> stop1({321,999,999});

    grid_key_dx<3> start2({322,0,0});
    grid_key_dx<3> stop2({999,999,999});

    auto it1 = sparseGrid.getIterator(start1,stop1);
    sparseGrid.template packRequest<0>(it1,req);

    auto it2 = sparseGrid.getIterator(start2,stop2);
    sparseGrid.template packRequest<0>(it2,req);

    sparseGrid.template packCalculate<0>(req,ctx);
    }

    sparseGrid.template deviceToHost<0>();


    size_t cnt = sparseGrid.countExistingElements();

    size_t tot = 8 +                // how many chunks
    		     sparseGrid.private_get_index_array().size()*16 + 8 +// store the scan + chunk indexes
    		     cnt*(sizeof(float) + 2); // how much data

    std::cout << __FILE__ << ":" << __LINE__ << "  To fix this" << std::endl;
//    BOOST_REQUIRE_EQUAL(req,tot);

    ////////////////////////////////// test something else

    req = 0;
    sparseGrid.packReset();

    {
    grid_key_dx<3> start1({0,0,0});
    grid_key_dx<3> stop1({999,999,999});

    auto it1 = sparseGrid.getIterator(start1,stop1);
    sparseGrid.template packRequest<0>(it1,req);

    auto it2 = sparseGrid.getIterator(start1,stop1);
    sparseGrid.template packRequest<0>(it2,req);

    sparseGrid.template packCalculate<0>(req,ctx);
    }


    tot = 8 +                // how many chunks
    		     sparseGrid.private_get_index_array().size()*16 + 8 + // store the scan + chunk indexes
    		     2*cnt*(sizeof(float) + 2); // how much data

//    BOOST_REQUIRE_EQUAL(req,tot);
}

BOOST_AUTO_TEST_CASE(sparsegridgpu_remove_test)
{
	size_t sz[] = {1000,1000,1000};

	constexpr int blockEdgeSize = 4;
	constexpr int dim = 3;

	typedef SparseGridGpu<dim, aggregate<float>, blockEdgeSize, 64, long int> SparseGridZ;

	SparseGridZ sparseGrid(sz);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// now create a 3D sphere

    grid_key_dx<3,int> start({256,256,256});

    dim3 gridSize(32,32,32);

    // Insert values on the grid
    sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
    CUDA_LAUNCH_DIM3((insertSphere3D_radius<0>),
            gridSize, dim3(SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_,1,1),
            sparseGrid.toKernel(), start,64, 56, 1);

    sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    // remove the center

    Box<3,unsigned int> remove_section1({310,0,0},{330,999,999});
    Box<3,unsigned int> remove_section2({0,310,0},{999,330,999});
    Box<3,unsigned int> remove_section3({0,0,310},{999,999,330});

    sparseGrid.remove(remove_section1);
    sparseGrid.remove(remove_section2);
    sparseGrid.remove(remove_section3);

    sparseGrid.removeAddUnpackFinalize<>(ctx);

    sparseGrid.deviceToHost<0>();

    // Check we have the sphere points with the exception of the box

    auto it = sparseGrid.getIterator();

    bool match = true;

    while (it.isNext())
    {
    	auto p = it.get();

    	Point<3,size_t> pt = p.toPoint();

    	// Calculate redius
    	float radius = sqrt((pt.get(0) - 320)*(pt.get(0) - 320) +
    						(pt.get(1) - 320)*(pt.get(1) - 320) +
    						(pt.get(2) - 320)*(pt.get(2) - 320));


    	if (radius < 55.99 || radius > 64.01)
    	{match &= false;}

    	if (remove_section1.isInside(pt) == true)
    	{match &= false;}

    	if (remove_section2.isInside(pt) == true)
    	{match &= false;}

    	if (remove_section3.isInside(pt) == true)
    	{match &= false;}

    	++it;
    }

    BOOST_REQUIRE_EQUAL(match,true);
}

template<typename SG_type>
void pack_unpack_test(SG_type & sparseGridDst, SG_type & sparseGridSrc,
		Box<3,size_t> & box1_dst,
		Box<3,size_t> & box2_dst,
		Box<3,size_t> & box3_dst,
		Box<3,size_t> & box4_dst,
		mgpu::ofp_context_t & ctx,
		bool test_pack)
{
    Box<3,size_t> box1_src({256,256,256},{273,390,390});
    Box<3,size_t> box2_src({320,256,256},{337,390,390});

    // And two vertical sections

    Box<3,size_t> box3_src({256,256,256},{273,320,390});
    Box<3,size_t> box4_src({320,256,256},{337,320,390});

    // Now we calculate the memory required to pack

    sparseGridSrc.packReset();

    size_t req = 0;
    auto sub_it = sparseGridSrc.getIterator(box1_src.getKP1(),box1_src.getKP2());
    sparseGridSrc.template packRequest<0,1>(sub_it,req);

    sub_it = sparseGridSrc.getIterator(box2_src.getKP1(),box2_src.getKP2());
    sparseGridSrc.template packRequest<0,1>(sub_it,req);

    sub_it = sparseGridSrc.getIterator(box3_src.getKP1(),box3_src.getKP2());
    sparseGridSrc.template packRequest<0,1>(sub_it,req);

    sub_it = sparseGridSrc.getIterator(box4_src.getKP1(),box4_src.getKP2());
    sparseGridSrc.template packRequest<0,1>(sub_it,req);

    sparseGridSrc.template packCalculate<0,1>(req,ctx);

    CudaMemory mem;
    mem.resize(req);

	// Create an object of preallocated memory for properties
	ExtPreAlloc<CudaMemory> & prAlloc_prp = *(new ExtPreAlloc<CudaMemory>(req,mem));

	prAlloc_prp.incRef();

	// Pack information
	Pack_stat sts;

    sub_it = sparseGridSrc.getIterator(box1_src.getKP1(),box1_src.getKP2());
    sparseGridSrc.template pack<0,1>(prAlloc_prp,sub_it,sts);

    sub_it = sparseGridSrc.getIterator(box2_src.getKP1(),box2_src.getKP2());
    sparseGridSrc.template pack<0,1>(prAlloc_prp,sub_it,sts);

    sub_it = sparseGridSrc.getIterator(box3_src.getKP1(),box3_src.getKP2());
    sparseGridSrc.template pack<0,1>(prAlloc_prp,sub_it,sts);

    sub_it = sparseGridSrc.getIterator(box4_src.getKP1(),box4_src.getKP2());
    sparseGridSrc.template pack<0,1>(prAlloc_prp,sub_it,sts);


	sparseGridSrc.template packFinalize<0,1>(prAlloc_prp,sts);

	// Now we analyze the package

	if (test_pack == true)
	{
		size_t ncnk = *(size_t *)mem.getPointer();
		BOOST_REQUIRE_EQUAL(ncnk,1107);
		size_t actual_offset = ncnk*sizeof(size_t) + sizeof(size_t) + 2*3*sizeof(int);
		mem.deviceToHost(actual_offset + ncnk*sizeof(unsigned int),actual_offset + ncnk*sizeof(unsigned int) + sizeof(unsigned int));
		unsigned int n_pnt = *(unsigned int *)((unsigned char *)mem.getPointer() + actual_offset + ncnk*sizeof(unsigned int));
		BOOST_REQUIRE_EQUAL(n_pnt,41003);

		actual_offset += align_number(sizeof(size_t),(ncnk+1)*sizeof(unsigned int));
		actual_offset += align_number(sizeof(size_t),n_pnt*(16));
		actual_offset += align_number(sizeof(size_t),n_pnt*sizeof(short int));


		ncnk = *(size_t *)((unsigned char *)mem.getPointer() + actual_offset);
		BOOST_REQUIRE_EQUAL(ncnk,1420);
		actual_offset += ncnk*sizeof(size_t) + sizeof(size_t) + 2*3*sizeof(int);
		mem.deviceToHost(actual_offset + ncnk*sizeof(unsigned int),actual_offset + ncnk*sizeof(unsigned int) + sizeof(unsigned int));
		n_pnt = *(unsigned int *)((unsigned char *)mem.getPointer() + actual_offset + ncnk*sizeof(unsigned int));
		BOOST_REQUIRE_EQUAL(n_pnt,54276);

		actual_offset += align_number(sizeof(size_t),(ncnk+1)*sizeof(unsigned int));
		actual_offset += align_number(sizeof(size_t),n_pnt*(16));
		actual_offset += align_number(sizeof(size_t),n_pnt*sizeof(short int));

		ncnk = *(size_t *)((unsigned char *)mem.getPointer() + actual_offset);
		BOOST_REQUIRE_EQUAL(ncnk,610);
		actual_offset += ncnk*sizeof(size_t) + sizeof(size_t) + 2*3*sizeof(int);
		mem.deviceToHost(actual_offset + ncnk*sizeof(unsigned int),actual_offset + ncnk*sizeof(unsigned int) + sizeof(unsigned int));
		n_pnt = *(unsigned int *)((unsigned char *)mem.getPointer() + actual_offset + ncnk*sizeof(unsigned int));
		BOOST_REQUIRE_EQUAL(n_pnt,20828);

		actual_offset += align_number(sizeof(size_t),(ncnk+1)*sizeof(unsigned int));
		actual_offset += align_number(sizeof(size_t),n_pnt*(16));
		actual_offset += align_number(sizeof(size_t),n_pnt*sizeof(short int));

		ncnk = *(size_t *)((unsigned char *)mem.getPointer() + actual_offset);
		BOOST_REQUIRE_EQUAL(ncnk,739);
		actual_offset += ncnk*sizeof(size_t) + sizeof(size_t) + 2*3*sizeof(int);
		mem.deviceToHost(actual_offset + ncnk*sizeof(unsigned int),actual_offset + ncnk*sizeof(unsigned int) + sizeof(unsigned int));
		n_pnt = *(unsigned int *)((unsigned char *)mem.getPointer() + actual_offset + ncnk*sizeof(unsigned int));
		BOOST_REQUIRE_EQUAL(n_pnt,27283);
	}

	prAlloc_prp.reset();

	Unpack_stat ps;

	sparseGridDst.removeAddUnpackReset();

	// sub-grid where to unpack
	auto sub2 = sparseGridDst.getIterator(box1_dst.getKP1(),box1_dst.getKP2());
	sparseGridDst.remove(box1_dst);
	sparseGridDst.template unpack<0,1>(prAlloc_prp,sub2,ps,ctx);

	sub2 = sparseGridDst.getIterator(box2_dst.getKP1(),box2_dst.getKP2());
	sparseGridDst.remove(box2_dst);
	sparseGridDst.template unpack<0,1>(prAlloc_prp,sub2,ps,ctx);

	sub2 = sparseGridDst.getIterator(box3_dst.getKP1(),box3_dst.getKP2());
	sparseGridDst.remove(box3_dst);
	sparseGridDst.template unpack<0,1>(prAlloc_prp,sub2,ps,ctx);

	sub2 = sparseGridDst.getIterator(box4_dst.getKP1(),box4_dst.getKP2());
	sparseGridDst.remove(box4_dst);
	sparseGridDst.template unpack<0,1>(prAlloc_prp,sub2,ps,ctx);

	sparseGridDst.template removeAddUnpackFinalize<0,1>(ctx);

	sparseGridDst.template deviceToHost<0,1>();
}

BOOST_AUTO_TEST_CASE(sparsegridgpu_pack_unpack)
{
	size_t sz[] = {1000,1000,1000};

	constexpr int blockEdgeSize = 4;
	constexpr int dim = 3;

    Box<3,size_t> box1_dst({256,256,256},{273,390,390});
    Box<3,size_t> box2_dst({274,256,256},{291,390,390});

    Box<3,size_t> box3_dst({300,256,256},{317,320,390});
    Box<3,size_t> box4_dst({320,256,256},{337,320,390});

	typedef SparseGridGpu<dim, aggregate<float,float[3]>, blockEdgeSize, 64, long int> SparseGridZ;

	SparseGridZ sparseGridSrc(sz);
	SparseGridZ sparseGridDst(sz);
	mgpu::ofp_context_t ctx;
	sparseGridSrc.template setBackgroundValue<0>(0);
	sparseGridDst.template setBackgroundValue<0>(0);

	// now create a 3D sphere

    grid_key_dx<3,int> start({256,256,256});

    dim3 gridSize(32,32,32);

    // Insert values on the grid
    sparseGridSrc.setGPUInsertBuffer(gridSize,dim3(1));
    CUDA_LAUNCH_DIM3((insertSphere3D_radiusV<0>),
            gridSize, dim3(SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_,1,1),
            sparseGridSrc.toKernel(), start,64, 56, 1);

    sparseGridSrc.flush < smax_< 0 >, smax_<1> > (ctx, flush_type::FLUSH_ON_DEVICE);

    // Now we pack two vertical sections

	pack_unpack_test(sparseGridDst,sparseGridSrc,
					 box1_dst,box2_dst,
					 box3_dst,box4_dst,
					ctx,true);

	sparseGridDst.template deviceToHost<0,1>();

	int cnt1 = 0;
	int cnt2 = 0;
	int cnt3 = 0;
	int cnt4 = 0;

	auto it = sparseGridDst.getIterator();

	bool match = true;

	while (it.isNext())
	{
		auto p = it.get();

		auto pt = p.toPoint();

		if (box1_dst.isInside(pt) == true)
		{
			++cnt1;

		    const long int x = (long int)pt.get(0) - (start.get(0) + gridSize.x / 2 * blockEdgeSize);
		    const long int y = (long int)pt.get(1) - (start.get(1) + gridSize.y / 2 * blockEdgeSize);
		    const long int z = (long int)pt.get(2) - (start.get(2) + gridSize.z / 2 * blockEdgeSize);

		    float radius = sqrt((float) (x*x + y*y + z*z));

		    bool is_active = radius < 64 && radius > 56;

		    if (is_active == true)
		    {match &= true;}
		    else
		    {match &= false;}
		}
		else if (box2_dst.isInside(pt) == true)
		{
			++cnt2;

		    const long int x = (long int)pt.get(0) - (start.get(0) - 46 + gridSize.x / 2 * blockEdgeSize);
		    const long int y = (long int)pt.get(1) - (start.get(1) + gridSize.y / 2 * blockEdgeSize);
		    const long int z = (long int)pt.get(2) - (start.get(2) + gridSize.z / 2 * blockEdgeSize);

		    float radius = sqrt((float) (x*x + y*y + z*z));

		    bool is_active = radius < 64 && radius > 56;

		    if (is_active == true)
		    {match &= true;}
		    else
		    {match &= false;}
		}
		else if (box3_dst.isInside(pt) == true)
		{
			++cnt3;

		    const long int x = (long int)pt.get(0) - (start.get(0) + 44 + gridSize.x / 2 * blockEdgeSize);
		    const long int y = (long int)pt.get(1) - (start.get(1) + gridSize.y / 2 * blockEdgeSize);
		    const long int z = (long int)pt.get(2) - (start.get(2) + gridSize.z / 2 * blockEdgeSize);

		    float radius = sqrt((float) (x*x + y*y + z*z));

		    bool is_active = radius < 64 && radius > 56;

		    if (is_active == true)
		    {match &= true;}
		    else
		    {match &= false;}
		}
		else if (box4_dst.isInside(pt) == true)
		{
			++cnt4;

		    const long int x = (long int)pt.get(0) - (start.get(0) + gridSize.x / 2 * blockEdgeSize);
		    const long int y = (long int)pt.get(1) - (start.get(1) + gridSize.y / 2 * blockEdgeSize);
		    const long int z = (long int)pt.get(2) - (start.get(2) + gridSize.z / 2 * blockEdgeSize);

		    float radius = sqrt((float) (x*x + y*y + z*z));

		    bool is_active = radius < 64 && radius > 56;

		    if (is_active == true)
		    {match &= true;}
		    else
		    {match &= false;}
		}

		++it;
	}

	BOOST_REQUIRE_EQUAL(match,true);
	BOOST_REQUIRE_EQUAL(cnt1,41003);
	BOOST_REQUIRE_EQUAL(cnt2,54276);
	BOOST_REQUIRE_EQUAL(cnt3,20828);
	BOOST_REQUIRE_EQUAL(cnt4,27283);

	// Now we remove even points

    // Insert values on the grid
    sparseGridSrc.setGPUInsertBuffer(gridSize,dim3(1));
    CUDA_LAUNCH_DIM3((removeSphere3D_even_radiusV<0>),
            gridSize, dim3(SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_*SparseGridZ::blockEdgeSize_,1,1),
            sparseGridSrc.toKernel(), start,64, 56, 1);

    pack_unpack_test(sparseGridDst,sparseGridSrc,
			 	 	 box1_dst,box2_dst,
			 	 	 box3_dst,box4_dst,
    				ctx,false);

	sparseGridDst.template deviceToHost<0,1>();

	cnt1 = 0;
	cnt2 = 0;
	cnt3 = 0;
	cnt4 = 0;

	auto it2 = sparseGridDst.getIterator();

	match = true;

	while (it2.isNext())
	{
		auto p = it2.get();

		auto pt = p.toPoint();

		if (box1_dst.isInside(pt) == true)
		{
			++cnt1;

			const long int x = (long int)pt.get(0) - (start.get(0) + gridSize.x / 2 * blockEdgeSize);
			const long int y = (long int)pt.get(1) - (start.get(1) + gridSize.y / 2 * blockEdgeSize);
			const long int z = (long int)pt.get(2) - (start.get(2) + gridSize.z / 2 * blockEdgeSize);

			float radius = sqrt((float) (x*x + y*y + z*z));

			bool is_active = radius < 64 && radius > 56;

			if (is_active == true)
			{match &= true;}
			else
			{match &= false;}
		}
		else if (box2_dst.isInside(pt) == true)
		{
			++cnt2;

			const long int x = (long int)pt.get(0) - (start.get(0) - 46 + gridSize.x / 2 * blockEdgeSize);
			const long int y = (long int)pt.get(1) - (start.get(1) + gridSize.y / 2 * blockEdgeSize);
			const long int z = (long int)pt.get(2) - (start.get(2) + gridSize.z / 2 * blockEdgeSize);

			float radius = sqrt((float) (x*x + y*y + z*z));

			bool is_active = radius < 64 && radius > 56;

			if (is_active == true)
			{match &= true;}
			else
			{match &= false;}
		}
		else if (box3_dst.isInside(pt) == true)
		{
			++cnt3;

			const long int x = (long int)pt.get(0) - (start.get(0) + 44 + gridSize.x / 2 * blockEdgeSize);
			const long int y = (long int)pt.get(1) - (start.get(1) + gridSize.y / 2 * blockEdgeSize);
			const long int z = (long int)pt.get(2) - (start.get(2) + gridSize.z / 2 * blockEdgeSize);

			float radius = sqrt((float) (x*x + y*y + z*z));

			bool is_active = radius < 64 && radius > 56;

			if (is_active == true)
			{match &= true;}
			else
			{match &= false;}
		}
		else if (box4_dst.isInside(pt) == true)
		{
			++cnt4;

			const long int x = (long int)pt.get(0) - (start.get(0) + gridSize.x / 2 * blockEdgeSize);
			const long int y = (long int)pt.get(1) - (start.get(1) + gridSize.y / 2 * blockEdgeSize);
			const long int z = (long int)pt.get(2) - (start.get(2) + gridSize.z / 2 * blockEdgeSize);

			float radius = sqrt((float) (x*x + y*y + z*z));

			bool is_active = radius < 64 && radius > 56;

			if (is_active == true)
			{match &= true;}
			else
			{match &= false;}
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
	BOOST_REQUIRE_EQUAL(cnt1,20520);
	BOOST_REQUIRE_EQUAL(cnt2,27152);
	BOOST_REQUIRE_EQUAL(cnt3,10423);
	BOOST_REQUIRE_EQUAL(cnt4,13649);
}

#if defined(OPENFPM_DATA_ENABLE_IO_MODULE) || defined(PERFORMANCE_TEST)

BOOST_AUTO_TEST_CASE(testSparseGridGpuOutput3DHeatStencil)
{
	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 4;
	typedef aggregate<float, float> AggregateT;

	size_t sz[3] = {512,512,512};
//        dim3 gridSize(128,128,128);
	dim3 gridSize(32,32,32);

	grid_smb<dim, blockEdgeSize> blockGeometry(sz);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64, long int> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	///// Insert sparse content, a set of 3 hollow spheres /////
	// Sphere 1
	grid_key_dx<3,int> start1({256,256,256});
	sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
	CUDA_LAUNCH_DIM3((insertSphere3D<0>),
					 gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
					 sparseGrid.toKernel(), start1, 64, 32, 1);
	cudaDeviceSynchronize();
	sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
	cudaDeviceSynchronize();

	// Sphere 2
	grid_key_dx<3,int> start2({192,192,192});
	sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
	CUDA_LAUNCH_DIM3((insertSphere3D<0>),
					 gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
					 sparseGrid.toKernel(), start2, 64, 44, 1);
	cudaDeviceSynchronize();
	sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
	cudaDeviceSynchronize();

	// Sphere 3
	grid_key_dx<3,int> start3({340,192,192});
	sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
	CUDA_LAUNCH_DIM3((insertSphere3D<0>),
					 gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
					 sparseGrid.toKernel(), start3, 20, 15, 1);
	cudaDeviceSynchronize();
	sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
	cudaDeviceSynchronize();
	///// /////

	sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
	sparseGrid.tagBoundaries(ctx);

	// Now apply some boundary conditions
	sparseGrid.template applyBoundaryStencils<BoundaryStencilSetXRescaled<dim,0,0>>(
			192, 384,
			0.0, 10.0);
	cudaDeviceSynchronize();

	// Now apply the laplacian operator
//        const unsigned int maxIter = 1000;
	const unsigned int maxIter = 100;
	for (unsigned int iter=0; iter<maxIter; ++iter)
	{
		for (int innerIter=0; innerIter<10; ++innerIter)
		{
			sparseGrid.applyStencils<HeatStencil<dim, 0, 1>>(STENCIL_MODE_INPLACE, 0.1);
			cudaDeviceSynchronize();
			sparseGrid.applyStencils<HeatStencil<dim, 1, 0>>(STENCIL_MODE_INPLACE, 0.1);
			cudaDeviceSynchronize();
		}
	}

	sparseGrid.deviceToHost<0,1>();
	sparseGrid.write("SparseGridGPU_output3DHeatStencil.vtk");
}

BOOST_AUTO_TEST_CASE(testSparseGridGpuOutput)
{
	constexpr unsigned int dim = 2;
	constexpr unsigned int blockEdgeSize = 8;
	typedef aggregate<float> AggregateT;

	size_t sz[2] = {1000000,1000000};
	dim3 gridSize(128,128);

	grid_smb<dim, blockEdgeSize> blockGeometry(sz);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64, long int> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	grid_key_dx<2,int> start({500000,500000});

	// Insert values on the grid
	sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
	CUDA_LAUNCH_DIM3((insertSphere<0>),gridSize, dim3(blockEdgeSize*blockEdgeSize,1),sparseGrid.toKernel(), start, 512, 256, 1);
	sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
	sparseGrid.tagBoundaries(ctx);

	sparseGrid.template deviceToHost<0>();

	sparseGrid.write("SparseGridGPU_output.vtk");
}

BOOST_AUTO_TEST_CASE(testSparseGridGpuOutput3D)
{
	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 4;
	typedef aggregate<float> AggregateT;

	size_t sz[3] = {512,512,512};
//        dim3 gridSize(128,128,128);
	dim3 gridSize(32,32,32);

	grid_smb<dim, blockEdgeSize> blockGeometry(sz);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, 64, long int> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	grid_key_dx<3,int> start({256,256,256});

	// Insert values on the grid
	sparseGrid.setGPUInsertBuffer(gridSize,dim3(1));
	CUDA_LAUNCH_DIM3((insertSphere3D<0>),
			gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
			sparseGrid.toKernel(), start, 64, 56, 1);
	sparseGrid.flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	cudaDeviceSynchronize();

	sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!
	sparseGrid.tagBoundaries(ctx);

	sparseGrid.template applyBoundaryStencils<BoundaryStencilSetX<dim,0,0>>();

	cudaDeviceSynchronize();

	sparseGrid.template deviceToHost<0>();

	sparseGrid.write("SparseGridGPU_output3D.vtk");

	bool test = compare("SparseGridGPU_output3D.vtk","test_data/SparseGridGPU_output3D_test.vtk");
	BOOST_REQUIRE_EQUAL(true,test);
}


#endif

BOOST_AUTO_TEST_SUITE_END()
