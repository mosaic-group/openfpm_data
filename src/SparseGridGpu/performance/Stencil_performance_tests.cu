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

extern char * test_dir;

// Property tree
struct report_sparse_grid_tests
{
	boost::property_tree::ptree graphs;
};

report_sparse_grid_tests report_sparsegrid_funcs;



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
            unsigned int offset,
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

        sparseGrid.loadGhostBlock<p_src>(dataBlockLoad,dataBlockId, enlargedBlock);

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
void testStencilHeat_perf(unsigned int i)
{
    auto testName = "In-place stencil";
    unsigned int gridEdgeSize = 1024;
    typedef HeatStencil<SparseGridZ::dims,0,1> StencilT;

    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").stencil");

    report_sparsegrid_funcs.graphs.put(base + ".dim",2);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",8);
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
    std::cout << "Grid: " << gridEdgeSize*SparseGridZ::blockEdgeSize_ << "x" << gridEdgeSize*SparseGridZ::blockEdgeSize_ << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\tStencil: " << mean_gf << " dev:" << deviation_gf << " s" << std::endl;
    std::cout << "Throughput: " << std::endl << "\t " << gElemS << " GElem/s " << std::endl << "\t " << gFlopsS << " GFlops/s" << std::endl;

    report_sparsegrid_funcs.graphs.put(base + ".Gflops.mean",mean_gf);
    report_sparsegrid_funcs.graphs.put(base +".Gflops.dev",deviation_gf);
    report_sparsegrid_funcs.graphs.put(base + ".time.mean",mean_tm);
    report_sparsegrid_funcs.graphs.put(base +".time.dev",deviation_tm);
}

BOOST_AUTO_TEST_SUITE(performance)

BOOST_AUTO_TEST_SUITE(SparseGridGpu_test)

BOOST_AUTO_TEST_CASE(testStencilHeat)
{
	constexpr unsigned int dim = 2;
	constexpr unsigned int blockEdgeSize = 8;

	typedef aggregate<float,float> AggregateT;
	constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

	report_sparsegrid_funcs.graphs.put("performance.SparseGridGpu(0).stencil.test.name","StencilN");

	testStencilHeat_perf<SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize>>(0);
}

BOOST_AUTO_TEST_CASE(testStencilHeatZ)
{
	constexpr unsigned int dim = 2;
	constexpr unsigned int blockEdgeSize = 8;

	typedef aggregate<float,float> AggregateT;
    constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;

    report_sparsegrid_funcs.graphs.put("performance.SparseGridGpu(1).stencil.test.name","StencilZ");

	testStencilHeat_perf<SparseGridGpu_z<dim, AggregateT, blockEdgeSize, chunkSize>>(1);
}

void testInsertStencil(unsigned int gridEdgeSize, unsigned int i)
{
	auto testName = "Insert stencil";
	constexpr unsigned int dim = 2;
	constexpr unsigned int blockEdgeSize = 8;
	constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
	typedef aggregate<float,float> AggregateT;
	typedef HeatStencil<dim,0,1> StencilT;

	unsigned int iterations = 10;

    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").insertStencil");

    report_sparsegrid_funcs.graphs.put(base + ".dim",2);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",8);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.x",gridEdgeSize*blockEdgeSize);
    report_sparsegrid_funcs.graphs.put(base + ".gridSize.y",gridEdgeSize*blockEdgeSize);

	openfpm::vector<double> measures;

	dim3 gridSize(gridEdgeSize, gridEdgeSize);
	dim3 blockSize(blockEdgeSize, blockEdgeSize);
	grid_smb<dim, blockEdgeSize> blockGeometry(gridSize);
	SparseGridGpu<dim, AggregateT, blockEdgeSize, chunkSize> sparseGrid(blockGeometry);
	mgpu::ofp_context_t ctx;
	sparseGrid.template setBackgroundValue<0>(0);

	// Initialize the grid
	sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
	CUDA_LAUNCH_DIM3((insertConstantValue<0>),gridSize, blockSize,sparseGrid.toKernel(), 0);
	sparseGrid.flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.setGPUInsertBuffer(gridSize, dim3(1));
	dim3 sourcePt(gridSize.x * blockEdgeSize / 2, gridSize.y * blockEdgeSize / 2, 0);
	insertOneValue<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), sourcePt, 100);
	sparseGrid.flush < sRight_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

	sparseGrid.findNeighbours(); // Pre-compute the neighbours pos for each block!

	unsigned long long numElements = gridEdgeSize*blockEdgeSize*gridEdgeSize*blockEdgeSize;

	for (unsigned int iter=0; iter<5; ++iter)
	{
		sparseGrid.applyStencils<StencilT>(STENCIL_MODE_INSERT, 0.1);
		sparseGrid.template flush<smax_<0>>(ctx, flush_type::FLUSH_ON_DEVICE);
	}

	for (unsigned int iter=0; iter<iterations; ++iter)
	{
		cudaDeviceSynchronize();

		timer ts;
		ts.start();

		sparseGrid.applyStencils<StencilT>(STENCIL_MODE_INSERT, 0.1);
		sparseGrid.template flush<smax_<0>>(ctx, flush_type::FLUSH_ON_DEVICE);

		cudaDeviceSynchronize();

		ts.stop();

		float gElemS = numElements / (1e6 * ts.getwct()) * StencilT::flops;

		measures.add(gElemS);
	}


	double mean = 0;
	double deviation = 0;
	standard_deviation(measures,mean,deviation);
	// All times above are in ms

    report_sparsegrid_funcs.graphs.put(base + ".GElems.mean",mean);
    report_sparsegrid_funcs.graphs.put(base +".GElems.dev",deviation);

	std::cout << "Test: " << testName << "\n";
	std::cout << "Grid: " << gridEdgeSize*blockEdgeSize << "x" << gridEdgeSize*blockEdgeSize << "\n";
	std::cout << "Iterations: " << iterations << "\n";
	std::cout << "Throughput:\n\t" << mean << " MElem/s dev: " << deviation << " GElem/s\n";
}

BOOST_AUTO_TEST_CASE(testStencilHeatInsert)
{
	testInsertStencil(256,0);
	testInsertStencil(512,1);
}

void testInsertSingle(unsigned int gridEdgeSize, unsigned int i)
{
	auto testName = "Insert (one chunk per element)";
	constexpr unsigned int dim = 2;
	constexpr unsigned int blockEdgeSize = 8;
	constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
	typedef aggregate<float> AggregateT;

	unsigned int iterations = 10;
	bool prePopulateGrid = true;

    std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").insertSingle");

    report_sparsegrid_funcs.graphs.put(base + ".dim",2);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",8);
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
		sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
		cudaDeviceSynchronize();
		///
	}

	for (unsigned int iter=0; iter<5; ++iter)
	{
		auto offset = 0;
		sparseGrid.setGPUInsertBuffer(gridSize, blockSize);
		insertValues2D<0> << < gridSize, blockSize >> > (sparseGrid.toKernel(), offset, offset);
		sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
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
		sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
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
	std::cout << "Grid: " << gridEdgeSize*blockEdgeSize << "x" << gridEdgeSize*blockEdgeSize << "\n";
	std::cout << "Iterations: " << iterations << "\n";
	std::cout << "Throughput:\n\t" << mean << "M/s" << "\n";
}

BOOST_AUTO_TEST_CASE(testInsert)
{
	testInsertSingle(64,0);
	testInsertSingle(128,1);
}

void test_insert_block(unsigned int gridEdgeSize, unsigned int i)
{
	auto testName = "Insert (one chunk per block)";
	constexpr unsigned int dim = 2;
	constexpr unsigned int blockEdgeSize = 8;
	constexpr unsigned int chunkSize = IntPow<blockEdgeSize,dim>::value;
	typedef aggregate<float> AggregateT;

	std::string base("performance.SparseGridGpu(" + std::to_string(i) + ").insert");

	report_sparsegrid_funcs.graphs.put(base + ".name","Block insert");
    report_sparsegrid_funcs.graphs.put(base + ".dim",2);
    report_sparsegrid_funcs.graphs.put(base + ".blockSize",8);
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
		sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);
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
		sparseGrid.flush < smax_ < 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

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
	std::cout << "Grid: " << gridEdgeSize*blockEdgeSize << "x" << gridEdgeSize*blockEdgeSize << "\n";
	std::cout << "Iterations: " << iterations << "\n";
	std::cout << "\tInsert: " << mean << " dev: " << deviation << " s" << std::endl;
	std::cout << "Throughput:\n\t" << mean << " MElem/s\n";
}

BOOST_AUTO_TEST_CASE(testInsertBlocked)
{
	test_insert_block(128,0);
	test_insert_block(256,1);
}

BOOST_AUTO_TEST_CASE(write_teport)
{
	report_sparsegrid_funcs.graphs.put("graphs.graph(0).type","line");
	report_sparsegrid_funcs.graphs.add("graphs.graph(0).title","SparseGridGPU stencil performance");
	report_sparsegrid_funcs.graphs.add("graphs.graph(0).x.title","test");
	report_sparsegrid_funcs.graphs.add("graphs.graph(0).y.title","Gflops");
	report_sparsegrid_funcs.graphs.add("graphs.graph(0).y.data(0).source","performance.SparseGridGpu(#).stencil.Gflops.mean");
	report_sparsegrid_funcs.graphs.add("graphs.graph(0).x.data(0).source","performance.SparseGridGpu(#).stencil.test.name");
	report_sparsegrid_funcs.graphs.add("graphs.graph(0).y.data(0).title","line");

	report_sparsegrid_funcs.graphs.put("graphs.graph(1).type","line");
	report_sparsegrid_funcs.graphs.add("graphs.graph(1).title","SparseGridGPU insert blocked performance");
	report_sparsegrid_funcs.graphs.add("graphs.graph(1).x.title","size");
	report_sparsegrid_funcs.graphs.add("graphs.graph(1).y.title","Milion inserts");
	report_sparsegrid_funcs.graphs.add("graphs.graph(1).y.data(0).source","performance.SparseGridGpu(#).insert.Minsert.mean");
	report_sparsegrid_funcs.graphs.add("graphs.graph(1).x.data(0).source","performance.SparseGridGpu(#).insert.gridSize.x");
	report_sparsegrid_funcs.graphs.add("graphs.graph(1).y.data(0).title","line");

	report_sparsegrid_funcs.graphs.put("graphs.graph(2).type","line");
	report_sparsegrid_funcs.graphs.add("graphs.graph(2).title","SparseGridGPU insert single performance");
	report_sparsegrid_funcs.graphs.add("graphs.graph(2).x.title","size");
	report_sparsegrid_funcs.graphs.add("graphs.graph(2).y.title","Milion inserts");
	report_sparsegrid_funcs.graphs.add("graphs.graph(2).y.data(0).source","performance.SparseGridGpu(#).insertSingle.Minsert.mean");
	report_sparsegrid_funcs.graphs.add("graphs.graph(2).x.data(0).source","performance.SparseGridGpu(#).insertSingle.gridSize.x");
	report_sparsegrid_funcs.graphs.add("graphs.graph(2).y.data(0).title","line");

	report_sparsegrid_funcs.graphs.put("graphs.graph(3).type","line");
	report_sparsegrid_funcs.graphs.add("graphs.graph(3).title","SparseGridGPU insert single performance");
	report_sparsegrid_funcs.graphs.add("graphs.graph(3).x.title","size");
	report_sparsegrid_funcs.graphs.add("graphs.graph(3).y.title","Milion inserts");
	report_sparsegrid_funcs.graphs.add("graphs.graph(3).y.data(0).source","performance.SparseGridGpu(#).insertStencil.GElems.mean");
	report_sparsegrid_funcs.graphs.add("graphs.graph(3).x.data(0).source","performance.SparseGridGpu(#).insertStencil.gridSize.x");
	report_sparsegrid_funcs.graphs.add("graphs.graph(3).y.data(0).title","line");

	boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
	boost::property_tree::write_xml("SparseGridGpu_performance.xml", report_sparsegrid_funcs.graphs,std::locale(),settings);

	std::string file_xml_ref(test_dir);
	file_xml_ref += std::string("/openfpm_pdata/SparseGridGpu_performance_ref.xml");

	GoogleChart cg;

	StandardXMLPerformanceGraph("SparseGridGpu_performance.xml",file_xml_ref,cg);

	addUpdtateTime(cg,1);
	cg.write("SparseGridGpu_performance.html");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
