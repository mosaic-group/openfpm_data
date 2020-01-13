/*
 * SparseGridGpu_util_test.cuh
 *
 *  Created on: Sep 9, 2019
 *      Author: i-bird
 */

#ifndef SPARSEGRIDGPU_UTIL_TEST_CUH_
#define SPARSEGRIDGPU_UTIL_TEST_CUH_

#include "SparseGridGpu/tests/utils/SparseGridGpu_testKernels.cuh"

template<unsigned int p, typename SparseGridType>
__global__ void getValues2D(SparseGridType sparseGrid, const int offsetX=0, const int offsetY=0)
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

    auto value = sparseGrid.template get<p>(coord);
    value++;

    // Compiler avoid warning
    x++;
    y++;
}

template<unsigned int p, typename SparseGridType>
__global__ void getValuesNeighbourhood2D(SparseGridType sparseGrid, const int offsetX=0, const int offsetY=0)
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
    --x; --y;
    grid_key_dx<SparseGridType::d, int> coord({x, y});
    for (int i=0; i < 9; ++i)
    {
        auto value = sparseGrid.template get<p>(coord);
        coord.set_d(0, x + i%3);
        coord.set_d(1, y + i/3);
        value++;
    }
}

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

    int posX = blockIdx.x * blockDim.x + threadIdx.x + sOffsetX;
    int posY = blockIdx.y * blockDim.y + threadIdx.y + sOffsetY;
    const unsigned int offsetX = posX % blockEdgeSize;
    const unsigned int offsetY = posY % blockEdgeSize;

    const unsigned int offset = offsetY * blockEdgeSize + offsetX;

	grid_key_dx<SparseGridType::d, int> blockCoord({posX / blockEdgeSize, posY / blockEdgeSize});
	auto encap = sparseGrid.insertBlock(sparseGrid.getBlockLinId(blockCoord));

    encap.template get<p>()[offset] = posX*posX * posY*posY;
    BlockMapGpu_ker<>::setExist(encap.template get<pMask>()[offset]);

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

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    grid_key_dx<SparseGridType::d, size_t> coord({x, y, z});

    auto pos = sparseGrid.getLinId(coord);
    unsigned int dataBlockId = pos / BlockT::size;
    unsigned int offset = pos % BlockT::size;

    auto encap = sparseGrid.template insertBlock<chunksPerBlock>(dataBlockId,BlockT::size);

    encap.template get<p>()[offset] = value;
    BlockMapGpu_ker<>::setExist(encap.template get<pMask>()[offset]);

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
	typedef NNStar<dim> stencil_type;

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

        sparseGrid.template loadGhostBlock<p_src>(dataBlockLoad, dataBlockIdPos, enlargedBlock);

        __syncthreads();

        decltype(sparseGrid.getLinIdInEnlargedBlock(0)) linId = 0;
        ScalarT res = 0;

        if (isActive)
        {
            const auto coord = sparseGrid.getCoordInEnlargedBlock(offset);
            // const auto linId = sparseGrid.getLinIdInEnlargedBlock(offset);
            linId = sparseGrid.getLinIdInEnlargedBlock(offset);
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
            // enlargedBlock[linId] = cur + dt * laplacian;
            res = cur + dt * laplacian;
        }
        __syncthreads();
        if (isActive)
        {
            enlargedBlock[linId] = res;
        }
        __syncthreads();
        sparseGrid.template storeBlock<p_dst>(dataBlockStore, enlargedBlock);
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
        sparseGrid.template flush <smax_<0>> (ctx, flush_type::FLUSH_ON_DEVICE);
    }
};

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

    static constexpr unsigned int flops = 2*27;

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

        sparseGrid.loadGhostBlock<p_src>(dataBlockLoad,dataBlockIdPos,enlargedBlock);

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

    static constexpr unsigned int flops = 2*27;

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

                    	auto pos = sparseGrid.template getNNPoint<stencil_type>(dataBlockIdPos, offset, key);

                    	tot += sparseGrid.template get<p_src>(pos) * cc.coeff[i][j][k];
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
void testConv3x3x3_perf(std::string testName)
{
	constexpr unsigned int dim = 3;
	constexpr unsigned int blockEdgeSize = 8;
	typedef aggregate<float,float> AggregateT;

	unsigned int iterations = 100;

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
	            gridSize, dim3(blockEdgeSize*blockEdgeSize*blockEdgeSize,1,1),
	            sparseGrid.toKernel(), start,128, 56, 1);

	sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    sparseGrid.template deviceToHost<0>(); // NECESSARY as count takes place on Host!
    auto existingElements = sparseGrid.countExistingElements();
    auto boundaryElements = sparseGrid.countBoundaryElements();
    unsigned long long numElements = existingElements - boundaryElements;

	sparseGrid.template findNeighbours<NNFull<3>>(); // Pre-compute the neighbours pos for each block!

	sparseGrid.template setNNType<NNFull<dim>>();
	sparseGrid.template tagBoundaries<NNFull<3>>(ctx);

    openfpm::vector<double> measures_gf;
    openfpm::vector<double> measures_tm;

	for (unsigned int iter=0; iter<iterations; ++iter)
	{
		hipDeviceSynchronize();
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

		timer ts;
		ts.start();

		sparseGrid.template applyStencils<Conv3x3x3<dim,0,1>>(STENCIL_MODE_INPLACE,cc);

		hipDeviceSynchronize();
		ts.stop();

		measures_tm.add(ts.getwct());

	    float gElemS = numElements / (1e9 * ts.getwct());
	    float gFlopsS = gElemS * Conv3x3x3<dim,0,1>::flops;

		measures_gf.add(gFlopsS);
	}

    double mean_tm = 0;
    double deviation_tm = 0;
    standard_deviation(measures_tm,mean_tm,deviation_tm);

    double mean_gf = 0;
    double deviation_gf = 0;
    standard_deviation(measures_gf,mean_gf,deviation_gf);

    std::cout << "Test: " << testName << std::endl;
    std::cout << "Block: " << SparseGridZ::blockEdgeSize_
              << "x" << SparseGridZ::blockEdgeSize_
              << "x" << SparseGridZ::blockEdgeSize_
              << std::endl;
    std::cout << "Grid: " << sz[0] << "x" << sz[1] << "x" << sz[2] << std::endl;

    double dataOccupancyMean, dataOccupancyDev;
    sparseGrid.measureBlockOccupancy(dataOccupancyMean, dataOccupancyDev);
    std::cout << "Data Occupancy: " << dataOccupancyMean << " dev:" << dataOccupancyDev << std::endl;

    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\tConvolution3x3x3: " << mean_tm << " dev:" << deviation_tm << " s" << std::endl;
    std::cout << "Throughput: " << std::endl << "\t " << mean_gf << " GFlops/s   dev: " <<  deviation_gf << "   GFlops/s " << std::endl;
}

template<typename SparseGridZ>
static void testConv3x3x3_no_shared_perf(std::string testName)
{
	unsigned int iterations = 100;

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
	            sparseGrid.toKernel(), start,128, 56, 1);

	sparseGrid.template flush < smax_< 0 >> (ctx, flush_type::FLUSH_ON_DEVICE);

    sparseGrid.template deviceToHost<0>(); // NECESSARY as count takes place on Host!
    auto existingElements = sparseGrid.countExistingElements();
    auto boundaryElements = sparseGrid.countBoundaryElements();
    unsigned long long numElements = existingElements - boundaryElements;

	sparseGrid.template findNeighbours<NNFull<3>>(); // Pre-compute the neighbours pos for each block!

	sparseGrid.template setNNType<NNFull<SparseGridZ::dims>>();
	sparseGrid.template tagBoundaries<NNFull<3>>(ctx,tag_boundaries::CALCULATE_EXISTING_POINTS);


    openfpm::vector<double> measures_gf;
    openfpm::vector<double> measures_tm;

	for (unsigned int iter=0; iter<iterations; ++iter)
	{
		hipDeviceSynchronize();
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

		timer ts;
		ts.start();

		sparseGrid.template applyStencils<Conv3x3x3_noshared<SparseGridZ::dims,0,1>>(STENCIL_MODE_INPLACE_NO_SHARED,cc);

		hipDeviceSynchronize();
		ts.stop();

		measures_tm.add(ts.getwct());

	    float gElemS = numElements / (1e9 * ts.getwct());
	    float gFlopsS = gElemS * Conv3x3x3_noshared<SparseGridZ::dims,0,1>::flops;

		measures_gf.add(gFlopsS);
	}

    double mean_tm = 0;
    double deviation_tm = 0;
    standard_deviation(measures_tm,mean_tm,deviation_tm);

    double mean_gf = 0;
    double deviation_gf = 0;
    standard_deviation(measures_gf,mean_gf,deviation_gf);

    std::cout << "Test: " << testName << std::endl;
    std::cout << "Block: " << SparseGridZ::blockEdgeSize_
              << "x" << SparseGridZ::blockEdgeSize_
              << "x" << SparseGridZ::blockEdgeSize_
              << std::endl;
    std::cout << "Grid: " << sz[0] << "x" << sz[1] << "x" << sz[2] << std::endl;

    double dataOccupancyMean, dataOccupancyDev;
    sparseGrid.measureBlockOccupancy(dataOccupancyMean, dataOccupancyDev);
    std::cout << "Data Occupancy: " << dataOccupancyMean << " dev:" << dataOccupancyDev << std::endl;

    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\tConvolution3x3x3: " << mean_tm << " dev:" << deviation_tm << " s" << std::endl;
    std::cout << "Throughput: " << std::endl << "\t " << mean_gf << " GFlops/s   dev: " <<  deviation_gf << "   GFlops/s " << std::endl;
}

template<unsigned int dim, unsigned int p_src, unsigned int p_dst>
struct HeatStencilGet
{
    typedef NNStar<dim> stencil_type;

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

        ScalarT res = 0;
        auto coord = pointCoord;

        if (isActive)
        {
            // ScalarT cur = sparseGrid.template get<p_src>(pointCoord);
            ScalarT cur = dataBlockLoad.template get<p_dst>()[offset];
            ScalarT laplacian = -2.0 * dim * cur; // The central part of the stencil

            for (int d = 0; d < dim; ++d)
            {
                auto locC = coord.get(d);
                coord.set_d(d, locC+1);
                auto nPlus = sparseGrid.template get<p_src>(coord);
                coord.set_d(d, locC-1);
                auto nMinus = sparseGrid.template get<p_src>(coord);
                laplacian += nMinus + nPlus;
                coord.set_d(d, locC);
            }
            res = cur + dt * laplacian;
        }
        __syncthreads();
        if (isActive)
        {
            dataBlockStore.template get<p_dst>()[offset] = res;
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

template<unsigned int dim, unsigned int p_src, unsigned int p_dst>
struct SkeletonStencil
{
    typedef NNStar<dim> stencil_type;

    // This is an example of a laplacian smoothing stencil to apply using the apply stencil facility of SparseGridGpu

    static constexpr unsigned int flops = 1;

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

        decltype(sparseGrid.getLinIdInEnlargedBlock(0)) linId = 0;
        ScalarT res = 0;

        if (isActive)
        {
            linId = sparseGrid.getLinIdInEnlargedBlock(offset);
            ScalarT cur = enlargedBlock[linId];

            res = cur + cur;
        }
        __syncthreads();
        if (isActive)
        {
            enlargedBlock[linId] = res;
        }
        __syncthreads();
        sparseGrid.storeBlock<p_dst>(dataBlockStore, enlargedBlock);
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




#endif /* SPARSEGRIDGPU_UTIL_TEST_CUH_ */
