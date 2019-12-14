//
// Created by tommaso on 15/8/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_TESTKERNELS_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_TESTKERNELS_CUH

/////////////////// BOUNDARY STENCILS ////////////////////////

template<unsigned int dim, unsigned int p_src, unsigned int p_dst>
struct BoundaryStencilSetX
{
    // This is an example of a boundary stencil setting the value to the same value as the x coordinate

    typedef NNStar<dim> stencil_type;

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
            bool applyStencilHere)
    {
        if (applyStencilHere)
        {
            dataBlockStore.template get<p_dst>()[offset] = pointCoord.get(0);
        }
    }
};

template<unsigned int dim, unsigned int p_src, unsigned int p_dst, typename ScalarT = float>
struct BoundaryStencilSetXRescaled
{
    // This is an example of a boundary stencil setting the value to the same value as the x coordinate

    typedef NNStar<dim> stencil_type;

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
            ScalarT minX, ScalarT maxX, ScalarT minValue, ScalarT maxValue)
    {
        if (applyStencilHere)
        {
            const ScalarT x = pointCoord.get(0);
            auto value = maxValue * (x - minX) / (maxX - minX - 1);
            if (x < minX)
            {
                value = minValue;
            }
            else if (x > maxX)
            {
                value = maxValue;
            }
            dataBlockStore.template get<p_dst>()[offset] = value;
        }
    }
};

/////////////////// KERNELS ////////////////////////

template<unsigned int p, typename SparseGridType, typename ValueType>
__global__ void insertSphere(SparseGridType sparseGrid, grid_key_dx<2,int> start, float r1, float r2, ValueType value)
{
    constexpr unsigned int pMask = SparseGridType::pMask;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p> BlockT;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, pMask> MaskBlockT;

    grid_key_dx<2,int> blk({
        blockIdx.x + start.get(0) / sparseGrid.getBlockEdgeSize(),
        blockIdx.y + start.get(1) / sparseGrid.getBlockEdgeSize()
    });

    unsigned int offset = threadIdx.x;

    __shared__ bool is_block_empty;

    if (threadIdx.x == 0 && threadIdx.y == 0)
    {is_block_empty = true;}

    sparseGrid.init();

    auto blockId = sparseGrid.getBlockLinId(blk);

    grid_key_dx<2,int> keyg;
    keyg = sparseGrid.getGlobalCoord(blk,offset);

    float radius = sqrt( (float)
            (keyg.get(0) - (start.get(0) + gridDim.x/2*SparseGridType::blockEdgeSize_))
            * (keyg.get(0) - (start.get(0) + gridDim.x/2*SparseGridType::blockEdgeSize_))
            + (keyg.get(1) - (start.get(1) + gridDim.y/2*SparseGridType::blockEdgeSize_))
            * (keyg.get(1) - (start.get(1) + gridDim.y/2*SparseGridType::blockEdgeSize_)) );

    bool is_active = radius < r1 && radius > r2;

    if (is_active == true)
    {
        is_block_empty = false;
    }

    __syncthreads();

    if (is_block_empty == false)
    {
        auto ec = sparseGrid.insertBlock(blockId);

        if ( is_active == true)
        {
            ec.template get<p>()[offset] = value;
            BlockMapGpu_ker<>::setExist(ec.template get<pMask>()[offset]);
        }
    }

    __syncthreads();

    sparseGrid.flush_block_insert();
}

template<unsigned int p, typename SparseGridType, typename ValueType>
__global__ void insertSphere3D(SparseGridType sparseGrid, grid_key_dx<3,int> start, float r1, float r2, ValueType value)
{
    constexpr unsigned int pMask = SparseGridType::pMask;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p> BlockT;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, pMask> MaskBlockT;

    grid_key_dx<3,int> blk({
                                   blockIdx.x + start.get(0) / sparseGrid.getBlockEdgeSize(),
                                   blockIdx.y + start.get(1) / sparseGrid.getBlockEdgeSize(),
                                   blockIdx.z + start.get(2) / sparseGrid.getBlockEdgeSize()});

    unsigned int offset = threadIdx.x;

    __shared__ bool is_block_empty;

    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {is_block_empty = true;}

    sparseGrid.init();

    auto blockId = sparseGrid.getBlockLinId(blk);

    grid_key_dx<3,int> keyg;
    keyg = sparseGrid.getGlobalCoord(blk,offset);

    const long int x = (long int)keyg.get(0) - (start.get(0) + gridDim.x / 2 * SparseGridType::blockEdgeSize_);
    const long int y = (long int)keyg.get(1) - (start.get(1) + gridDim.y / 2 * SparseGridType::blockEdgeSize_);
    const long int z = (long int)keyg.get(2) - (start.get(2) + gridDim.z / 2 * SparseGridType::blockEdgeSize_);

    float radius = sqrt((float) (x*x + y*y + z*z));

    bool is_active = radius < r1 && radius >= r2;

    if (is_active == true)
    {
        is_block_empty = false;
    }

    __syncthreads();

    if (is_block_empty == false)
    {
        auto ec = sparseGrid.insertBlock(blockId);

        if ( is_active == true)
        {
            ec.template get<p>()[offset] = value;
            BlockMapGpu_ker<>::setExist(ec.template get<pMask>()[offset]);
        }
    }

    __syncthreads();

    sparseGrid.flush_block_insert();
}

template<unsigned int p, typename SparseGridType, typename ValueType>
__global__ void insertSphere3D_radius(SparseGridType sparseGrid, grid_key_dx<3,int> start, float r1, float r2, ValueType value)
{
    constexpr unsigned int pMask = SparseGridType::pMask;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p> BlockT;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, pMask> MaskBlockT;

    grid_key_dx<3,int> blk({
                                   blockIdx.x + start.get(0) / sparseGrid.getBlockEdgeSize(),
                                   blockIdx.y + start.get(1) / sparseGrid.getBlockEdgeSize(),
                                   blockIdx.z + start.get(2) / sparseGrid.getBlockEdgeSize()});

    unsigned int offset = threadIdx.x;

    __shared__ bool is_block_empty;

    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {is_block_empty = true;}

    sparseGrid.init();

    auto blockId = sparseGrid.getBlockLinId(blk);

    grid_key_dx<3,int> keyg;
    keyg = sparseGrid.getGlobalCoord(blk,offset);

    const long int x = (long int)keyg.get(0) - (start.get(0) + gridDim.x / 2 * SparseGridType::blockEdgeSize_);
    const long int y = (long int)keyg.get(1) - (start.get(1) + gridDim.y / 2 * SparseGridType::blockEdgeSize_);
    const long int z = (long int)keyg.get(2) - (start.get(2) + gridDim.z / 2 * SparseGridType::blockEdgeSize_);

    float radius = sqrt((float) (x*x + y*y + z*z));

    bool is_active = radius < r1 && radius > r2;

    if (is_active == true)
    {
        is_block_empty = false;
    }

    __syncthreads();

    if (is_block_empty == false)
    {
        auto ec = sparseGrid.insertBlock(blockId);

        if ( is_active == true)
        {
            ec.template get<p>()[offset] = x+y+z;
            BlockMapGpu_ker<>::setExist(ec.template get<pMask>()[offset]);
        }
    }

    __syncthreads();

    sparseGrid.flush_block_insert();
}

template<unsigned int p, typename SparseGridType, typename ValueType>
__global__ void insertSphere3D_radiusV(SparseGridType sparseGrid, grid_key_dx<3,int> start, float r1, float r2, ValueType value)
{
    constexpr unsigned int pMask = SparseGridType::pMask;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p> BlockT;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, pMask> MaskBlockT;

    grid_key_dx<3,int> blk({
                                   blockIdx.x + start.get(0) / sparseGrid.getBlockEdgeSize(),
                                   blockIdx.y + start.get(1) / sparseGrid.getBlockEdgeSize(),
                                   blockIdx.z + start.get(2) / sparseGrid.getBlockEdgeSize()});

    unsigned int offset = threadIdx.x;

    __shared__ bool is_block_empty;

    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {is_block_empty = true;}

    sparseGrid.init();

    auto blockId = sparseGrid.getBlockLinId(blk);

    grid_key_dx<3,int> keyg;
    keyg = sparseGrid.getGlobalCoord(blk,offset);

    const long int x = (long int)keyg.get(0) - (start.get(0) + gridDim.x / 2 * SparseGridType::blockEdgeSize_);
    const long int y = (long int)keyg.get(1) - (start.get(1) + gridDim.y / 2 * SparseGridType::blockEdgeSize_);
    const long int z = (long int)keyg.get(2) - (start.get(2) + gridDim.z / 2 * SparseGridType::blockEdgeSize_);

    float radius = sqrt((float) (x*x + y*y + z*z));

    bool is_active = radius < r1 && radius > r2;

    if (is_active == true)
    {
        is_block_empty = false;
    }

    __syncthreads();

    if (is_block_empty == false)
    {
        auto ec = sparseGrid.insertBlock(blockId);

        if ( is_active == true)
        {
            ec.template get<p>()[offset] = x+y+z;
            ec.template get<p+1>()[0][offset] = x;
            ec.template get<p+1>()[1][offset] = y;
            ec.template get<p+1>()[2][offset] = z;
            BlockMapGpu_ker<>::setExist(ec.template get<pMask>()[offset]);
        }
    }

    __syncthreads();

    sparseGrid.flush_block_insert();
}

template<unsigned int p, typename SparseGridType, typename ValueType>
__global__ void removeSphere3D_even_radiusV(SparseGridType sparseGrid, grid_key_dx<3,int> start, float r1, float r2, ValueType value)
{
    constexpr unsigned int pMask = SparseGridType::pMask;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p> BlockT;
    typedef BlockTypeOf<typename SparseGridType::AggregateType, pMask> MaskBlockT;

    grid_key_dx<3,int> blk({
                                   blockIdx.x + start.get(0) / sparseGrid.getBlockEdgeSize(),
                                   blockIdx.y + start.get(1) / sparseGrid.getBlockEdgeSize(),
                                   blockIdx.z + start.get(2) / sparseGrid.getBlockEdgeSize()});

    unsigned int offset = threadIdx.x;

    auto blockId = sparseGrid.getBlockLinId(blk);

    grid_key_dx<3,int> keyg;
    keyg = sparseGrid.getGlobalCoord(blk,offset);

    const long int x = (long int)keyg.get(0) - (start.get(0) + gridDim.x / 2 * SparseGridType::blockEdgeSize_);
    const long int y = (long int)keyg.get(1) - (start.get(1) + gridDim.y / 2 * SparseGridType::blockEdgeSize_);
    const long int z = (long int)keyg.get(2) - (start.get(2) + gridDim.z / 2 * SparseGridType::blockEdgeSize_);

    float radius = sqrt((float) (x*x + y*y + z*z));

    bool is_active = radius < r1 && radius > r2 && (keyg.get(0) + keyg.get(1) + keyg.get(2)) % 2 == 0;

    if (is_active == true)
    {
    	sparseGrid.remove(keyg);
    }
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_TESTKERNELS_CUH
