//
// Created by tommaso on 11/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_KER_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_KER_CUH

#include "BlockMapGpu_ker.cuh"

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
class SparseGridGpu_ker : public BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>
{
private:
    GridSmT grid;

protected:
    const static unsigned char PADDING_BIT = 1;
    static constexpr unsigned int blockSize = BlockTypeOf<AggregateBlockT, 0>::size;
    // Good reason on why not to use constant memory for ghostLayer mapping is here https://stackoverflow.com/a/18021374
    unsigned int stencilSupportRadius;
    unsigned int *ghostLayerToThreadsMapping;

public:
    static constexpr unsigned int d = dim;

public:
    SparseGridGpu_ker(const openfpm::vector_sparse_gpu_ker<AggregateBlockT, indexT, layout_base> &blockMap,
                      GridSmT grid,
                      unsigned int stencilSupportRadius,
                      unsigned int *mapping,
                      unsigned int ghostLayerSize)
            : BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>(blockMap),
              grid(grid),
              stencilSupportRadius(stencilSupportRadius)
    {
        cudaMalloc(&ghostLayerToThreadsMapping, ghostLayerSize * sizeof(unsigned int));
        cudaMemcpy(ghostLayerToThreadsMapping, mapping, ghostLayerSize * sizeof(unsigned int), cudaMemcpyHostToDevice);
    }

    // Geometry
    template <typename CoordT>
    inline __device__ size_t getLinId(CoordT coord)
    {
        return grid.LinId(coord);
    }

    // Data management methods

    template<unsigned int p, typename CoordT>
    inline __device__ auto
    get(CoordT coord) const -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template get<p>(
            0));

    template<unsigned int p, typename CoordT>
    inline __device__ auto
    insert(CoordT coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template insert<p>(
            0));

    template<typename CoordT>
    inline __device__ unsigned int getBlockId(CoordT coord);

    template<typename CoordT>
    inline __device__ auto
    getBlock(CoordT coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlock(0));

    template<typename CoordT>
    inline __device__ auto
    insertBlock(CoordT coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock(0));

    // Load & Store aux functions for user kernels. To be used for loading to or writing from shared memory.

    /**
     * Load a data block into the inner part of a shared memory region.
     * The given shared memory region should be shaped as a dim-dimensional array and sized so that it
     * can contain the block plus the ghost layer around it.
     *
     * @tparam p The property to retrieve from global memory.
     * @tparam CoordT The coordinate type.
     * @tparam Shared The type of the shared memory region.
     * @param coord The coordinate of the block.
     * @param sharedRegionPtr The pointer to the shared memory region.
     */
//    template<unsigned int p, typename CoordT>
//    inline __device__ void
//    loadBlock(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion, unsigned int stencilSupportRadius);
    template<unsigned int p, typename CoordT>
    inline __device__ void
    loadBlock(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion);

    /**
     * Load a data block into the inner part of a shared memory region.
     * The given shared memory region should be shaped as a dim-dimensional array and sized so that it
     * can contain the block plus the ghost layer around it.
     *
     * @tparam CoordT The coordinate type.
     * @tparam props The set of properties to retrieve from global memory.
     * @param coord The coordinate of the block.
     * @param sharedRegionPtr The array of pointers to the shared memory regions, one for each property.
     */
//    template<unsigned int ... props, typename CoordT>
//    inline __device__ void loadBlock(CoordT coord, void *sharedRegionPtr[sizeof...(props)],
//                                     unsigned int stencilSupportRadius[sizeof...(props)]);
    template<unsigned int ... props, typename CoordT>
    inline __device__ void loadBlock(CoordT coord, void *sharedRegionPtr[sizeof...(props)]);

    /**
     * Load the ghost layer of a data block into the boundary part of a shared memory region.
     * The given shared memory region should be shaped as a dim-dimensional array and sized so that it
     * can contain the block plus the ghost layer around it.
     *
     * @tparam p The property to retrieve from global memory.
     * @tparam CoordT The coordinate type.
     * @tparam Shared The type of the shared memory region.
     * @param coord The coordinate of the block.
     * @param sharedRegionPtr The pointer to the shared memory region.
     */
    template<unsigned int p, typename CoordT>
    inline __device__ void
    loadGhost(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion);

    /**
     * Load the ghost layer of a data block into the boundary part of a shared memory region.
     * The given shared memory region should be shaped as a dim-dimensional array and sized so that it
     * can contain the block plus the ghost layer around it.
     *
     * @tparam CoordT The coordinate type.
     * @tparam Shared The type of the shared memory region.
     * @tparam props The set of properties to retrieve from global memory.
     * @param coord The coordinate of the block.
     * @param sharedRegionPtr The pointer to the shared memory region.
     */
    template<unsigned int ... props, typename CoordT>
    inline __device__ void loadGhost(CoordT coord, void *sharedRegionPtr[sizeof...(props)]);

    /**
     * Read a data block from the inner part of a shared memory region and store it in global memory.
     * The given shared memory region should be shaped as a dim-dimensional array and sized so that it
     * can contain the block plus the ghost layer around it.
     *
     * @tparam p The property to store into global memory.
     * @tparam CoordT The coordinate type.
     * @tparam Shared The type of the shared memory region.
     * @param coord The coordinate of the block.
     * @param sharedRegionPtr The pointer to the shared memory region.
     */
//    template<unsigned int p, typename CoordT>
//    inline __device__ void
//    storeBlock(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion, unsigned int stencilSupportRadius);
    template<unsigned int p, typename CoordT>
    inline __device__ void
    storeBlock(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion);

    /**
     * Read a data block from the inner part of a shared memory region and store it in global memory.
     * The given shared memory region should be shaped as a dim-dimensional array and sized so that it
     * can contain the block plus the ghost layer around it.
     *
     * @tparam CoordT The coordinate type.
     * @tparam Shared The type of the shared memory region.
     * @tparam props The set of properties to retrieve from global memory.
     * @param coord The coordinate of the block.
     * @param sharedRegionPtr The array of pointers to the shared memory regions, one for each property.
     */
//    template<unsigned int ... props, typename CoordT>
//    inline __device__ void storeBlock(CoordT coord, void *sharedRegionPtr[sizeof...(props)],
//                                      unsigned int stencilSupportRadius[sizeof...(props)]);
    template<unsigned int ... props, typename CoordT>
    inline __device__ void storeBlock(CoordT coord, void *sharedRegionPtr[sizeof...(props)]);

private:
    template<unsigned int edgeSize>
    void linToCoordWithOffset(unsigned int linId, unsigned int offset, unsigned int coord[dim])
    {
        for (unsigned int d = 0; d < dim; ++d)
        {
            coord[d] = linId % edgeSize;
            coord[d] += offset;
            linId /= edgeSize;
        }
    }

    template<unsigned int edgeSize>
    unsigned int coordToLin(unsigned int coord[dim])
    {
        unsigned int linId = coord[dim - 1];
        for (unsigned int d = dim - 2; d >= 0; --d)
        {
            linId *= edgeSize;
            linId += coord[d];
        }
        return linId;
    }

    template<unsigned int p, typename AggrWrapperT>
    inline __device__ void
//    __loadBlock(const AggrWrapperT &block, void *sharedRegionPtr[1],
//                unsigned int stencilSupportRadius[1])
    __loadBlock(const AggrWrapperT &block, void *sharedRegionPtr[1])
    {
        const unsigned int pos = threadIdx.x;
        //todo: Improve this version to allow multiple chunks per block!
        if (pos < blockSize)
        {
            // Convert pos into a linear id accounting for the ghost offsets
            unsigned int coord[dim];
            linToCoordWithOffset<blockEdgeSize>(pos, stencilSupportRadius, coord);
            const unsigned int linId = coordToLin<blockEdgeSize>(coord);
            // Actually load the data into the shared region
            *(
                    static_cast<ScalarTypeOf<AggregateBlockT, p> *>(sharedRegionPtr[0])
                    + linId
            ) = block.template get<p>(pos);
        }
    }

    template<unsigned int p, unsigned int ... props, typename AggrWrapperT>
    inline __device__ void
//    __loadBlock(const AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)],
//                unsigned int stencilSupportRadius[sizeof...(props)])
    __loadBlock(const AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)])
    {
//        __loadBlock<p>(block, sharedRegionPtr, stencilSupportRadius);
//        __loadBlock<props ...>(block, sharedRegionPtr + 1, stencilSupportRadius + 1);
        __loadBlock<p>(block, sharedRegionPtr);
        __loadBlock<props ...>(block, sharedRegionPtr + 1);
    }

    template<unsigned int p, typename AggrWrapperT>
    inline __device__ void
    __loadGhost(const AggrWrapperT &block, void *sharedRegionPtr[1])
    {
        const unsigned int pos = threadIdx.x;
        //todo: Improve this version to allow multiple chunks per block!
        if (pos < blockSize)
        {
            // Convert pos into a linear id accounting for the inner domain offsets
            const unsigned int linId = ghostLayerToThreadsMapping[pos];
            // Actually load the data into the shared region
            *(
                    static_cast<ScalarTypeOf<AggregateBlockT, p> *>(sharedRegionPtr[0])
                    + linId
            ) = block.template get<p>(pos);
        }
    }

    template<unsigned int p, unsigned int ... props, typename AggrWrapperT>
    inline __device__ void
    __loadGhost(const AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)])
    {
        __loadGhost<p>(block, sharedRegionPtr);
        __loadGhost<props ...>(block, sharedRegionPtr + 1);
    }

    template<unsigned int p, typename AggrWrapperT>
    inline __device__ void
//    __storeBlock(AggrWrapperT &block, void *sharedRegionPtr[1],
//                 unsigned int stencilSupportRadius[1])
    __storeBlock(AggrWrapperT &block, void *sharedRegionPtr[1])
    {
        const unsigned int pos = threadIdx.x;
        //todo: Improve this version to allow multiple chunks per block!
        if (pos < blockSize)
        {
            // Convert pos into a linear id accounting for the ghost offsets
            unsigned int coord[dim];
            linToCoordWithOffset<blockEdgeSize>(pos, stencilSupportRadius, coord);
            const unsigned int linId = coordToLin<blockEdgeSize>(coord);
            // Actually store the data from the shared region
            block.template get<p>(pos)
                    = *(
                    static_cast<ScalarTypeOf<AggregateBlockT, p> *>(sharedRegionPtr[0])
                    + linId);
        }
    }

    template<unsigned int p, unsigned int ... props, typename AggrWrapperT>
    inline __device__ void
//    __storeBlock(AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)],
//                 unsigned int stencilSupportRadius[sizeof...(props)])
    __storeBlock(AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)])
    {
//        __storeBlock<p>(block, sharedRegionPtr, stencilSupportRadius);
//        __storeBlock<props ...>(block, sharedRegionPtr + 1, stencilSupportRadius + 1);
        __storeBlock<p>(block, sharedRegionPtr);
        __storeBlock<props ...>(block, sharedRegionPtr + 1);
    }
};

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int p, typename CoordT>
inline __device__ auto SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::get(CoordT coord) const -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template get<p>(
        0))
{
    return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template get<p>(grid.LinId(coord));
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int p, typename CoordT>
inline __device__ auto SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::insert(CoordT coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template insert<p>(
        0))
{
    return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template insert<p>(grid.LinId(coord));
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<typename CoordT>
inline __device__ unsigned int SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::getBlockId(CoordT coord)
{
    return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlockId(grid.LinId(coord));
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<typename CoordT>
inline __device__ auto SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::getBlock(CoordT coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlock(0))
{
    return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlock(getBlockId(coord));
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<typename CoordT>
inline __device__ auto SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::insertBlock(CoordT coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock(0))
{
    return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock(getBlockId(coord));
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int p, typename CoordT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
//::loadBlock(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion, unsigned int stencilSupportRadius)
::loadBlock(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion)
{
    const auto block = getBlock(coord);
    __loadBlock<p>(block, sharedRegion, stencilSupportRadius);
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int p, typename CoordT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::loadGhost(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion)
{
    auto block = getBlock(coord);
    __loadGhost<p>(block, sharedRegion);
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int p, typename CoordT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
//::storeBlock(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion, unsigned int stencilSupportRadius)
::storeBlock(CoordT coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion)
{
    auto block = insertBlock(coord);
//    __storeBlock<p>(block, sharedRegion, stencilSupportRadius);
    __storeBlock<p>(block, sharedRegion);
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int ... props, typename CoordT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
//::loadBlock(CoordT coord, void *sharedRegionPtr[sizeof...(props)],
//            unsigned int stencilSupportRadius[sizeof...(props)])
::loadBlock(CoordT coord, void *sharedRegionPtr[sizeof...(props)])
{
    auto block = getBlock(coord);
//    __loadBlock<props ...>(block, sharedRegionPtr, stencilSupportRadius);
    __loadBlock<props ...>(block, sharedRegionPtr);
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int ... props, typename CoordT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::loadGhost(CoordT coord, void *sharedRegionPtr[sizeof...(props)])
{
    auto block = getBlock(coord);
    __loadGhost<props ...>(block, sharedRegionPtr);
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int ... props, typename CoordT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
//::storeBlock(CoordT coord, void *sharedRegionPtr[sizeof...(props)],
//            unsigned int stencilSupportRadius[sizeof...(props)])
::storeBlock(CoordT coord, void *sharedRegionPtr[sizeof...(props)])
{
    auto block = getBlock(coord);
//    __storeBlock<props ...>(block, sharedRegionPtr, stencilSupportRadius);
    __storeBlock<props ...>(block, sharedRegionPtr);
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KER_CUH
