//
// Created by tommaso on 11/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_KER_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_KER_CUH

#include <SparseGridGpu/Geometry/BlockGeometry.hpp>
#include "BlockMapGpu.hpp"

//todo Remove template param GridSmT and just use BlockGeometry
template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
class SparseGridGpu_ker : public BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>
{
private:
    BlockGeometry<dim, blockEdgeSize> grid;
    GridSmT blockWithGhostGrid;

protected:
    const static unsigned char PADDING_BIT = 1;
    static constexpr unsigned int blockSize = BlockTypeOf<AggregateBlockT, 0>::size;
    // Good reason on why not to use constant memory for ghostLayer mapping is here https://stackoverflow.com/a/18021374
    unsigned int ghostLayerSize;
    unsigned int *ghostLayerToThreadsMapping;

public:
    static constexpr unsigned int d = dim;
    unsigned int stencilSupportRadius;

public:
    SparseGridGpu_ker(const openfpm::vector_sparse_gpu_ker<AggregateBlockT, indexT, layout_base> &blockMap,
                      BlockGeometry<dim, blockEdgeSize> grid,
                      GridSmT extendedBlockGeometry,
                      unsigned int stencilSupportRadius,
                      unsigned int *mapping,
                      unsigned int ghostLayerSize)
            : BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>(blockMap),
              grid(grid),
              blockWithGhostGrid(extendedBlockGeometry),
              stencilSupportRadius(stencilSupportRadius),
              ghostLayerSize(ghostLayerSize)
    {
        cudaMalloc(&ghostLayerToThreadsMapping, ghostLayerSize * sizeof(unsigned int));
        cudaMemcpy(ghostLayerToThreadsMapping, mapping, ghostLayerSize * sizeof(unsigned int), cudaMemcpyHostToDevice);
    }

    // Geometry
    template<typename CoordT>
    inline __device__ size_t getLinId(const grid_key_dx<dim, CoordT> & coord) const
    {
        return grid.LinId(coord);
    }

    inline __device__ unsigned int size(unsigned int i)
    {
    	return grid.getSize()[i];
    }

    inline __device__ grid_key_dx<dim> getCoord(size_t linId) const
    {
        return grid.InvLinId(linId);

    }

    template<typename CoordT>
    inline __device__ size_t getBlockLinId(CoordT blockCoord) const
    {
        return grid.BlockLinId(blockCoord);
    }

    inline __device__ grid_key_dx<dim> getBlockCoord(size_t blockLinId) const
    {
        return grid.BlockInvLinId(blockLinId);
    }

    inline __device__ grid_key_dx<dim> getBlockBaseCoord(size_t blockLinId) const
    {
        return grid.InvLinId(blockLinId * blockSize);
    }

    inline __device__ grid_key_dx<dim> getNeighbour(grid_key_dx<dim> base, unsigned int dimension, char offset) const
    {
        grid_key_dx<dim> res = base;
        auto i = base.get(dimension) + offset;
        res.set_d(dimension, i);
        return res;
    }

    constexpr __device__ unsigned int getBlockEdgeSize() const
    {
        return blockEdgeSize;
    }

    inline __device__ unsigned int getEnlargedBlockSize() const
    {
        return std::pow(blockEdgeSize + 2*stencilSupportRadius, dim);
    }

    inline __device__ unsigned int posToEnlargedBlockPos(unsigned int pos) const
    {
        // Convert pos into a linear id accounting for the ghost offsets
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(pos, stencilSupportRadius, coord);
        const unsigned int linId = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
        return linId;
    }

    inline __device__ grid_key_dx<dim>
    getCoordInEnlargedBlock(const unsigned int offset) const
    {
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(offset, stencilSupportRadius, coord);
        return grid_key_dx<dim>(coord);
    }

    inline __device__ unsigned int
    getLinIdInEnlargedBlock(const unsigned int offset) const
    {
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(offset, stencilSupportRadius, coord);
        return coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
    }

    inline __device__ unsigned int
    getNeighbourLinIdInEnlargedBlock(grid_key_dx<dim> base, unsigned int dimension, char offset) const
    {
        grid_key_dx<dim> res = getNeighbour(base, dimension, offset);
        return coordToLin<blockEdgeSize>(res, stencilSupportRadius);
    }

    // Data management methods

    template<unsigned int p, typename CoordT>
    inline __device__ auto
    get(const grid_key_dx<dim, CoordT> & coord) const -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template get<p>(
            0));

    template<unsigned int p, typename CoordT>
    inline __device__ auto
    insert(const grid_key_dx<dim, CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template insert<p>(
            0));

    template<typename CoordT>
    inline __device__ unsigned int getBlockId(const grid_key_dx<dim, CoordT> & coord);

    template<typename CoordT>
    inline __device__ auto
    getBlock(const grid_key_dx<dim, CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlock(0));

    template<typename CoordT>
    inline __device__ auto
    insertBlock(const grid_key_dx<dim, CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock(0));

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
//    loadBlock(const grid_key_dx<dim, CoordT> coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion);

    template<unsigned int p, typename AggrWrapperT>
    inline __device__ void
    loadBlock(AggrWrapperT &block, ScalarTypeOf<AggregateBlockT, p> *sharedRegion);

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
//    inline __device__ void loadBlock(const grid_key_dx<dim, CoordT> coord, void *sharedRegionPtr[sizeof...(props)]);

    template<unsigned int ... props, typename AggrWrapperT>
    inline __device__ void loadBlock(AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)]);

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
    loadGhost(const grid_key_dx<dim, CoordT> & coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion);

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
    inline __device__ void loadGhost(const grid_key_dx<dim, CoordT> & coord, void *sharedRegionPtr[sizeof...(props)]);

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
//    storeBlock(const grid_key_dx<dim, CoordT> coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion);

    template<unsigned int p, typename AggrWrapperT>
    inline __device__ void
    storeBlock(AggrWrapperT &block, ScalarTypeOf<AggregateBlockT, p> *sharedRegion);

    template<unsigned int p, typename CoordT>
    inline __device__ void
    storeBlockInPlace(const grid_key_dx<dim, CoordT> & coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion);

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
//    inline __device__ void storeBlock(const grid_key_dx<dim, CoordT> coord, void *sharedRegionPtr[sizeof...(props)]);

    template<unsigned int ... props, typename AggrWrapperT>
    inline __device__ void storeBlock(AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)]);

    template<unsigned int ... props, typename CoordT>
    inline __device__ void storeBlockInPlace(const grid_key_dx<dim, CoordT> & coord, void *sharedRegionPtr[sizeof...(props)]);

    template <unsigned int p, typename CoordT>
    inline __device__ ScalarTypeOf<AggregateBlockT, p> & get(
            grid_key_dx<dim, CoordT> coord,
            Box<dim, indexT> sharedMemBox,
            ScalarTypeOf<AggregateBlockT, p> *sharedRegion)
    {
        //NOTE: Size of Box must be equal to size of shared region!
        //NOTE: Box must be square
        if (sharedMemBox.isInside(coord.toPoint()))
        {
            //Get data from shared mem
            auto one = coord;
            one.one();
            const auto boxDimensions = sharedMemBox.getKP2() - sharedMemBox.getKP1() + one; // The +1 is because upper bound is inclusive
            const auto relativeCoord = coord - sharedMemBox.getKP1();
            const auto locLinId = coordToLin(relativeCoord, boxDimensions);
            return sharedRegion[locLinId];
        }
        else
        {
            //Get data from global mem
            return get(coord);
        }
    }

    template<typename BitMaskT>
    inline static __device__ bool isPadding(const BitMaskT &bitMask)
    {
        return getBit(bitMask, PADDING_BIT);
    }

    template <typename keyIndexT>
    inline __device__ bool isPadding(grid_key_dx<dim, keyIndexT> coord) const
    {
        auto mask = getMask(grid.LinId(coord));
        return isPadding(mask);
    }

    template<typename BitMaskT>
    inline static __device__ void setPadding(BitMaskT &bitMask)
    {
        setBit(bitMask, PADDING_BIT);
    }

    template<typename BitMaskT>
    inline static __device__ void unsetPadding(BitMaskT &bitMask)
    {
        unsetBit(bitMask, PADDING_BIT);
    }

private:
    template<unsigned int edgeSize>
    inline __device__ void linToCoordWithOffset(const unsigned int linId, const unsigned int offset, unsigned int (&coord)[dim]) const
    {
        unsigned int linIdTmp = linId;
        for (unsigned int d = 0; d < dim; ++d)
        {
            coord[d] = linIdTmp % edgeSize;
            coord[d] += offset;
            linIdTmp /= edgeSize;
        }
    }

    template<unsigned int edgeSize>
    inline __device__ unsigned int coordToLin(const unsigned int (&coord)[dim], const unsigned int paddingSize = 0) const
    {
        unsigned int linId = coord[dim - 1];
        for (int d = dim - 2; d >= 0; --d)
        {
            linId *= edgeSize + 2 * paddingSize;
            linId += coord[d];
        }
        return linId;
    }

    template<unsigned int edgeSize>
    inline __device__ unsigned int coordToLin(const grid_key_dx<dim> &coord, const unsigned int paddingSize = 0) const
    {
        unsigned int linId = coord.get(dim - 1);
        for (int d = dim - 2; d >= 0; --d)
        {
            linId *= edgeSize + 2 * paddingSize;
            linId += coord.get(d);
        }
        return linId;
    }

    inline __device__ unsigned int coordToLin(const grid_key_dx<dim> & coord, grid_key_dx<dim> & blockDimensions) const
    {
        unsigned int linId = coord.get(dim - 1);
        for (int d = dim - 2; d >= 0; --d)
        {
            linId *= blockDimensions.get(d);
            linId += coord.get(d);
        }
        return linId;
    }

    template<unsigned int p, typename AggrWrapperT, typename SharedPtrT>
    inline __device__ void
    __loadBlock(const AggrWrapperT &block, SharedPtrT sharedRegionPtr)
    {
        typedef ScalarTypeOf<AggregateBlockT, p> ScalarT;

        const unsigned int pos = threadIdx.x;
        //todo: Improve this version to allow multiple chunks per block!
        if (pos < blockSize)
        {
            // Convert pos into a linear id accounting for the ghost offsets
            unsigned int coord[dim];
            linToCoordWithOffset<blockEdgeSize>(pos, stencilSupportRadius, coord);
            const unsigned int linId = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
            // Actually load the data into the shared region
            ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
            *(basePtr + linId) = block.template get<p>()[pos];
        }
    }

    template<unsigned int p, unsigned int ... props, typename AggrWrapperT>
    inline __device__ void
    __loadBlock(const AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)+1])
    {
        __loadBlock<p>(block, sharedRegionPtr);
        if (sizeof...(props) > 1)
        {
            __loadBlock<props ...>(block, sharedRegionPtr + 1);
        }
        else if (sizeof...(props) == 1)
        {
            __loadBlock<props ...>(block, *(sharedRegionPtr + 1));
        }
    }

    // NOTE: this must be called with linear thread grid, nice-to-have would be a generic converter (easy to do)
    // from dim3 to linear which would work under all possible launch params
    template<unsigned int p, typename SharedPtrT>
    inline __device__ void
    __loadGhost(const unsigned int blockId, SharedPtrT sharedRegionPtr)
    {
        typedef ScalarTypeOf<AggregateBlockT, p> ScalarT;

        const unsigned int rounds = ghostLayerSize % blockDim.x == 0 // ceiling
                                    ? ghostLayerSize / blockDim.x
                                    : 1 + ghostLayerSize / blockDim.x;
        for (int round = 0; round < rounds; ++round)
        {
            const unsigned int pos = round * blockDim.x + threadIdx.x;
            if (pos < ghostLayerSize)
            {
                // Convert pos into a linear id accounting for the inner domain offsets
                const unsigned int linId = ghostLayerToThreadsMapping[pos];
                // Now get linear offset wrt the first element of the block
                grid_key_dx<dim> localCoord = blockWithGhostGrid.InvLinId(linId);
                grid_key_dx<dim> elementCoord = getBlockBaseCoord(blockId);
                for (int i = 0; i < dim; ++i)
                {
                    localCoord.set_d(i, localCoord.get(i) - stencilSupportRadius);
                    elementCoord.set_d(i, elementCoord.get(i) - localCoord.get(i));
                }

                // Actually load the data into the shared region
                ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
                *(basePtr + linId) = get<p>(elementCoord);
            }
        }
    }

    template<unsigned int p, unsigned int ... props>
    inline __device__ void
    __loadGhost(const unsigned int blockId, void *sharedRegionPtr[sizeof...(props)+1])
    {
        __loadGhost<p>(blockId, sharedRegionPtr);
        if (sizeof...(props) > 1)
        {
            __loadGhost<props ...>(blockId, sharedRegionPtr + 1);
        }
        else if (sizeof...(props) == 1)
        {
            __loadGhost<props ...>(blockId, *(sharedRegionPtr + 1));
        }
    }

    template<unsigned int p, typename AggrWrapperT, typename SharedPtrT>
    inline __device__ void
    __storeBlock(AggrWrapperT &block, SharedPtrT sharedRegionPtr)
    {
        typedef ScalarTypeOf<AggregateBlockT, p> ScalarT;

        const unsigned int pos = threadIdx.x;
        //todo: Improve this version to allow multiple chunks per block!
        if (pos < blockSize)
        {
            // Convert pos into a linear id accounting for the ghost offsets
            unsigned int coord[dim];
            linToCoordWithOffset<blockEdgeSize>(pos, stencilSupportRadius, coord);
            const unsigned int linId = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
            // Actually store the data from the shared region
            ScalarT *basePtr = (ScalarT *)sharedRegionPtr;

            block.template get<p>()[pos] = *(basePtr + linId);
        }
    }

    template<unsigned int p, unsigned int ... props, typename AggrWrapperT>
    inline __device__ void
    __storeBlock(AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)+1])
    {
        __storeBlock<p>(block, sharedRegionPtr);
        if (sizeof...(props) > 1)
        {
            __storeBlock<props ...>(block, sharedRegionPtr + 1);
        }
        else if (sizeof...(props) == 1)
        {
            __storeBlock<props ...>(block, *(sharedRegionPtr + 1));
        }
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
::get(const grid_key_dx<dim, CoordT> & coord) const -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template get<p>(
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
::insert(const grid_key_dx<dim, CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template insert<p>(
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
inline __device__ unsigned int
SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::getBlockId(const grid_key_dx<dim, CoordT> & coord)
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
::getBlock(const grid_key_dx<dim, CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlock(0))
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
::insertBlock(const grid_key_dx<dim, CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock(0))
{
    return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock(getBlockId(coord));
}

//template<unsigned int dim,
//        unsigned int blockEdgeSize,
//        typename AggregateBlockT,
//        typename indexT,
//        template<typename> class layout_base,
//        typename GridSmT>
//template<unsigned int p, typename CoordT>
//inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
//::loadBlock(const grid_key_dx<dim, CoordT> coord, ScalarTypeOf <AggregateBlockT, p> *sharedRegion)
//{
//    //todo: Make this work well with multiples chunks per block or check not to get several chunks or dragons ahoy!
//    auto block = getBlock(coord);
//    __loadBlock<p>(block, sharedRegion);
//}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int p, typename AggrWrapperT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::loadBlock(AggrWrapperT &block, ScalarTypeOf <AggregateBlockT, p> *sharedRegion)
{
    //todo: Make this work well with multiples chunks per block or check not to get several chunks or dragons ahoy!
    __loadBlock<p>(block, sharedRegion);
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int p, typename CoordT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::loadGhost(const grid_key_dx<dim, CoordT> & coord, ScalarTypeOf <AggregateBlockT, p> *sharedRegion)
{
    //todo: Make this work well with multiples chunks per block or check not to get several chunks or dragons ahoy!
    auto blockLinId = getBlockId(coord);
    __loadGhost<p>(blockLinId, sharedRegion);
}

//template<unsigned int dim,
//        unsigned int blockEdgeSize,
//        typename AggregateBlockT,
//        typename indexT,
//        template<typename> class layout_base,
//        typename GridSmT>
//template<unsigned int p, typename CoordT>
//inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
//::storeBlock(const grid_key_dx<dim, CoordT> coord, ScalarTypeOf <AggregateBlockT, p> *sharedRegion)
//{
//    //todo: Make this work well with multiples chunks per block or check not to get several chunks or dragons ahoy!
//    auto block = insertBlock(coord);
//    __storeBlock<p>(block, sharedRegion);
//}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int p, typename AggrWrapperT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::storeBlock(AggrWrapperT &block, ScalarTypeOf <AggregateBlockT, p> *sharedRegion)
{
    //todo: Make this work well with multiples chunks per block or check not to get several chunks or dragons ahoy!
    __storeBlock<p>(block, sharedRegion);
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int p, typename CoordT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::storeBlockInPlace(const grid_key_dx<dim, CoordT> & coord, ScalarTypeOf <AggregateBlockT, p> *sharedRegion)
{
    //todo: Make this work well with multiples chunks per block or check not to get several chunks or dragons ahoy!
    auto & block = getBlock(coord);
    __storeBlock<p>(block, sharedRegion);
}

//template<unsigned int dim,
//        unsigned int blockEdgeSize,
//        typename AggregateBlockT,
//        typename indexT,
//        template<typename> class layout_base,
//        typename GridSmT>
//template<unsigned int ... props, typename CoordT>
//inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
//::loadBlock(const grid_key_dx<dim, CoordT> coord, void *sharedRegionPtr[sizeof...(props)])
//{
//    auto block = getBlock(coord);
//    __loadBlock<props ...>(block, sharedRegionPtr);
//}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int ... props, typename AggrWrapperT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::loadBlock(AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)])
{
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
::loadGhost(const grid_key_dx<dim, CoordT> & coord, void *sharedRegionPtr[sizeof...(props)])
{
    auto blockLinId = getBlockId(coord);
    __loadGhost<props ...>(blockLinId, sharedRegionPtr);
}

//template<unsigned int dim,
//        unsigned int blockEdgeSize,
//        typename AggregateBlockT,
//        typename indexT,
//        template<typename> class layout_base,
//        typename GridSmT>
//template<unsigned int ... props, typename CoordT>
//inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
//::storeBlock(const grid_key_dx<dim, CoordT> coord, void *sharedRegionPtr[sizeof...(props)])
//{
//    auto block = insertBlock(coord);
//    __storeBlock<props ...>(block, sharedRegionPtr);
//}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int ... props, typename AggrWrapperT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::storeBlock(AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)])
{
    __storeBlock<props ...>(block, sharedRegionPtr);
}

template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT>
template<unsigned int ... props, typename CoordT>
inline __device__ void SparseGridGpu_ker<dim, blockEdgeSize, AggregateBlockT, indexT, layout_base, GridSmT>
::storeBlockInPlace(const grid_key_dx<dim, CoordT> & coord, void *sharedRegionPtr[sizeof...(props)])
{
    auto block = getBlock(coord);
    __storeBlock<props ...>(block, sharedRegionPtr);
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KER_CUH
