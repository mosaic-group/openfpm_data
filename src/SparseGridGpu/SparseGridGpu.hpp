//
// Created by tommaso on 6/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_HPP
#define OPENFPM_PDATA_SPARSEGRIDGPU_HPP

#include <cstdlib>
#include <SparseGridGpu/BlockMapGpu.hpp>
#include <Grid/iterators/grid_skin_iterator.hpp>
#include <SparseGridGpu/Geometry/BlockGeometry.hpp>
#include "SparseGridGpu_ker.cuh"
#include "SparseGridGpu_kernels.cuh"

template<unsigned int dim, unsigned int blockEdgeSize, typename AggregateBlockT, unsigned int threadBlockSize = 128, typename indexT=int, template<typename> class layout_base=memory_traits_inte>
class SparseGridGpu : public BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>
{
private:
    BlockGeometry<dim, blockEdgeSize> gridGeometry;
    grid_sm<dim, int> extendedBlockGeometry;
    unsigned int stencilSupportRadius;
    unsigned int ghostLayerSize;
    unsigned int *ghostLayerToThreadsMapping;

protected:
    static constexpr unsigned int blockSize = BlockTypeOf<AggregateBlockT, 0>::size;

public:
    template<typename dim3T>
    inline static int dim3SizeToInt(dim3T d)
    {
        return d.x * d.y * d.z;
    }

    inline static int dim3SizeToInt(unsigned int d)
    {
        return d;
    }

private:
    void computeSizeOfGhostLayer()
    {
        unsigned int term1 = 1;
        for (int i = 0; i < dim; ++i)
        {
            term1 *= blockEdgeSize + 2 * stencilSupportRadius;
        }
        unsigned int term2 = 1;
        for (int i = 0; i < dim; ++i)
        {
            term2 *= blockEdgeSize;
        }
        ghostLayerSize = term1 - term2;
    }

    void allocateGhostLayerMapping()
    {
        ghostLayerToThreadsMapping = new unsigned int[ghostLayerSize];
    }

    void computeGhostLayerMapping()
    {
//        std::cout << "blockEdgeSize=" << blockEdgeSize << ", stencilSupportRadius=" << stencilSupportRadius << std::endl; //debug
        size_t dimensions[dim],
                origin[dim],
                innerDomainBegin[dim], innerDomainEnd[dim],
                outerBoxBegin[dim], outerBoxEnd[dim],
                bc[dim];
        for (int i = 0; i < dim; ++i)
        {
            dimensions[i] = blockEdgeSize + 2 * stencilSupportRadius;
            origin[i] = 0;
            innerDomainBegin[i] = stencilSupportRadius - 1;
            innerDomainEnd[i] = dimensions[i] - stencilSupportRadius;
            outerBoxBegin[i] = origin[i];
            outerBoxEnd[i] = dimensions[i];
            bc[i] = NON_PERIODIC;
        }
        grid_sm<dim, void> enlargedGrid;
        enlargedGrid.setDimensions(dimensions);
        Box<dim, size_t> outerBox(outerBoxBegin, outerBoxEnd);
        Box<dim, size_t> innerBox(innerDomainBegin, innerDomainEnd);

        grid_skin_iterator_bc<dim> gsi(enlargedGrid, innerBox, outerBox, bc);

        unsigned int i = 0;
        unsigned int limit = ghostLayerSize;
        while (gsi.isNext())
        {
            auto coord = gsi.get();
//            std::cout << "coord[" << i << "] = {" << coord.get(0) << "," << coord.get(1) << "}" << std::endl; //debug
            assert(i < ghostLayerSize);
            ghostLayerToThreadsMapping[i] = enlargedGrid.LinId(coord);
            ++i;
            ++gsi;
        }
        assert(i == limit);
//        std::cout << "i=" << i << ", limit=" << limit << std::endl; //debug
    }

public:
    SparseGridGpu(BlockGeometry<dim, blockEdgeSize> gridGeometry, unsigned int stencilSupportRadius = 1)
            : gridGeometry(gridGeometry),
              stencilSupportRadius(stencilSupportRadius)
    {
        computeSizeOfGhostLayer();
        allocateGhostLayerMapping();
        computeGhostLayerMapping();

        size_t extBlockDims[dim];
        for (int d=0; d<dim; ++d)
        {
            extBlockDims[d] = blockEdgeSize + stencilSupportRadius;
        }
        extendedBlockGeometry.setDimensions(extBlockDims);
    };

    virtual ~SparseGridGpu();

    SparseGridGpu_ker
            <
                    dim,
                    blockEdgeSize,
                    typename BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                    indexT,
                    layout_base,
                    decltype(gridGeometry),
                    decltype(extendedBlockGeometry)
            > toKernel()
    {
        SparseGridGpu_ker
                <
                        dim,
                        blockEdgeSize,
                        typename BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                        indexT,
                        layout_base,
                        decltype(gridGeometry),
                        decltype(extendedBlockGeometry)
                > toKer(
                BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::blockMap.toKernel(),
                gridGeometry,
                extendedBlockGeometry,
                stencilSupportRadius,
                ghostLayerToThreadsMapping,
                ghostLayerSize);
        return toKer;
    }

    // Geometry
    template<typename CoordT>
    inline size_t getLinId(CoordT &coord)
    {
        return gridGeometry.LinId(coord);
    }

    inline grid_key_dx<dim> getCoord(size_t linId)
    {
        return gridGeometry.InvLinId(linId);

    }

    // Data management methods

    template<unsigned int p, typename CoordT>
    auto
    get(CoordT &coord) const -> ScalarTypeOf<AggregateBlockT, p>;

    template<unsigned int p, typename CoordT>
    auto
    insert(CoordT &coord) -> ScalarTypeOf<AggregateBlockT, p> &;

    template<typename dim3T>
    void setGPUInsertBuffer(dim3T nBlock, dim3T nSlot)
    {
        BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>
        ::setGPUInsertBuffer(
                dim3SizeToInt(nBlock),
                dim3SizeToInt(nSlot)
        );
    }

    // Stencil-related methods
    void tagBoundaries();

    template<typename stencil>
    void applyStencils();

    template<typename stencil, typename ... otherStencils>
    void applyStencils();
};

template<unsigned int dim, unsigned int blockEdgeSize, typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
void SparseGridGpu<dim, blockEdgeSize, AggregateBlockT, threadBlockSize, indexT, layout_base>::tagBoundaries()
{
    //todo: Here iterate on all existing elements and tag those which are at the edge
    // Notes:
    //  1) On host we need some kind of iterator to be provided from the blockMap
    //  2) On GPU we need a way to get the device arrays of block keys and blocks from the underlying data structure
    auto indexBuffer = BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
    auto dataBuffer = BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();
    unsigned int gridSize = indexBuffer.size()%threadBlockSize==0
                                ? indexBuffer.size() / threadBlockSize
                                : 1 + indexBuffer.size() / threadBlockSize;
    SparseGridGpuKernels::tagBoundaries<<<gridSize, threadBlockSize>>>(indexBuffer.toKernel(), dataBuffer.toKernel(), this->toKernel());
}

template<unsigned int dim, unsigned int blockEdgeSize, typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<typename stencil>
void SparseGridGpu<dim, blockEdgeSize, AggregateBlockT, threadBlockSize, indexT, layout_base>::applyStencils()
{
    //todo: Apply the given stencil on all elements which are not boundary-tagged
}

template<unsigned int dim, unsigned int blockEdgeSize, typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<typename stencil, typename ... otherStencils>
void SparseGridGpu<dim, blockEdgeSize, AggregateBlockT, threadBlockSize, indexT, layout_base>::applyStencils()
{
    applyStencils<stencil>();
    applyStencils<otherStencils ...>();
}

template<unsigned int dim, unsigned int blockEdgeSize, typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<unsigned int p, typename CoordT>
auto SparseGridGpu<dim, blockEdgeSize, AggregateBlockT, threadBlockSize, indexT, layout_base>::
get(CoordT & coord) const -> ScalarTypeOf<AggregateBlockT, p>
{
    return BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::template get<p>(gridGeometry.LinId(coord));
}

template<unsigned int dim, unsigned int blockEdgeSize, typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<unsigned int p, typename CoordT>
auto SparseGridGpu<dim, blockEdgeSize, AggregateBlockT, threadBlockSize, indexT, layout_base>::insert(
        CoordT &coord) -> ScalarTypeOf<AggregateBlockT, p> &
{
    return BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::template insert<p>(gridGeometry.LinId(coord));
}

template<unsigned int dim, unsigned int blockEdgeSize, typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
SparseGridGpu<dim, blockEdgeSize, AggregateBlockT, threadBlockSize, indexT, layout_base>::~SparseGridGpu()
{
    delete[](ghostLayerToThreadsMapping);
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_HPP
