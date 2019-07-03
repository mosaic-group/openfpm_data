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
#include "Iterators/SparseGridGpu_iterator_sub.hpp"

// todo: Move all the following utils into some proper file inside TemplateUtils

template<unsigned int dim>
struct default_edge
{
	typedef boost::mpl::int_<2> type;
};


template<>
struct default_edge<1>
{
	typedef boost::mpl::int_<64> type;
};

template<>
struct default_edge<2>
{
	typedef boost::mpl::int_<8> type;
};

template<>
struct default_edge<3>
{
	typedef boost::mpl::int_<4> type;
};

template<unsigned int dim, unsigned int blockEdgeSize, typename ... aggr_list>
struct aggregate_transform_datablock_impl
{
	typedef aggregate<DataBlock<aggr_list,IntPow<blockEdgeSize,dim>::value>...> type;
};

template<unsigned int dim, unsigned int blockEdgeSize, typename aggr>
struct aggregate_convert
{
};

template<unsigned int dim, unsigned int blockEdgeSize, typename ... types>
struct aggregate_convert<dim,blockEdgeSize,aggregate<types ...>>
{
	typedef typename aggregate_transform_datablock_impl<dim,blockEdgeSize,types ...>::type type;
};

/////////////

enum StencilMode
{
    STENCIL_MODE_INSERT = 0,
    STENCIL_MODE_INPLACE = 1
};

template<unsigned int dim,
		 typename AggregateT,
		 unsigned int blockEdgeSize = default_edge<dim>::type::value,
		 unsigned int threadBlockSize = 128,
		 typename indexT=int,
		 template<typename> class layout_base=memory_traits_inte>
class SparseGridGpu : public BlockMapGpu<typename aggregate_convert<dim,blockEdgeSize,AggregateT>::type, threadBlockSize, indexT, layout_base>
{
private:
    const static unsigned char PADDING_BIT = 1;
	typedef typename aggregate_convert<dim,blockEdgeSize,AggregateT>::type AggregateBlockT;
    BlockGeometry<dim, blockEdgeSize> gridGeometry;
    grid_sm<dim, int> extendedBlockGeometry;
    grid_sm<dim, int> gridSize;
    unsigned int stencilSupportRadius;
    unsigned int ghostLayerSize;
    unsigned int *ghostLayerToThreadsMapping;

protected:
    static constexpr unsigned int blockSize = BlockTypeOf<AggregateBlockT, 0>::size;

public:
    static constexpr unsigned int dims = dim;

    typedef grid_key_dx<dim, indexT> base_key;

    /*! \brief return the size of the grid
     *
     * \return Return the size of the grid
     *
     */
    inline size_t size() const
    {
        return gridGeometry.size();
    }

    /*! \brief This is a meta-function return which type of sub iterator a grid produce
     *
     * \return the type of the sub-grid iterator
     *
     */
    template <typename stencil = no_stencil>
    static SparseGridGpu_iterator_sub<dim> type_of_subiterator()
    {
        return SparseGridGpu_iterator_sub<dim>();
    }

    /*! \brief This is a meta-function return which type of iterator a grid produce
     *
     * \return the type of the sub-grid iterator
     *
     */
    static SparseGridGpu_iterator type_of_iterator()
    {
        return SparseGridGpu_iterator();
    }

    template<typename dim3T>
    inline static int dim3SizeToInt(dim3T d)
    {
        return d.x * d.y * d.z;
    }

    inline static int dim3SizeToInt(size_t d)
    {
        return d;
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
        while (gsi.isNext())
        {
            auto coord = gsi.get();
            assert(i < ghostLayerSize);
            ghostLayerToThreadsMapping[i] = enlargedGrid.LinId(coord);
            ++i;
            ++gsi;
        }
        assert(i == ghostLayerSize);
    }

    void initialize(const size_t (& res)[dim])
    {
    	// calculate the number of blocks on each direction
    	size_t blocks[dim];

    	for (size_t i = 0 ; i < dim ; i++)
    	{
    		blocks[i] = res[i] / blockEdgeSize;

    		if (res[i] % blockEdgeSize != 0)
    		{blocks[i]++;}
    	}

    	gridGeometry = BlockGeometry<dim,blockEdgeSize>(blocks);

        computeSizeOfGhostLayer();
        allocateGhostLayerMapping();
        computeGhostLayerMapping();

        size_t extBlockDims[dim];
        for (int d=0; d<dim; ++d)
        {
            extBlockDims[d] = blockEdgeSize + stencilSupportRadius;
        }
        extendedBlockGeometry.setDimensions(extBlockDims);
        gridSize.setDimensions(res);
    }

    template <typename stencil, typename... Args>
    void applyStencilInPlace(Args... args)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;

        // NOTE: Here we want to work only on one data chunk per block!
        constexpr unsigned int chunksPerBlock = 1;
        const unsigned int localThreadBlockSize = dataChunkSize * chunksPerBlock;
        const unsigned int threadGridSize = numScalars % localThreadBlockSize == 0
                                            ? numScalars / localThreadBlockSize
                                            : 1 + numScalars / localThreadBlockSize;

        SparseGridGpuKernels::applyStencil
                <dim,
                BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::pMask,
                stencil>
                <<<threadGridSize, localThreadBlockSize>>>(
                        indexBuffer.toKernel(),
                        dataBuffer.toKernel(),
                        this->toKernel(),
                        true,
                        args...);
    }

    //todo: the applyInsert should also allocate the gpu insert buffer, initialize it, etc
    template <typename stencil, typename... Args>
    void applyStencilInsert(Args... args)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;

        // NOTE: Here we want to work only on one data chunk per block!
        constexpr unsigned int chunksPerBlock = 1;
        const unsigned int localThreadBlockSize = dataChunkSize * chunksPerBlock;
        const unsigned int threadGridSize = numScalars % localThreadBlockSize == 0
                                      ? numScalars / localThreadBlockSize
                                      : 1 + numScalars / localThreadBlockSize;

        setGPUInsertBuffer(threadGridSize, chunksPerBlock); // Here 1 is because we process 1 chunk per block!
        this->initializeGPUInsertBuffer();

        SparseGridGpuKernels::applyStencil
                <dim,
                BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::pMask,
                stencil>
                <<<threadGridSize, localThreadBlockSize>>>(
                        indexBuffer.toKernel(),
                        dataBuffer.toKernel(),
                        this->toKernel(),
                        true,
                        args...);

        mgpu::ofp_context_t ctx;
        stencil::flush(*this, ctx);
    }

public:

    SparseGridGpu()
	:stencilSupportRadius(1)
    {};

    /*! \brief resize the SparseGrid
     *
     * \param res indicate the resolution in each dimension
     *
     */
    void resize(size_t (& res)[dim])
    {
    	initialize(res);
    }

    /*! \brief
     *
     *
     */
    SparseGridGpu(BlockGeometry<dim, blockEdgeSize> & gridGeometry, unsigned int stencilSupportRadius = 1)
            : gridGeometry(gridGeometry),
              stencilSupportRadius(stencilSupportRadius)
    {
    	initialize(gridGeometry.getSize());
    };

    virtual ~SparseGridGpu();

    SparseGridGpu_ker
            <
                    dim,
                    blockEdgeSize,
                    typename BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                    indexT,
                    layout_base,
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

    inline ite_gpu<dim> getGridGPUIterator(const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop, size_t n_thr = threadBlockSize)
    {
    	return gridSize.getGPUIterator(start,stop,n_thr);
    }

    // Data management methods

    template<unsigned int p, typename CoordT>
    auto
    get(const CoordT &coord) const -> const ScalarTypeOf<AggregateBlockT, p> &;

    template<unsigned int p, typename CoordT>
    auto
    insert(const CoordT &coord) -> ScalarTypeOf<AggregateBlockT, p> &;

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

    template<typename stencil, typename... Args>
    void applyStencils(StencilMode mode = STENCIL_MODE_INSERT, Args... args);

    template<typename stencil1, typename stencil2, typename ... otherStencils, typename... Args>
    void applyStencils(StencilMode mode = STENCIL_MODE_INSERT, Args... args);

    template<typename BitMaskT>
    inline static bool isPadding(BitMaskT &bitMask)
    {
        return getBit(bitMask, PADDING_BIT);
    }

    template<typename BitMaskT>
    inline static void setPadding(BitMaskT &bitMask)
    {
        setBit(bitMask, PADDING_BIT);
    }

    template<typename BitMaskT>
    inline static void unsetPadding(BitMaskT &bitMask)
    {
        unsetBit(bitMask, PADDING_BIT);
    }
};

template<unsigned int dim, typename AggregateT,
		 unsigned int blockEdgeSize, unsigned int threadBlockSize,
		 typename indexT, template<typename> class layout_base>
void SparseGridGpu<dim, AggregateT, blockEdgeSize, threadBlockSize, indexT, layout_base>::tagBoundaries()
{
    // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
    auto & indexBuffer = BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
    auto & dataBuffer = BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

    const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
    unsigned int numScalars = indexBuffer.size() * dataChunkSize;

    if (numScalars == 0) return;

    // NOTE: Here we want to work only on one data chunk per block!

    unsigned int localThreadBlockSize = dataChunkSize;
    unsigned int threadGridSize = numScalars % dataChunkSize == 0
                                ? numScalars / dataChunkSize
                                : 1 + numScalars / dataChunkSize;
    if (stencilSupportRadius == 1)
    {
        SparseGridGpuKernels::tagBoundaries<
                dim,
                1,
                BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::pMask>
                <<< threadGridSize, localThreadBlockSize >>> (indexBuffer.toKernel(), dataBuffer.toKernel(), this->toKernel());
    }
    else if (stencilSupportRadius == 2)
    {
        SparseGridGpuKernels::tagBoundaries<
                dim,
                2,
                BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::pMask>
                <<<threadGridSize, localThreadBlockSize>>>(indexBuffer.toKernel(), dataBuffer.toKernel(), this->toKernel());
    }
    else if (stencilSupportRadius == 0)
    {
        SparseGridGpuKernels::tagBoundaries<
                dim,
                0,
                BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::pMask>
                <<<threadGridSize, localThreadBlockSize>>>(indexBuffer.toKernel(), dataBuffer.toKernel(), this->toKernel());
    }
    else
    {
        //todo: KABOOOOOOM!!!
        std::cout << __FILE__ << ":" << __LINE__ << " error: stencilSupportRadius supported only up to 2, passed: " << stencilSupportRadius << std::endl;

    }
}

template<unsigned int dim, typename AggregateT,
		 unsigned int blockEdgeSize, unsigned int threadBlockSize,
		 typename indexT, template<typename> class layout_base>
template<typename stencil, typename... Args>
void SparseGridGpu<dim, AggregateT, blockEdgeSize, threadBlockSize, indexT, layout_base>
        ::applyStencils(StencilMode mode, Args... args)
{
    // Apply the given stencil on all elements which are not boundary-tagged
    // The idea is to have this function launch a __global__ function (made by us) on all existing blocks
    // then this kernel checks if elements exist && !padding and on those it calls the user-provided
    // __device__ functor. The mode of the stencil application is used basically to choose how we load the block
    // that we pass to the user functor as storeBlock: in case of Insert, we get the block through an insert (and
    // we also call the necessary aux functions); in case of an In-place we just get the block from the data buffer.
    switch (mode)
    {
        case STENCIL_MODE_INPLACE:
            applyStencilInPlace<stencil>(args...);
            break;
        case STENCIL_MODE_INSERT:
            applyStencilInsert<stencil>(args...);
            break;
    }
}

template<unsigned int dim,typename AggregateT,
		 unsigned int blockEdgeSize, unsigned int threadBlockSize,
		 typename indexT, template<typename> class layout_base>
template<typename stencil1, typename stencil2, typename ... otherStencils, typename... Args>
void SparseGridGpu<dim, AggregateT, blockEdgeSize, threadBlockSize, indexT, layout_base>
        ::applyStencils(StencilMode mode, Args... args)
{
    applyStencils<stencil1>(mode, args...);
    applyStencils<stencil2, otherStencils ...>(mode, args...);
}

template<unsigned int dim, typename AggregateT,
		 unsigned int blockEdgeSize, unsigned int threadBlockSize,
		 typename indexT, template<typename> class layout_base>
template<unsigned int p, typename CoordT>
auto SparseGridGpu<dim, AggregateT, blockEdgeSize, threadBlockSize, indexT, layout_base>::
get(const CoordT & coord) const -> const ScalarTypeOf<AggregateBlockT, p> &
{
    return BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::template get<p>(gridGeometry.LinId(coord));
}

template<unsigned int dim, typename AggregateT, unsigned int blockEdgeSize, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<unsigned int p, typename CoordT>
auto SparseGridGpu<dim, AggregateT, blockEdgeSize, threadBlockSize, indexT, layout_base>::insert(
        const CoordT &coord) -> ScalarTypeOf<AggregateBlockT, p> &
{
    return BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::template insert<p>(gridGeometry.LinId(coord));
}

template<unsigned int dim, typename AggregateT, unsigned int blockEdgeSize, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
SparseGridGpu<dim, AggregateT, blockEdgeSize, threadBlockSize, indexT, layout_base>::~SparseGridGpu()
{
    delete[](ghostLayerToThreadsMapping);
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_HPP
