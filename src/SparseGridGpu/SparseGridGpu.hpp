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
    grid_sm<dim,void> gridSize;
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

    template<typename stencil>
    void applyStencils();

    template<typename stencil, typename ... otherStencils>
    void applyStencils();
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
template<typename stencil>
void SparseGridGpu<dim, AggregateT, blockEdgeSize, threadBlockSize, indexT, layout_base>::applyStencils()
{
    //todo: Apply the given stencil on all elements which are not boundary-tagged
}

template<unsigned int dim,typename AggregateT,
		 unsigned int blockEdgeSize, unsigned int threadBlockSize,
		 typename indexT, template<typename> class layout_base>
template<typename stencil, typename ... otherStencils>
void SparseGridGpu<dim, AggregateT, blockEdgeSize, threadBlockSize, indexT, layout_base>::applyStencils()
{
    applyStencils<stencil>();
    applyStencils<otherStencils ...>();
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
