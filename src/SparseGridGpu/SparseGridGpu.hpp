//
// Created by tommaso on 6/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_HPP
#define OPENFPM_PDATA_SPARSEGRIDGPU_HPP

#include <cstdlib>
#include <SparseGridGpu/BlockMapGpu.hpp>
#include <Grid/iterators/grid_skin_iterator.hpp>
#include <SparseGridGpu/Geometry/grid_smb.hpp>
#include "SparseGridGpu_ker.cuh"
#include "SparseGridGpu_kernels.cuh"
#include "Iterators/SparseGridGpu_iterator_sub.hpp"
#include "Geometry/grid_zmb.hpp"
#include "util/stat/common_statistics.hpp"
#ifdef OPENFPM_DATA_ENABLE_IO_MODULE
#include "VTKWriter/VTKWriter.hpp"
#endif

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

template<typename aggr>
struct aggregate_add
{
};

template<typename ... types>
struct aggregate_add<aggregate<types ...>>
{
	typedef aggregate<types ..., unsigned char> type;
};

/////////////

enum StencilMode
{
    STENCIL_MODE_INSERT = 0,
    STENCIL_MODE_INPLACE = 1
};

struct NNFull
{
	static const int nNN = 27;
};

struct NNStar
{
	static const int nNN = 8;
};

template<unsigned int nNN_, unsigned int nLoop_>
struct ct_par
{
	static const unsigned int nNN = nNN_;
	static const unsigned int nLoop = nLoop_;
};


template<typename Tsrc,typename Tdst>
class copy_prop_to_vector_block
{
	//! source
	Tsrc src;

	//! destination
	Tdst dst;

	size_t pos;

	unsigned int bPos;

public:

	copy_prop_to_vector_block(Tsrc src, Tdst dst,size_t pos, size_t bPos)
	:src(src),dst(dst),pos(pos),bPos(bPos)
	{}

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename std::remove_reference<decltype(dst.template get<T::value>())>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(src.template get<T::value>()[bPos],dst.template get<T::value>());
	}

};

template<unsigned int dim,
		 typename AggregateT,
		 unsigned int blockEdgeSize = default_edge<dim>::type::value,
		 unsigned int threadBlockSize = 128,
		 typename indexT=long int,
		 template<typename> class layout_base=memory_traits_inte,
		 typename linearizer = grid_smb<dim, blockEdgeSize>>
class SparseGridGpu : public BlockMapGpu<
        typename aggregate_convert<dim,blockEdgeSize,AggregateT>::type,
        threadBlockSize, indexT, layout_base>
{
private:
    const static unsigned char PADDING_BIT = 1;
	typedef typename aggregate_convert<dim,blockEdgeSize,AggregateT>::type AggregateBlockT;
    linearizer gridGeometry;
    grid_sm<dim, int> extendedBlockGeometry;
    grid_sm<dim, int> gridSize;
    unsigned int stencilSupportRadius;
    unsigned int ghostLayerSize;

    //! For stencil in a block-wise computation we have to load blocks + ghosts area. The ghost area live in neighborhood blocks
    //! For example the left ghost margin live in the right part of the left located neighborhood block, the right margin live in the
    //! left part of the of the right located neighborhood block, the top ...
    //! The first index indicate the index of the point in the block + ghost area, the second index indicate the correspondent neighborhood
    //! index (in a star like 0 mean negative x 1 positive x, 1 mean negative y and so on)
    openfpm::vector_gpu<aggregate<short int,short int>> ghostLayerToThreadsMapping;

    openfpm::vector_gpu<aggregate<indexT>> nn_blocks;

protected:
    static constexpr unsigned int blockSize = BlockTypeOf<AggregateBlockT, 0>::size;
    typedef AggregateBlockT AggregateInternalT;

public:
    static constexpr unsigned int dims = dim;
    static constexpr unsigned int blockEdgeSize_ = blockEdgeSize;

    typedef linearizer grid_info;

    template<typename Tfunc> using layout_mfunc = memory_traits_inte<Tfunc>;

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

    template<typename ... v_reduce>
    void flush(mgpu::ofp_context_t &context, flush_type opt = FLUSH_ON_HOST)
    {
        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
                ::template flush<v_reduce ...>(context, opt);
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
    	ghostLayerToThreadsMapping.resize(ghostLayerSize);
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
            mem_id linId = enlargedGrid.LinId(coord);
            // Mapping
            ghostLayerToThreadsMapping.template get<gt>(i) = linId;
            // Now compute the neighbour position to use
//            grid_key_dx<dim, int> localCoord = extendedBlockGeometry.InvLinId(linId);
            int neighbourNum = -1;
            int ctr = 0;
            for (int j = 0; j < dim; ++j)
            {
                int c = static_cast<int>(coord.get(j)) - static_cast<int>(stencilSupportRadius);
                if (c < 0)
                {
                    neighbourNum = 2*j;
                    ++ctr;
                }
                else if (c >= blockEdgeSize)
                {
                    neighbourNum = 2*j + 1;
                    ++ctr;
                }
            }
            if (ctr > 1) // If we are on a "corner"
            {
                neighbourNum = 0;
            }
            ghostLayerToThreadsMapping.template get<nt>(i) = neighbourNum;
            //
            ++i;
            ++gsi;
        }
        assert(i == ghostLayerSize);

        ghostLayerToThreadsMapping.template hostToDevice<gt,nt>();
    }

    void initialize(const size_t (& res)[dim])
    {
    	gridGeometry = linearizer(res);

        computeSizeOfGhostLayer();
        allocateGhostLayerMapping();
        computeGhostLayerMapping();

        size_t extBlockDims[dim];
        for (int d=0; d<dim; ++d)
        {
            extBlockDims[d] = blockEdgeSize + 2*stencilSupportRadius;
        }
        extendedBlockGeometry.setDimensions(extBlockDims);
        gridSize.setDimensions(res);
    }

    // This is for testing stencil application in Host, remember to deviceToHost the sparse grid before doing this!
    template <typename stencil, typename... Args>
    void applyStencilInPlaceHost(Args... args)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;

        SparseGridGpuKernels::applyStencilInPlaceHost
                <dim,
                        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
                        stencil>(
                         indexBuffer,
                         dataBuffer,
                         *this,
                         args...);
    }

    template <typename stencil, typename... Args>
    void applyStencilInPlace(Args... args)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;

        // NOTE: Here we want to work only on one data chunk per block!
        constexpr unsigned int chunksPerBlock = 1;
        const unsigned int localThreadBlockSize = dataChunkSize * chunksPerBlock;
        const unsigned int threadGridSize = numScalars % localThreadBlockSize == 0
                                            ? numScalars / localThreadBlockSize
                                            : 1 + numScalars / localThreadBlockSize;

        constexpr unsigned int nLoop = UIntDivCeil<(IntPow<blockEdgeSize + 2, dim>::value - IntPow<blockEdgeSize, dim>::value), (blockSize * chunksPerBlock)>::value; // todo: This works only for stencilSupportSize==1

        CUDA_LAUNCH_DIM3((SparseGridGpuKernels::applyStencilInPlace
                <dim,
                BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
                stencil>),
                threadGridSize, localThreadBlockSize,
                        indexBuffer.toKernel(),
                        dataBuffer.toKernel(),
                        this->template toKernelNN<stencil::stencil_type::nNN, nLoop>(),
                        args...);
    }

    //todo: the applyInsert should also allocate the gpu insert buffer, initialize it, etc
    template <typename stencil, typename... Args>
    void applyStencilInsert(Args... args)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;

        // NOTE: Here we want to work only on one data chunk per block!
        constexpr unsigned int chunksPerBlock = 1;
        const unsigned int localThreadBlockSize = dataChunkSize * chunksPerBlock;
        const unsigned int threadGridSize = numScalars % localThreadBlockSize == 0
                                      ? numScalars / localThreadBlockSize
                                      : 1 + numScalars / localThreadBlockSize;

        constexpr unsigned int nLoop = UIntDivCeil<(IntPow<blockEdgeSize + 2, dim>::value - IntPow<blockEdgeSize, dim>::value), (blockSize * chunksPerBlock)>::value; // todo: This works only for stencilSupportSize==1

        setGPUInsertBuffer(threadGridSize, chunksPerBlock);
//        setGPUInsertBuffer(threadGridSize, localThreadBlockSize);

        CUDA_LAUNCH_DIM3((SparseGridGpuKernels::applyStencilInsert
                <dim,
                BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
                stencil>),
                threadGridSize, localThreadBlockSize,
                        indexBuffer.toKernel(),
                        dataBuffer.toKernel(),
                        this->template toKernelNN<stencil::stencil_type::nNN, nLoop>(),
                        args...);

//        cudaDeviceSynchronize();
        mgpu::ofp_context_t ctx;
        stencil::flush(*this, ctx);
    }

    template <typename stencil, typename... Args>
    void applyBoundaryStencilInPlace(Args... args)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;

        // NOTE: Here we want to work only on one data chunk per block!
        constexpr unsigned int chunksPerBlock = 1;
        const unsigned int localThreadBlockSize = dataChunkSize * chunksPerBlock;
        const unsigned int threadGridSize = numScalars % localThreadBlockSize == 0
                                            ? numScalars / localThreadBlockSize
                                            : 1 + numScalars / localThreadBlockSize;

        constexpr unsigned int nLoop = UIntDivCeil<(IntPow<blockEdgeSize + 2, dim>::value - IntPow<blockEdgeSize, dim>::value), (blockSize * chunksPerBlock)>::value; // todo: This works only for stencilSupportSize==1

        CUDA_LAUNCH_DIM3((SparseGridGpuKernels::applyBoundaryStencilInPlace
                <dim,
                        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
                        stencil>),
                         threadGridSize, localThreadBlockSize,
                         indexBuffer.toKernel(),
                         dataBuffer.toKernel(),
                         this->template toKernelNN<stencil::stencil_type::nNN, nLoop>(),
                         args...);
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
    SparseGridGpu(linearizer & gridGeometry, unsigned int stencilSupportRadius = 1)
            : gridGeometry(gridGeometry),
              stencilSupportRadius(stencilSupportRadius)
    {
    	initialize(gridGeometry.getSize());
    };

    SparseGridGpu_ker
            <
                    dim,
                    blockEdgeSize,
                    typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                    ct_par<0,1>,
                    indexT,
                    layout_base,
                    decltype(extendedBlockGeometry),
                    linearizer
            > toKernel()
    {
        SparseGridGpu_ker
                <
                        dim,
                        blockEdgeSize,
                        typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                        ct_par<0,1>,
                        indexT,
                        layout_base,
                        decltype(extendedBlockGeometry),
                        linearizer
                > toKer(
                BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.toKernel(),
                gridGeometry,
                extendedBlockGeometry,
                stencilSupportRadius,
                ghostLayerToThreadsMapping.toKernel(),
                nn_blocks.toKernel(),
                ghostLayerSize);
        return toKer;
    }

    template<unsigned int nNN, unsigned int nLoop>
    SparseGridGpu_ker
            <
                    dim,
                    blockEdgeSize,
                    typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                    ct_par<nNN,nLoop>,
                    indexT,
                    layout_base,
                    decltype(extendedBlockGeometry),
                    linearizer
            > toKernelNN()
    {
        SparseGridGpu_ker
                <
                        dim,
                        blockEdgeSize,
                        typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                        ct_par<nNN,nLoop>,
                        indexT,
                        layout_base,
                        decltype(extendedBlockGeometry),
                        linearizer
                > toKer(
                BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.toKernel(),
                gridGeometry,
                extendedBlockGeometry,
                stencilSupportRadius,
                ghostLayerToThreadsMapping.toKernel(),
                nn_blocks.toKernel(),
                ghostLayerSize);
        return toKer;
    }

    //

    constexpr static unsigned int getBlockEdgeSize()
    {
        return blockEdgeSize;
    }

    constexpr unsigned int getBlockSize() const
    {
        return blockSize;
    }

    // Geometry
    template<typename CoordT>
    inline size_t getLinId(CoordT &coord)
    {
        return gridGeometry.LinId(coord);
    }

    inline grid_key_dx<dim, int> getCoord(size_t linId)
    {
        return gridGeometry.InvLinId(linId);
    }

    inline ite_gpu<dim> getGridGPUIterator(const grid_key_dx<dim, int> & start, const grid_key_dx<dim, int> & stop, size_t n_thr = threadBlockSize)
    {
    	return gridSize.getGPUIterator(start,stop,n_thr);
    }

    // Data management methods

    template<unsigned int p, typename CoordT>
    auto get(const CoordT &coord) const -> const ScalarTypeOf<AggregateBlockT, p> &
    {
        return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::template get<p>(gridGeometry.LinId(coord));
    }

    template<unsigned int p, typename CoordT>
    auto insert(const CoordT &coord) -> ScalarTypeOf<AggregateBlockT, p> &
    {
        return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::template insert<p>(gridGeometry.LinId(coord));
    }

    template<typename dim3T>
    void setGPUInsertBuffer(dim3T nBlock, dim3T nSlot)
    {
        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
        ::setGPUInsertBuffer(
                dim3SizeToInt(nBlock),
                dim3SizeToInt(nSlot)
        );
    }

    template<typename stencil_type = NNStar>
    void tagBoundaries()
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;

        // NOTE: Here we want to work only on one data chunk per block!

        unsigned int localThreadBlockSize = dataChunkSize;
        unsigned int threadGridSize = numScalars % dataChunkSize == 0
                                    ? numScalars / dataChunkSize
                                    : 1 + numScalars / dataChunkSize;

        constexpr unsigned int nLoop = UIntDivCeil<(IntPow<blockEdgeSize + 2, dim>::value - IntPow<blockEdgeSize, dim>::value), (blockSize * 1)>::value; // todo: This works only for stencilSupportSize==1
//        constexpr unsigned int nLoop = IntPow<blockEdgeSize + 2, dim>::value; // todo: This works only for stencilSupportSize==1

        if (stencilSupportRadius == 1)
        {
            CUDA_LAUNCH_DIM3((SparseGridGpuKernels::tagBoundaries<
                    dim,
                    1,
                    BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>),
                    threadGridSize, localThreadBlockSize,indexBuffer.toKernel(), dataBuffer.toKernel(), this->template toKernelNN<stencil_type::nNN, nLoop>(), nn_blocks.toKernel());
        }
        else if (stencilSupportRadius == 2)
        {
        	CUDA_LAUNCH_DIM3((SparseGridGpuKernels::tagBoundaries<
                    dim,
                    2,
                    BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>),
                    threadGridSize, localThreadBlockSize,indexBuffer.toKernel(), dataBuffer.toKernel(), this->template toKernelNN<stencil_type::nNN, nLoop>(), nn_blocks.toKernel());
        }
        else if (stencilSupportRadius == 0)
        {
        	CUDA_LAUNCH_DIM3((SparseGridGpuKernels::tagBoundaries<
                    dim,
                    0,
                    BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>),
                    threadGridSize, localThreadBlockSize,indexBuffer.toKernel(), dataBuffer.toKernel(), this->template toKernelNN<stencil_type::nNN, nLoop>(), nn_blocks.toKernel());
        }
        else
        {
            //todo: KABOOOOOOM!!!
            std::cout << __FILE__ << ":" << __LINE__ << " error: stencilSupportRadius supported only up to 2, passed: " << stencilSupportRadius << std::endl;

        }
        cudaDeviceSynchronize();
    }

    template<typename NNtype = NNStar>
    void findNeighbours()
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();

        const unsigned int numBlocks = indexBuffer.size();
        const unsigned int numScalars = numBlocks * NNtype::nNN;
        nn_blocks.resize(numScalars);

        if (numScalars == 0) return;

        // NOTE: Here we want to work only on one data chunk per block!

        unsigned int localThreadBlockSize = NNtype::nNN;

        unsigned int threadGridSize = numScalars % localThreadBlockSize == 0
                                      ? numScalars / localThreadBlockSize
                                      : 1 + numScalars / localThreadBlockSize;
        CUDA_LAUNCH_DIM3((SparseGridGpuKernels::findNeighbours<dim,NNtype::nNN>),
                threadGridSize, localThreadBlockSize,indexBuffer.toKernel(), this->toKernel(),nn_blocks.toKernel());
        cudaDeviceSynchronize();
    }

    size_t countExistingElements()
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        constexpr unsigned int pMask = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask;
        typedef typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT BAggregateT;
        typedef BlockTypeOf<BAggregateT, pMask> MaskBlockT;
        constexpr unsigned int blockSize = MaskBlockT::size;
        const auto bufferSize = indexBuffer.size();

        size_t numExistingElements = 0;

        for (size_t blockId=0; blockId<bufferSize; ++blockId)
        {
            auto dataBlock = dataBuffer.get(blockId); // Avoid binary searches as much as possible
            for (size_t elementId=0; elementId<blockSize; ++elementId)
            {
                const auto curMask = dataBlock.template get<pMask>()[elementId];

                if (this->exist(curMask))
                {
                    ++numExistingElements;
                }
            }
        }

        return numExistingElements;
    }

    size_t countBoundaryElements()
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        constexpr unsigned int pMask = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask;
        typedef typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT BAggregateT;
        typedef BlockTypeOf<BAggregateT, pMask> MaskBlockT;
        constexpr unsigned int blockSize = MaskBlockT::size;
        const auto bufferSize = indexBuffer.size();

        size_t numBoundaryElements = 0;

        for (size_t blockId=0; blockId<bufferSize; ++blockId)
        {
            auto dataBlock = dataBuffer.get(blockId); // Avoid binary searches as much as possible
            for (size_t elementId=0; elementId<blockSize; ++elementId)
            {
                const auto curMask = dataBlock.template get<pMask>()[elementId];

                if (this->exist(curMask) && this->isPadding(curMask))
                {
                    ++numBoundaryElements;
                }
            }
        }

        return numBoundaryElements;
    }

    // Also count mean+stdDev of occupancy of existing blocks
    void measureBlockOccupancy(double &mean, double &deviation)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        constexpr unsigned int pMask = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask;
        typedef typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT BAggregateT;
        typedef BlockTypeOf<BAggregateT, pMask> MaskBlockT;
        constexpr unsigned int blockSize = MaskBlockT::size;
        const auto bufferSize = indexBuffer.size();

        openfpm::vector<double> measures;

        for (size_t blockId=0; blockId<bufferSize; ++blockId)
        {
            auto dataBlock = dataBuffer.get(blockId); // Avoid binary searches as much as possible
            size_t numElementsInBlock = 0;
            for (size_t elementId=0; elementId<blockSize; ++elementId)
            {
                const auto curMask = dataBlock.template get<pMask>()[elementId];

                if (this->exist(curMask))
                {
                    ++numElementsInBlock;
                }
            }
            double blockOccupancy = static_cast<double>(numElementsInBlock)/blockSize;
            measures.add(blockOccupancy);
        }

        standard_deviation(measures, mean, deviation);
    }

    template<typename stencil, typename... Args>
    void applyStencilsHost(StencilMode mode = STENCIL_MODE_INPLACE, Args... args)
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
                applyStencilInPlaceHost<stencil>(args...);
                break;
            case STENCIL_MODE_INSERT:
//                applyStencilInsertHost<stencil>(args...);
                std::cout << "No insert mode stencil application available on Host!" << std::endl;
                break;
        }
    }
    template<typename stencil1, typename stencil2, typename ... otherStencils, typename... Args>
    void applyStencilsHost(StencilMode mode = STENCIL_MODE_INSERT, Args... args)
    {
        applyStencilsHost<stencil1>(mode, args...);
        applyStencilsHost<stencil2, otherStencils ...>(mode, args...);
    }

    //todo: Move implems into a functor for compile time choice of stencil mode
    template<typename stencil, typename... Args>
    void applyStencils(StencilMode mode = STENCIL_MODE_INSERT, Args... args)
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
    template<typename stencil1, typename stencil2, typename ... otherStencils, typename... Args>
    void applyStencils(StencilMode mode = STENCIL_MODE_INSERT, Args... args)
    {
        applyStencils<stencil1>(mode, args...);
        applyStencils<stencil2, otherStencils ...>(mode, args...);
    }

    //todo: Move implems into a functor for compile time choice of stencil mode
    template<typename stencil, typename... Args>
    void applyBoundaryStencils(Args... args)
    {
        // Apply the given stencil on all elements which are not boundary-tagged
        // The idea is to have this function launch a __global__ function (made by us) on all existing blocks
        // then this kernel checks if elements exist && !padding and on those it calls the user-provided
        // __device__ functor. The mode of the stencil application is used basically to choose how we load the block
        // that we pass to the user functor as storeBlock: in case of Insert, we get the block through an insert (and
        // we also call the necessary aux functions); in case of an In-place we just get the block from the data buffer.
        StencilMode mode = STENCIL_MODE_INPLACE;
        switch (mode)
        {
            case STENCIL_MODE_INPLACE:
                applyBoundaryStencilInPlace<stencil>(args...);
                break;
            case STENCIL_MODE_INSERT:
//                applyBoundaryStencilInsert<stencil>(args...);
                std::cout << "No insert mode stencil application available on BOUNDARY!" << std::endl;
                break;
        }
    }
    template<typename stencil1, typename stencil2, typename ... otherStencils, typename... Args>
    void applyBoundaryStencils(Args... args)
    {
        applyBoundaryStencils<stencil1>(args...);
        applyBoundaryStencils<stencil2, otherStencils ...>(args...);
    }

    template<typename BitMaskT>
    inline static bool isPadding(BitMaskT &bitMask)
    {
        return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
        ::getBit(bitMask, PADDING_BIT);
    }

    template<typename BitMaskT>
    inline static void setPadding(BitMaskT &bitMask)
    {
        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
        ::setBit(bitMask, PADDING_BIT);
    }

    template<typename BitMaskT>
    inline static void unsetPadding(BitMaskT &bitMask)
    {
        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
        ::unsetBit(bitMask, PADDING_BIT);
    }

    template<unsigned int p>
    void print_vct_add_data()
    {
    	typedef BlockMapGpu<
    	        typename aggregate_convert<dim,blockEdgeSize,AggregateT>::type,
    	        threadBlockSize, indexT, layout_base> BMG;

    	auto & bM = BMG::blockMap.private_get_vct_add_data();
    	auto & vI = BMG::blockMap.private_get_vct_add_index();
    	bM.template deviceToHost<p>();
    	vI.template deviceToHost<0>();

    	std::cout << "vct_add_data: " << std::endl;

    	for (size_t i = 0 ; i < bM.size() ; i++)
    	{
    		std::cout << i << "  index: " << vI.template get<0>(i) << "  BlockData: " << std::endl;
    		for (size_t j = 0 ; j < blockSize ; j++)
    		{
    			std::cout << (int)bM.template get<p>(i)[j] << "  ";
    		}

    		std::cout << std::endl;
    	}
    }

#ifdef OPENFPM_DATA_ENABLE_IO_MODULE

	/*! \brief write the sparse grid into VTK
	 *
	 * \param out VTK output
	 *
	 */
	template<typename Tw = float> bool write(const std::string & output)
	{
		file_type ft = file_type::BINARY;

		auto & bm = this->private_get_blockMap();

		auto & index = bm.getIndexBuffer();
		auto & data = bm.getDataBuffer();

		openfpm::vector<Point<dim,Tw>> tmp_pos;
		openfpm::vector<typename aggregate_add<AggregateT>::type> tmp_prp;

		// copy position and properties

		auto it = index.getIterator();

		while(it.isNext())
		{
			auto key = it.get();

			Point<dim,Tw> p;

			for (size_t i = 0 ; i < gridGeometry.getBlockSize() ; i++)
			{
				if (data.template get<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>(key)[i] != 0)
				{
					// Get the block index
					grid_key_dx<dim,int> keyg = gridGeometry.InvLinId(index.template get<0>(key),i);

					for (size_t k = 0 ; k < dim ; k++)
					{p.get(k) = keyg.get(k);}

					tmp_pos.add(p);

					tmp_prp.add();
					copy_prop_to_vector_block<decltype(data.get_o(key)),decltype(tmp_prp.last())>
					cp(data.get_o(key),tmp_prp.last(),key,i);

					boost::mpl::for_each_ref< boost::mpl::range_c<int,0,AggregateT::max_prop> >(cp);

					tmp_prp.last().template get<AggregateT::max_prop>() = data.template get<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>(key)[i];
				}
			}

			++it;
		}

		// VTKWriter for a set of points
		VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim,Tw>>, openfpm::vector<typename aggregate_add<AggregateT>::type>>, VECTOR_POINTS> vtk_writer;
		vtk_writer.add(tmp_pos,tmp_prp,tmp_pos.size());

		openfpm::vector<std::string> prp_names;

		// Write the VTK file
		return vtk_writer.write(output,prp_names,"sparse_grid","",ft);
	}

#endif
};

template<unsigned int dim,
		 typename AggregateT,
		 unsigned int blockEdgeSize = default_edge<dim>::type::value,
		 unsigned int threadBlockSize = 128,
		 typename indexT=long int,
		 template<typename> class layout_base=memory_traits_inte,
		 typename linearizer = grid_zmb<dim, blockEdgeSize>>
using SparseGridGpu_z = SparseGridGpu<dim,AggregateT,blockEdgeSize,threadBlockSize,indexT,layout_base,linearizer>;

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_HPP
