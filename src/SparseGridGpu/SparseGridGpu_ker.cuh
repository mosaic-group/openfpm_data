//
// Created by tommaso on 11/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_KER_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_KER_CUH

#include <SparseGridGpu/Geometry/grid_smb.hpp>
#include "BlockMapGpu.hpp"
#include "SparseGridGpu_ker_util.hpp"

template<typename indexT>
struct block_offset
{
	indexT pos;
	indexT off;
};

//todo Remove template param GridSmT and just use BlockGeometry
template<unsigned int dim,
        unsigned int blockEdgeSize,
        typename AggregateBlockT,
        typename ct_params,
        typename indexT,
        template<typename> class layout_base,
        typename GridSmT,
        typename linearizer>
class SparseGridGpu_ker : public BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>
{
private:
    linearizer grid;
    GridSmT blockWithGhostGrid;

protected:
    const static unsigned char PADDING_BIT = 1;
    static constexpr unsigned int blockSize = BlockTypeOf<AggregateBlockT, 0>::size;
    unsigned int ghostLayerSize;
    openfpm::vector_gpu_ker<aggregate<short int,short int>,memory_traits_inte> ghostLayerToThreadsMapping;
    openfpm::vector_gpu_ker<aggregate<indexT>,memory_traits_inte> nn_blocks;
    openfpm::vector_gpu_ker<aggregate<indexT>,memory_traits_inte> buffPnt;

public:
    static constexpr unsigned int d = dim;
    static constexpr unsigned int blockEdgeSize_ = blockEdgeSize;
    unsigned int stencilSupportRadius;
    typedef AggregateBlockT AggregateBlockType;
    typedef indexT indexT_;

	//! Indicate this structure has a function to check the device pointer
	typedef int yes_has_check_device_pointer;

public:

	/*! \brief constructor
	 *
	 * this constructor is in general called in the function toKernel() of SparseGridGpu
	 *
	 */
    SparseGridGpu_ker(const openfpm::vector_sparse_gpu_ker<AggregateBlockT, indexT, layout_base> &blockMap,
                      linearizer & grid,
                      GridSmT extendedBlockGeometry,
                      unsigned int stencilSupportRadius,
                      openfpm::vector_gpu_ker<aggregate<short int,short int>,memory_traits_inte> ghostLayerToThreadsMapping,
                      openfpm::vector_gpu_ker<aggregate<indexT>,memory_traits_inte> nn_blocks,
                      openfpm::vector_gpu_ker<aggregate<indexT>,memory_traits_inte> buffPnt,
                      unsigned int ghostLayerSize)
            : BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>(blockMap),
              grid(grid),
              blockWithGhostGrid(extendedBlockGeometry),
              stencilSupportRadius(stencilSupportRadius),
              ghostLayerSize(ghostLayerSize),
              ghostLayerToThreadsMapping(ghostLayerToThreadsMapping),
              nn_blocks(nn_blocks),
              buffPnt(buffPnt)
    {}

    /*! \brief Get the coordinate of the block and the offset id inside the block it give the global coordinate
     *
     * \param blockCoord block coordinate
     * \param offset in the block
     *
     * \return the global coordinate position
     *
     */
    template<typename CoordT>
    __device__ __host__ inline grid_key_dx<dim,CoordT> getGlobalCoord(const grid_key_dx<dim, CoordT> & blockCoord, unsigned int offset)
    {
    	return grid.getGlobalCoord(blockCoord,offset);
    }

    /*! \brief Linearization of  global coordinates
     *
     * \param coord point coordinates
     *
     * \return the linearized position
     *
     */
    template<typename CoordT>
    inline __device__ size_t getLinId(const grid_key_dx<dim, CoordT> & coord) const
    {
        return grid.LinId(coord);
    }

    /*! \brief Size of the sparse grid  in each direction
     *
     * \param i direction
     *
     * \return the linearized position
     *
     */
    inline __device__ unsigned int size(unsigned int i)
    {
    	return grid.getSize()[i];
    }

    /*! \brief The inversion of getLinId
     *
     * \param linearized index
     *
     * \return the point coordinate
     *
     */
    inline __device__ grid_key_dx<dim, int> getCoord(size_t linId) const
    {
        return grid.InvLinId(linId);

    }

    /*! \brief Given a point to insert, return the block-id and offset of that point
     *
     * \param p Point the thread is processing
     * \param blk Block id the point is falling into
     * \param offset linearized index within the block
     *
     * \return true if the point is inactive
     *
     */
    template<typename ite_type>
    inline __device__ bool getInsertBlockOffset(const ite_type & itd, const grid_key_dx<dim, int> & p, grid_key_dx<dim, int> & blk, int & offset)
    {
    	int accu = 1;
    	offset = 0;

    	bool active = true;

    	for (int i = 0 ; i < dim ; i++)
    	{
    		blk.set_d(i,p.get(i) / getBlockEdgeSize());
    		offset += (p.get(i) % getBlockEdgeSize()) * accu;
    		accu *= getBlockEdgeSize();
    		active = active && (p.get(i) >= (itd.start.get(i) + itd.start_base.get(i))) && (p.get(i) <= itd.stop.get(i));
    	}

    	return active;
    }

    /*! \brief Linearization of  block coordinates
     *
     * \param blockCoord block coordinates
     *
     * \return the linearized index
     *
     */
    template<typename CoordT>
    inline __device__ size_t getBlockLinId(CoordT blockCoord) const
    {
        return grid.BlockLinId(blockCoord);
    }

    /*! \brief The inversion of getBlockLinId
     *
     * \param blockLinId linearized block index
     *
     * \return the block coordinates
     *
     */
    inline __device__ grid_key_dx<dim, int> getBlockCoord(size_t blockLinId) const
    {
        return grid.BlockInvLinId(blockLinId);
    }

    /*! \brief Given a linearized block index it return the coordinated of the lower-left point in 2D or
     *         in general the origin of the block in global coordinates
     *
     * \param blockLinId linearized block index
     *
     * \return the block coordinates
     *
     */
    inline __device__ grid_key_dx<dim, int> getBlockBaseCoord(size_t blockLinId) const
    {
        return grid.InvLinId(blockLinId * blockSize);
    }

    inline __device__ grid_key_dx<dim, int> getNeighbour(grid_key_dx<dim, int> base, unsigned int dimension, char offset) const
    {
        grid_key_dx<dim, int> res = base;
        auto i = base.get(dimension) + offset;
        res.set_d(dimension, i);
        return res;
    }

    /*! \brief Return the size of the block edge size
     *
     * \return the block edge size
     *
     */
    constexpr static __device__ unsigned int getBlockEdgeSize()
    {
        return blockEdgeSize;
    }


    /*! \brief Return the size of the block
     *
     * \return the block size
     *
     */
    constexpr __device__ unsigned int getBlockSize() const
    {
        return blockSize;
    }

    /*! \brief Return the size of the block + ghost needed to apply the stencil
     *
     * \return the block size + ghost
     *
     */
    inline __device__ unsigned int getEnlargedBlockSize() const
    {
        return std::pow(blockEdgeSize + 2*stencilSupportRadius, dim);
    }

    /*! \brief Get the neighborhood point in one direction
     *
     * \warning This function assume you do not move more than one block edge size in one direction
     *
     * \param pos data block position
     * \param offset inside the data block
     * \param position where to move
     *
     */
    template<typename NN_type, typename indexT2>
    inline __device__  block_offset<indexT2> getNNPoint(openfpm::sparse_index<unsigned int> pos,
    										unsigned int offset,
    										const grid_key_dx<dim,indexT2> & mov)
    {
    	block_offset<indexT2> bof;

    	grid_key_dx<dim,indexT2> coord;

    	for (int i = 0 ; i < dim ; i++)
    	{
    		coord.set_d(i,mov.get(i) + offset % blockEdgeSize);
    		offset /= blockEdgeSize;
    	}

    	unsigned int NN_index = 0;
    	unsigned int offset_nn = 0;

    	bool out = NN_type::template getNNindex_offset<blockEdgeSize>(coord,NN_index,offset_nn);

    	// Calculate internal coordinates

    	indexT nnb = pos.id;

    	if (out == true)
    	{nnb = nn_blocks.template get<0>(NN_index + NN_type::nNN*pos.id);}

    	bof.pos = nnb;
    	bof.off = offset_nn;

    	return bof;
    }

    inline __device__ unsigned int posToEnlargedBlockPos(unsigned int pos) const
    {
        // Convert pos into a linear id accounting for the ghost offsets
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(pos, stencilSupportRadius, coord);
        const unsigned int linId = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
        return linId;
    }

    inline __device__ grid_key_dx<dim,int>
    getCoordInEnlargedBlock(const unsigned int offset) const
    {
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(offset, stencilSupportRadius, coord);
        return grid_key_dx<dim, int>(coord);
    }

    inline __device__ unsigned int
    getLinIdInEnlargedBlock(const unsigned int offset) const
    {
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(offset, stencilSupportRadius, coord);
        return coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
    }

    template<typename Coordtype>
    inline __device__ unsigned int
    getNeighbourLinIdInEnlargedBlock(const grid_key_dx<dim, Coordtype> & base, grid_key_dx<dim, Coordtype> & offsets) const
    {
        grid_key_dx<dim, int> res = base + offsets;
        return coordToLin<blockEdgeSize>(res, stencilSupportRadius);
    }

    template<typename Coordtype>
    inline __device__ unsigned int
    getNeighbourLinIdInEnlargedBlock(const grid_key_dx<dim,Coordtype> & base, unsigned int dimension, char offset) const
    {
        grid_key_dx<dim, int> res = getNeighbour(base, dimension, offset);
        return coordToLin<blockEdgeSize>(res, stencilSupportRadius);
    }

    inline __device__ bool
    getIfBoundaryElementInEnlargedBlock(const grid_key_dx<dim, int> coordInEnlargedBlock, char (&boundaryDirection)[dim])
    {
        bool isBoundary = false;
        for (int d=0; d<dim; ++d)
        {
            const auto v = coordInEnlargedBlock.get(d);
            if (v==stencilSupportRadius)
            {
                boundaryDirection[d] = -1;
                isBoundary = true;
            }
            else if (v==stencilSupportRadius+blockEdgeSize-1)
            {
                boundaryDirection[d] = 1;
                isBoundary = true;
            }
            else
            {
                boundaryDirection[d] = 0;
            }
        }
        return isBoundary;
    }

    // Data management methods

    template<unsigned int p, typename CoordT>
    inline __device__ auto
    get(const grid_key_dx<dim, CoordT> & coord) const -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template get<p>(
            0))
    {
    	return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template get<p>(grid.LinId(coord));
    }

    // Data management methods

    template<typename CoordT>
    inline __device__ void
    get_sparse(const grid_key_dx<dim, CoordT> & coord, unsigned int & dataBlockPos, unsigned int & offset) const
    {
    	return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::get_sparse(grid.LinId(coord),dataBlockPos,offset);
    }

    /*! \brief Access the grid point
     *
     * \param coord point
     *
     * \return a reference to the data point
     *
     */
    template<unsigned int p, typename CoordT>
    inline __device__ auto
    get(const block_offset<CoordT> & coord) const -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::blockMap.template get_ele<p>(coord.pos)[coord.off])
    {
    	return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::blockMap.template get_ele<p>(coord.pos)[coord.off];
    }

    /*! \brief Access the grid point
     *
     * \param coord point
     *
     * \return a reference to the data point
     *
     */
    template<unsigned int p, typename CoordT>
    inline __device__ auto
    get(const block_offset<CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::blockMap.template get_ele<p>(coord.pos)[coord.off])
    {
    	return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::blockMap.template get_ele<p>(coord.pos)[coord.off];
    }

    template<unsigned int p, typename CoordT>
    inline __device__ auto
    insert(const grid_key_dx<dim, CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template insert<p>(0))
    {
    	return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::template insert<p>(grid.LinId(coord));
    }

    template<typename CoordT>
    inline __device__ unsigned int getBlockId(const grid_key_dx<dim, CoordT> & coord)
    {
        // todo: check this because it's bugged! maybe?
        return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlockId(grid.LinId(coord));
    }

    template<typename CoordT>
    inline __device__ unsigned int getOffset(const grid_key_dx<dim, CoordT> & coord)
    {
        return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getOffset(grid.LinId(coord));
    }

    template<typename CoordT>
    inline __device__ auto
    getBlock(const grid_key_dx<dim, CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlock(0))
    {
        return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlock(getBlockId(coord));
    }

    inline __device__ auto
    getBlock(const unsigned int blockLinId) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlock(0))
    {
        return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::getBlock(blockLinId);
    }

    template<unsigned int chunksPerBlocks = 1,typename CoordT>
    inline __device__ auto
    insertBlock(const grid_key_dx<dim, CoordT> & coord) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock(0))
    {
        return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock(getBlockId(coord));
    }

    template<unsigned int chunksPerBlocks = 1>
    inline __device__ auto
    insertBlock(const indexT blockLinId, const unsigned int stride = 8192) -> decltype(BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock(0))
    {
        return BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>::insertBlock<chunksPerBlocks>(blockLinId,stride);
    }

    /*! \brief Return the buffer of points
     *
     * \return the buffer of points
     *
     */
    inline __device__ auto getPointBuffer() -> decltype(buffPnt) &
    {
    	return buffPnt;
    }

    // Load & Store aux functions for user kernels. To be used for loading to or writing from shared memory.

    /**
     * Load a data block into the inner part of a shared memory region.
     * The given shared memory region should be shaped as a dim-dimensional array and sized so that it
     * can contain the block plus the ghost layer around it.
     *
     * \tparam p The property to retrieve from global memory.
     * \tparam CoordT The coordinate type.
     * \tparam Shared The type of the shared memory region.
     * \param coord The coordinate of the block.
     * \param sharedRegionPtr The pointer to the shared memory region.
     */
    template<unsigned int p, typename AggrWrapperT>
    inline __device__ void
    loadBlock(AggrWrapperT &block, ScalarTypeOf<AggregateBlockT, p> *sharedRegion)
    {
        //todo: Make this work well with multiples chunks per block or check not to get several chunks or dragons ahoy!
        __loadBlock<p>(block, sharedRegion);
    }

    /**
     * Load a data block into the inner part of a shared memory region.
     * The given shared memory region should be shaped as a dim-dimensional array and sized so that it
     * can contain the block plus the ghost layer around it.
     *
     * \tparam CoordT The coordinate type.
     * \tparam props The set of properties to retrieve from global memory.
     * \param coord The coordinate of the block.
     * \param sharedRegionPtr The array of pointers to the shared memory regions, one for each property.
     */
    template<unsigned int ... props, typename AggrWrapperT>
    inline __device__ void loadBlock(AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)])
    {
        __loadBlock<props ...>(block, sharedRegionPtr);
    }


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
    template<unsigned int p , typename AggrWrapperT , typename CoordT>
    inline __device__ void
    loadGhostBlock(const AggrWrapperT & dataBlockLoad,const grid_key_dx<dim, CoordT> & coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion)
    {
        auto blockLinId = getBlockId(coord);
        __loadGhostBlock<p>(dataBlockLoad,blockLinId, sharedRegion);
    }

    template<unsigned int p, typename AggrWrapperT>
    inline __device__ void
    loadGhostBlock(const AggrWrapperT & dataBlockLoad, const openfpm::sparse_index<unsigned int> blockLinId, ScalarTypeOf<AggregateBlockT, p> *sharedRegion)
    {
        __loadGhostBlock<p>(dataBlockLoad,blockLinId, sharedRegion);
    }

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
    inline __device__ void loadGhost(const grid_key_dx<dim, CoordT> & coord, const int * neighboursPos, void *sharedRegionPtr[sizeof...(props)])
    {
        auto blockLinId = getBlockId(coord);
        __loadGhost<props ...>(blockLinId, neighboursPos, sharedRegionPtr);
    }


    template<unsigned int ... props>
    inline __device__ void loadGhost(const unsigned int blockLinId, const int * neighboursPos, void *sharedRegionPtr[sizeof...(props)])
    {
        __loadGhost<props ...>(blockLinId, neighboursPos, sharedRegionPtr);
    }

    template<unsigned int p, typename AggrWrapperT>
    inline __device__ void
    storeBlock(AggrWrapperT &block, ScalarTypeOf<AggregateBlockT, p> *sharedRegion)
    {
        //todo: Make this work well with multiples chunks per block or check not to get several chunks or dragons ahoy!
        __storeBlock<p>(block, sharedRegion);
    }

    template<unsigned int p, typename CoordT>
    inline __device__ void
    storeBlockInPlace(const grid_key_dx<dim, CoordT> & coord, ScalarTypeOf<AggregateBlockT, p> *sharedRegion)
    {
        //todo: Make this work well with multiples chunks per block or check not to get several chunks or dragons ahoy!
        auto & block = getBlock(coord);
        __storeBlock<p>(block, sharedRegion);
    }

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
    template<unsigned int ... props, typename AggrWrapperT>
    inline __device__ void storeBlock(AggrWrapperT &block, void *sharedRegionPtr[sizeof...(props)])
    {
        __storeBlock<props ...>(block, sharedRegionPtr);
    }

    template<unsigned int ... props, typename CoordT>
    inline __device__ void storeBlockInPlace(const grid_key_dx<dim, CoordT> & coord, void *sharedRegionPtr[sizeof...(props)])
    {
        auto block = getBlock(coord);
        __storeBlock<props ...>(block, sharedRegionPtr);
    }

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

    template<typename NNtype>
    inline __device__ indexT getNeighboursPos(const indexT blockId, const unsigned int offset)
    {
        //todo: also do the full neighbourhood version, this is just cross
        auto blockCoord = getBlockCoord(blockId);

        return NNtype::template getNNpos<indexT>(blockCoord,this->blockMap,*this,offset);
    }

#ifdef SE_CLASS1

		/*! \brief Check if the device pointer is owned by this structure
		 *
		 * \return a structure pointer check with information about the match
		 *
		 */
		pointer_check check_device_pointer(void * ptr)
		{
			pointer_check pc;

			pc = ghostLayerToThreadsMapping.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("ghostLayerToThreadsMapping overflow : ") + "\n" + pc.match_str;
				return pc;
			}

			pc = nn_blocks.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("nn_blocks overflow: ") + "\n" + pc.match_str;
				return pc;
			}

			pc = ((BlockMapGpu_ker<AggregateBlockT, indexT, layout_base> *)this)->check_device_pointer(ptr);

			return pc;
		}

#endif

private:


    template<unsigned int p, typename AggrWrapperT, typename SharedPtrT>
    inline __device__ void
    __loadBlock(const AggrWrapperT &block, SharedPtrT sharedRegionPtr)
    {
        typedef ScalarTypeOf<AggregateBlockT, p> ScalarT;

        const unsigned int pos = threadIdx.x;
        //todo: Improve this version to allow multiple chunks per block!

        // Convert pos into a linear id accounting for the ghost offsets
        unsigned int coord[dim];
        linToCoordWithOffset<blockEdgeSize>(pos, stencilSupportRadius, coord);
        const unsigned int linId = coordToLin<blockEdgeSize>(coord, stencilSupportRadius);
        // Actually load the data into the shared region
        //ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
        sharedRegionPtr[linId] = block.template get<p>()[pos];
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
    __loadGhostNoNN(const unsigned int blockId, SharedPtrT * sharedRegionPtr)
    {
        typedef ScalarTypeOf<AggregateBlockT, p> ScalarT;

        const unsigned int edge = blockEdgeSize + 2*stencilSupportRadius;

        grid_key_dx<dim, int> localCoord;
        grid_key_dx<dim, int> elementCoord;

        for (int pos = threadIdx.x; pos < ghostLayerSize; pos += blockDim.x)
        {
                // Convert pos into a linear id accounting for the inner domain offsets
                const unsigned int linId = ghostLayerToThreadsMapping.template get<0>(pos);
                // Now get linear offset wrt the first element of the block
                elementCoord = getBlockBaseCoord(blockId);
                unsigned int ctr = linId;
                for (int i = 0; i < dim; ++i)
                {
                    int v = (ctr % edge) - stencilSupportRadius;
                    ctr /= edge;
                    elementCoord.set_d(i, elementCoord.get(i) + v);
                }

                // Actually load the data into the shared region
                ScalarT *basePtr = (ScalarT *)sharedRegionPtr;
                *(basePtr + linId) = get<p>(elementCoord);
        }
    }

    /*! \brief Load the ghost area in the shared region
     *
     * \param blockId Index of the center block
     * \param neighboursPos neighborhood blocks around blockId (if neighbourPos is null loadGhost)
     * \param sharedRegionPtr shared region
     *
     */
    template<unsigned int p, typename AggrWrapperT ,typename SharedPtrT>
    inline __device__ void
    __loadGhostBlock(const AggrWrapperT &block, const openfpm::sparse_index<unsigned int> blockId, SharedPtrT * sharedRegionPtr)
    {
    	loadGhostBlock_impl<ct_params::nLoop,dim,AggrWrapperT,p,ct_params,blockEdgeSize>::load(block,
    															   sharedRegionPtr,
    															   ghostLayerToThreadsMapping,
    															   nn_blocks,
    															   this->blockMap,
    															   stencilSupportRadius,
    															   ghostLayerSize,
    															   blockId.id);
    }

    template<unsigned int p, unsigned int ... props>
    inline __device__ void
    __loadGhost(const unsigned int blockId, const int * neighboursPos, void *sharedRegionPtr[sizeof...(props)+1])
    {
        __loadGhost<p>(blockId, neighboursPos, sharedRegionPtr);
        if (sizeof...(props) > 1)
        {
            __loadGhost<props ...>(blockId, neighboursPos, sharedRegionPtr + 1);
        }
        else if (sizeof...(props) == 1)
        {
            __loadGhost<props ...>(blockId, neighboursPos, *(sharedRegionPtr + 1));
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


#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KER_CUH
