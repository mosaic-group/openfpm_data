#ifndef BLOCK_MAP_GPU_HPP_
#define BLOCK_MAP_GPU_HPP_

#include "Vector/map_vector_sparse.hpp"
#include "BlockMapGpu_ker.cuh"
#include "BlockMapGpu_kernels.cuh"
#include "DataBlock.cuh"
#include <set>

template<typename BlockT, typename T>
struct AggregateAppend
{
};

template<typename BlockT, typename ... list>
struct AggregateAppend<BlockT, aggregate<list ...>>
{
    typedef aggregate<list..., BlockT> type;
};

template<typename AggregateT, unsigned int p>
using BlockTypeOf = typename std::remove_reference<typename boost::fusion::result_of::at_c<typename AggregateT::type, p>::type>::type;

template<typename AggregateT, unsigned int p>
using ScalarTypeOf = typename std::remove_reference<typename boost::fusion::result_of::at_c<typename AggregateT::type, p>::type>::type::scalarType;

template<typename AggregateBlockT, unsigned int threadBlockSize=128, typename indexT=long int, template<typename> class layout_base=memory_traits_inte>
class BlockMapGpu
{
private:
    typedef BlockTypeOf<AggregateBlockT, 0> BlockT0;

#ifdef SE_CLASS1

    //! Indicate if the setGPUInsertBuffer has been called
    bool is_setGPUInsertBuffer = false;

    //! Indicate if the initializeGPUInsertBuffer has been called
    bool is_initializeGPUInsertBuffer = false;

#endif

protected:
    const static unsigned char EXIST_BIT = 0;
    typedef typename AggregateAppend<DataBlock<unsigned char, BlockT0::size>, AggregateBlockT>::type AggregateInternalT;
    static const unsigned int pMask = AggregateInternalT::max_prop_real - 1;
    openfpm::vector_sparse_gpu_block<
            AggregateInternalT,
            BlockMapGpuFunctors::BlockFunctor<threadBlockSize>
            > blockMap;

public:
    typedef AggregateBlockT AggregateType;

    BlockMapGpu() = default;

	/*! \brief Get the background value
	 *
	 * \return background value
	 *
	 */
	auto getBackgroundValue() -> decltype(blockMap.getBackground())
	{
		return blockMap.getBackground();
	}

//    auto get(unsigned int linId) const -> decltype(blockMap.get(0));

    template<unsigned int p>
    auto get(unsigned int linId) const -> const ScalarTypeOf<AggregateBlockT, p> &
    {
        typedef BlockTypeOf<AggregateBlockT, p> BlockT;
        unsigned int blockId = linId / BlockT::size;
        unsigned int offset = linId % BlockT::size;
        auto aggregate = blockMap.get(blockId);
        auto &block = aggregate.template get<p>();
    	auto &mask = aggregate.template get<pMask>();
    	// Now check if the element actually exists
    	if (exist(mask[offset]))
    	{
    		return block[offset];
    	}
    	else
    	{
    		return blockMap.template getBackground<p>()[offset];
    	}
    }

    /*! \brief insert data, host version
     *
     * \tparam property id
     *
     * \param linId linearized id block + local linearization
     *
     * \return a reference to the data
     *
     */
    template<unsigned int p>
    auto insert(unsigned int linId) -> ScalarTypeOf<AggregateBlockT, p> &
    {
        typedef BlockTypeOf<AggregateBlockT, p> BlockT;
        unsigned int blockId = linId / BlockT::size;
        unsigned int offset = linId % BlockT::size;
        auto aggregate = blockMap.insert(blockId);
        auto &block = aggregate.template get<p>();
        auto &mask = aggregate.template get<pMask>();
        setExist(mask[offset]);
        return block[offset];
    }

    /*! \brief insert a block + flush, host version
     *
     * \tparam property id
     *
     * \param linId linearized id block
     *
     * \return a reference to the block data
     *
     */
    template<unsigned int p>
    auto insertBlockFlush(size_t blockId) -> decltype(blockMap.insertFlush(blockId).template get<p>())
    {
        typedef BlockTypeOf<AggregateBlockT, p> BlockT;

        auto aggregate = blockMap.insertFlush(blockId);
        auto &block = aggregate.template get<p>();

        return block;
    }

    /*! \brief insert a block + flush, host version
     *
     * \param linId linearized id block
     *
     * \return a reference to the block data
     *
     */
    auto insertBlockFlush(size_t blockId) -> decltype(blockMap.insertFlush(blockId))
    {
        return blockMap.insertFlush(blockId);
    }

    BlockMapGpu_ker<AggregateInternalT, indexT, layout_base> toKernel()
    {
        BlockMapGpu_ker<AggregateInternalT, indexT, layout_base> toKer(blockMap.toKernel());
        return toKer;
    }

    template<unsigned int ... prp>
    void deviceToHost()
    {
        blockMap.template deviceToHost<prp..., pMask>();
    }

    void deviceToHost();

    template<unsigned int ... prp>
    void hostToDevice();

    void hostToDevice();

    /*! \Brief Before inser any element you have to call this function to initialize the insert buffer
     *
     * \param nBlock number of blocks the insert buffer has
     * \param nSlot maximum number of insertion each thread block does
     *
     */
    void setGPUInsertBuffer(int nBlock, int nSlot)
    {
        // Prealloc the insert buffer on the underlying sparse vector
        blockMap.setGPUInsertBuffer(nBlock, nSlot);
        initializeGPUInsertBuffer();

#ifdef SE_CLASS1
        is_setGPUInsertBuffer = true;
#endif
    }

	/*! \brief In case we manually set the added index buffer and the add data buffer we have to call this
	 *         function before flush
	 *
	 *
	 */
	void preFlush()
	{
		blockMap.preFlush();
	}

    void initializeGPUInsertBuffer()
    {
        //todo: Test if it's enough to just initialize masks to 0, without any background value
        // Initialize the blocks to background
        auto & insertBuffer = blockMap.getGPUInsertBuffer();
        typedef BlockTypeOf<AggregateInternalT, pMask> BlockType; // Here assuming that all block types in the aggregate have the same size!
        constexpr unsigned int chunksPerBlock = 1; // Floor is good here...

        if (insertBuffer.size() != 0)
        {
        	CUDA_LAUNCH_DIM3((BlockMapGpuKernels::initializeInsertBuffer<pMask, chunksPerBlock>),insertBuffer.size()/chunksPerBlock, chunksPerBlock*BlockType::size,
                insertBuffer.toKernel());
        }

    #ifdef SE_CLASS1
            is_initializeGPUInsertBuffer = true;
    #endif
    }

    template<typename ... v_reduce>
    void flush(mgpu::ofp_context_t &context, flush_type opt = FLUSH_ON_HOST)
    {
#ifdef SE_CLASS1

    	if (is_setGPUInsertBuffer == false || is_initializeGPUInsertBuffer == false)
    	{std::cout << __FILE__ << ":" << __LINE__ << " error setGPUInsertBuffer you must call before doing any insertion " << std::endl;}
#endif

        blockMap.template flush<v_reduce ... >(context, opt);
    }

    /*! \brief set the background for property p
     *
     * \tparam p property p
     *
     */
    template<unsigned int p>
    void setBackgroundValue(ScalarTypeOf<AggregateBlockT, p> backgroundValue)
    {
        // NOTE: Here we assume user only passes Blocks and not scalars in the templated aggregate type
        typedef BlockTypeOf<AggregateInternalT, p> BlockT;
        typedef BlockTypeOf<AggregateInternalT, pMask> BlockM;

        BlockT bP;
        BlockM bM;

        for (unsigned int i = 0; i < BlockT::size; ++i)
        {
            bP[i] = backgroundValue;
            bM[i] = 0;
        }

        blockMap.template setBackground<p>(bP);
        blockMap.template setBackground<pMask>(bM);
    }

    template<typename BitMaskT>
	inline static bool getBit(const BitMaskT &bitMask, unsigned char pos)
	{
		return (bitMask>>pos)&1U;
	}

	template<typename BitMaskT>
	inline static bool setBit(BitMaskT &bitMask, unsigned char pos)
	{
		return bitMask |= 1U<<pos;
	}

	template<typename BitMaskT>
	inline static bool unsetBit(BitMaskT &bitMask, unsigned char pos)
	{
		return bitMask &= !(1U<<pos);
	}

    template<typename BitMaskT>
    inline static bool exist(BitMaskT &bitMask)
    {
        return getBit(bitMask, EXIST_BIT);
    }

    template<typename BitMaskT>
    inline static void setExist(BitMaskT &bitMask)
    {
        setBit(bitMask, EXIST_BIT);
    }

    template<typename BitMaskT>
    inline static void unsetExist(BitMaskT &bitMask)
    {
        unsetBit(bitMask, EXIST_BIT);
    }


    /*! \brief Return internal structure block map
     *
     * \return the blockMap
     *
     */
    decltype(blockMap) & private_get_blockMap_non_const()
	{
    	return blockMap;
	}

    /*! \brief Return internal structure block map
     *
     * \return the blockMap
     *
     */
    decltype(blockMap) & private_get_blockMap()
	{
    	return blockMap;
	}

    /*! \brief Return internal structure block map
     *
     * \return the blockMap
     *
     */
    const decltype(blockMap) & private_get_blockMap() const
	{
    	return blockMap;
	}
};

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
void BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::deviceToHost()
{
    blockMap.template deviceToHost<pMask>();
}

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<unsigned int ... prp>
void BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::hostToDevice()
{
    blockMap.template hostToDevice<prp..., pMask>();
}

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
void BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::hostToDevice()
{
    blockMap.template hostToDevice<pMask>();
}

//template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
//template<unsigned int p>
//void BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::setBackgroundValue(
//        ScalarTypeOf<AggregateBlockT, p> backgroundValue)


#endif /* BLOCK_MAP_GPU_HPP_ */
