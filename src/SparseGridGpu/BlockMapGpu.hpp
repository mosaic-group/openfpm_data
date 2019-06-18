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

template<typename AggregateBlockT, unsigned int threadBlockSize=128, typename indexT=int, template<typename> class layout_base=memory_traits_inte>
class BlockMapGpu
{
private:
    typedef BlockTypeOf<AggregateBlockT, 0> BlockT0;

protected:
    typedef typename AggregateAppend<DataBlock<unsigned char, BlockT0::size>, AggregateBlockT>::type AggregateInternalT;
    static const unsigned int pMask = AggregateInternalT::max_prop_real - 1;
    openfpm::vector_sparse_gpu<
            AggregateInternalT,
            openfpm::VECTOR_SPARSE_BLOCK,
            BlockMapGpuFunctors::BlockFunctor<threadBlockSize>
            > blockMap;

public:
    typedef AggregateBlockT AggregateType;

public:
    BlockMapGpu() = default;

//    auto get(unsigned int linId) const -> decltype(blockMap.get(0));

    template<unsigned int p>
    auto get(unsigned int linId) const -> ScalarTypeOf<AggregateBlockT, p>;

//    auto insert(unsigned int linId) -> decltype(blockMap.insert(0));

    template<unsigned int p>
    auto insert(unsigned int linId) -> ScalarTypeOf<AggregateBlockT, p> &;

    BlockMapGpu_ker<AggregateInternalT, indexT, layout_base> toKernel()
    {
        BlockMapGpu_ker<AggregateInternalT, indexT, layout_base> toKer(blockMap.toKernel());
        return toKer;
    }

    template<unsigned int ... prp>
    void deviceToHost();

    void setGPUInsertBuffer(int nBlock, int nSlot);

    template<unsigned int p>
    void initializeGPUInsertBuffer();

    template<typename ... v_reduce>
    void flush(mgpu::ofp_context_t &context, flush_type opt = FLUSH_ON_HOST, int i = 0);

    template<unsigned int p>
    void setBackground(ScalarTypeOf<AggregateBlockT, p> backgroundValue);
};

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<unsigned int p>
auto
BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::get(unsigned int linId) const -> ScalarTypeOf<AggregateBlockT, p>
{
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    unsigned int blockId = linId / BlockT::size;
    unsigned int offset = linId % BlockT::size;
    auto &block = blockMap.template get<p>(blockId);
    return block[offset];
}

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<unsigned int p>
auto
BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::insert(unsigned int linId) -> ScalarTypeOf<AggregateBlockT, p> &
{
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    unsigned int blockId = linId / BlockT::size;
    unsigned int offset = linId % BlockT::size;
    auto aggregate = blockMap.insert(blockId);
    auto &block = aggregate.template get<p>();
    auto &mask = aggregate.template get<pMask>();
    block.setElement(mask[offset]);
    return block[offset];
}

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<unsigned int ... prp>
void BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::deviceToHost()
{
    blockMap.template deviceToHost<prp..., pMask>();
}

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
void BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::setGPUInsertBuffer(int nBlock, int nSlot)
{
    // Prealloc the insert buffer on the underlying sparse vector
    blockMap.setGPUInsertBuffer(nBlock, nSlot);
}

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<unsigned int p>
void BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::initializeGPUInsertBuffer()
{
    //todo: Test if it's enough to just initialize masks to 0, without any background value
    // Initialize the blocks to background
    auto & insertBuffer = blockMap.getGPUInsertBuffer();
    std::cout << "initializeGPUInsertBuffer :: insertBuffer.size() = " << insertBuffer.size() << std::endl; //debug
    typedef BlockTypeOf<AggregateBlockT, p> BlockType; // Here assuming that all block types in the aggregate have the same size!
    constexpr unsigned int chunksPerBlock = threadBlockSize / BlockType::size; // Floor is good here...
    BlockMapGpuKernels::initializeInsertBuffer<p, pMask, chunksPerBlock> <<< insertBuffer.size()/chunksPerBlock, chunksPerBlock*BlockType::size >>>(
            insertBuffer.toKernel(),
            blockMap.template getBackground<p>()[0]
                    );
}

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<typename ... v_reduce>
void BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::flush(mgpu::ofp_context_t &context, flush_type opt, int i)
{
    blockMap.template flush<v_reduce ..., sBitwiseOr_<pMask>>(context, opt, i); // This is the one to use ideally
}

template<typename AggregateBlockT, unsigned int threadBlockSize, typename indexT, template<typename> class layout_base>
template<unsigned int p>
void BlockMapGpu<AggregateBlockT, threadBlockSize, indexT, layout_base>::setBackground(
        ScalarTypeOf<AggregateBlockT, p> backgroundValue)
{
    // NOTE: Here we assume user only passes Blocks and not scalars in the templated aggregate type
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    blockMap.template getBackground<pMask>() = 0;
    for (unsigned int i = 0; i < BlockT::size; ++i)
    {
        blockMap.template getBackground<p>()[i] = backgroundValue;
    }
}

#endif /* BLOCK_MAP_GPU_HPP_ */