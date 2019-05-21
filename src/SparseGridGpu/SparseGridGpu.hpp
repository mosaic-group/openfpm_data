#ifndef SPARSE_GRID_GPU_HPP_
#define SPARSE_GRID_GPU_HPP_

#include "Vector/map_vector_sparse.hpp"
#include "SparseGridGpu_ker.cuh"
#include "SparseGridGpu_kernels.cuh"
#include "DataBlock.cuh"
#include <set>

template<typename AggregateT, unsigned int p>
using BlockTypeOf = typename std::remove_reference<typename boost::fusion::result_of::at_c<typename AggregateT::type, p>::type>::type;

template<typename AggregateT, unsigned int p>
using ScalarTypeOf = typename std::remove_reference<typename boost::fusion::result_of::at_c<typename AggregateT::type, p>::type>::type::scalarType;

template<typename AggregateBlockT, typename indexT=int, template<typename> class layout_base=memory_traits_inte>
class SparseGridGpu
{
private:
    openfpm::vector_sparse_gpu<AggregateBlockT> blockMap;

public:
    typedef AggregateBlockT AggregateType;

public:
    SparseGridGpu() = default;

    template<unsigned int p>
    auto get(unsigned int linId) const -> ScalarTypeOf<AggregateBlockT, p>;

    template<unsigned int p>
    auto insert(unsigned int linId) -> ScalarTypeOf<AggregateBlockT, p> &;

    SparseGridGpu_ker<AggregateBlockT, indexT, layout_base> toKernel()
    {
        SparseGridGpu_ker<AggregateBlockT, indexT, layout_base> toKer(blockMap.toKernel());
        return toKer;
    }

    template<unsigned int ... prp>
    void deviceToHost();

    void setGPUInsertBuffer(int nBlock, int nSlot);

    template<unsigned int p>
//    void initializeGPUInsertBuffer(unsigned int gridSize, unsigned int blockSize);
    void initializeGPUInsertBuffer();

    template<typename ... v_reduce>
    void flush(mgpu::ofp_context_t &context, flush_type opt = FLUSH_ON_HOST, int i = 0);

    template<unsigned int p>
    void setBackground(ScalarTypeOf<AggregateBlockT, p> backgroundValue);
};

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
auto
SparseGridGpu<AggregateBlockT, indexT, layout_base>::get(unsigned int linId) const -> ScalarTypeOf<AggregateBlockT, p>
{
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    unsigned int blockId = linId / BlockT::size;
    unsigned int offset = linId % BlockT::size;
    auto &block = blockMap.template get<p>(blockId);
    return block[offset];
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
auto
SparseGridGpu<AggregateBlockT, indexT, layout_base>::insert(unsigned int linId) -> ScalarTypeOf<AggregateBlockT, p> &
{
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    unsigned int blockId = linId / BlockT::size;
    unsigned int offset = linId % BlockT::size;
    auto &block = blockMap.template insert<p>(blockId);
    block.setElement(offset);
    return block[offset];
}

//template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
//SparseGridGpu_ker<AggregateBlockT, indexT, layout_base> SparseGridGpu<AggregateBlockT, indexT, layout_base>::toKernel()
//{
//    SparseGridGpu_ker<AggregateBlockT, indexT, layout_base> toKer(blockMap.toKernel());
//    return toKer;
//}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int ... prp>
void SparseGridGpu<AggregateBlockT, indexT, layout_base>::deviceToHost()
{
    blockMap.template deviceToHost<prp...>();
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
void SparseGridGpu<AggregateBlockT, indexT, layout_base>::setGPUInsertBuffer(int nBlock, int nSlot)
{
    // Prealloc the insert buffer on the underlying sparse vector
    blockMap.setGPUInsertBuffer(nBlock, nSlot);
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
//void SparseGridGpu<AggregateBlockT, indexT, layout_base>::initializeGPUInsertBuffer(unsigned int gridSize, unsigned int blockSize)
void SparseGridGpu<AggregateBlockT, indexT, layout_base>::initializeGPUInsertBuffer()
{
    // Initialize the blocks to background
    auto & insertBuffer = blockMap.getGPUInsertBuffer();
    typedef BlockTypeOf<AggregateBlockT, 0> BlockType; // Here assuming that all block types in the aggregate have the same size!
//    SparseGridGpuKernels::initializeInsertBuffer<p> <<< gridSize * BlockType::size, blockSize >>>(
    SparseGridGpuKernels::initializeInsertBuffer<p> <<< insertBuffer.size()/2, 2*BlockType::size >>>(
            insertBuffer.toKernel(),
            blockMap.template getBackground<p>()[0]
                    );
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<typename ... v_reduce>
void SparseGridGpu<AggregateBlockT, indexT, layout_base>::flush(mgpu::ofp_context_t &context, flush_type opt, int i)
{
    blockMap.template flush<v_reduce ...>(context, opt, i);
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
void SparseGridGpu<AggregateBlockT, indexT, layout_base>::setBackground(
        ScalarTypeOf<AggregateBlockT, p> backgroundValue)
{
    // NOTE: Here we assume user only passes Blocks and not scalars in the templated aggregate type
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    blockMap.template getBackground<p>().existBitMask = 0;
    for (unsigned int i = 0; i < BlockT::size; ++i)
    {
        blockMap.template getBackground<p>()[i] = backgroundValue;
    }
}

#endif /* SPARSE_GRID_GPU_HPP_ */