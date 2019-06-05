#ifndef SPARSE_GRID_GPU_KER_CUH_
#define SPARSE_GRID_GPU_KER_CUH_

#include "util/cuda_util.hpp"
#include <cstdlib>
#include "Vector/map_vector_sparse.hpp"
#include "DataBlock.cuh"

template<typename AggregateT, unsigned int p>
using BlockTypeOf = typename std::remove_reference<typename boost::fusion::result_of::at_c<typename AggregateT::type, p>::type>::type;

template<typename AggregateT, unsigned int p>
using ScalarTypeOf = typename std::remove_reference<typename boost::fusion::result_of::at_c<typename AggregateT::type, p>::type>::type::scalarType;

template <typename AggregateT>
struct LastPOf
{
    static const unsigned int value = AggregateT::max_prop_real - 1;
};

template <typename AggregateT, unsigned int pMask>
struct InsertBlockWrapper
{
    AggregateT aggregate;

    InsertBlockWrapper() = default;

    InsertBlockWrapper(AggregateT aggregate) : aggregate(aggregate) {}

    InsertBlockWrapper(const InsertBlockWrapper<AggregateT, pMask> &other)
    {
        aggregate = other.aggregate;
    }

    InsertBlockWrapper<AggregateT, pMask> &operator=(const InsertBlockWrapper<AggregateT, pMask> &other)
    {
        aggregate = other.aggregate;
        return *this;
    }

    template <unsigned int p>
    inline auto get() -> decltype(aggregate.template get<p>())
    {
        return aggregate.template get<p>();
    }

    inline auto getMask() -> decltype(aggregate.template get<pMask>())
    {
        return aggregate.template get<pMask>();
    }
};

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
class SparseGridGpu_ker
{
private:
    openfpm::vector_sparse_gpu_ker<AggregateBlockT, indexT, layout_base> blockMap;

public:
    static const unsigned int pMask = AggregateBlockT::max_prop_real - 1;
    typedef AggregateBlockT AggregateType;
    typedef InsertBlockWrapper<AggregateBlockT, pMask> InsertBlockWrapperType;

public:
    SparseGridGpu_ker(openfpm::vector_sparse_gpu_ker<AggregateBlockT, indexT, layout_base> blockMap)
            : blockMap(blockMap) {};

    template<unsigned int p>
    inline __device__ auto get(unsigned int linId) const -> ScalarTypeOf<AggregateBlockT, p>;

    template<unsigned int p>
    inline __device__ auto get(unsigned int blockId, unsigned int offset) const -> ScalarTypeOf<AggregateBlockT, p>;

    template<unsigned int p>
    inline __device__ auto insert(unsigned int linId) -> ScalarTypeOf<AggregateBlockT, p>&;

    template<unsigned int p>
    inline __device__ auto insert(unsigned int blockId, unsigned int offset) -> ScalarTypeOf<AggregateBlockT, p>&;

//    inline __device__ auto insertBlock(unsigned int blockId) -> BlockTypeOf<AggregateBlockT, p>*;
    inline __device__ auto insertBlock(unsigned int blockId) -> decltype(blockMap.insert(0));

    inline __device__ void init()
    {
        blockMap.init();
    }

    inline __device__ void flush_block_insert()
    {
        blockMap.flush_block_insert();
    }
};

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ auto SparseGridGpu_ker<AggregateBlockT, indexT, layout_base>
::get(unsigned int linId) const -> ScalarTypeOf<AggregateBlockT, p>
{
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    unsigned int blockId = linId / BlockT::size;
    unsigned int offset = linId % BlockT::size;
    return get<p>(blockId, offset);
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ auto SparseGridGpu_ker<AggregateBlockT, indexT, layout_base>
::get(unsigned int blockId, unsigned int offset) const -> ScalarTypeOf<AggregateBlockT, p>
{
    auto &block = blockMap.template get<p>(blockId);
    return block[offset];
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ auto SparseGridGpu_ker<AggregateBlockT, indexT, layout_base>
::insert(unsigned int linId) -> ScalarTypeOf<AggregateBlockT, p>&
{
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    unsigned int blockId = linId / BlockT::size;
    unsigned int offset = linId % BlockT::size;
    return insert<p>(blockId, offset);
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ auto SparseGridGpu_ker<AggregateBlockT, indexT, layout_base>
::insert(unsigned int blockId, unsigned int offset) -> ScalarTypeOf<AggregateBlockT, p>&
{
    auto &block = blockMap.template insert<p>(blockId);
    auto &mask = blockMap.template insert<pMask>(blockId);
    block.setElement(mask[offset]);
    return block[offset];
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
inline __device__ auto SparseGridGpu_ker<AggregateBlockT, indexT, layout_base>
        ::insertBlock(unsigned int blockId) -> decltype(blockMap.insert(0))
{
    return blockMap.insert(blockId);
}

#endif /* SPARSE_GRID_GPU_KER_CUH_ */