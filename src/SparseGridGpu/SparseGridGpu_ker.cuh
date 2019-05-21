#ifndef SPARSE_GRID_GPU_KER_CUH_
#define SPARSE_GRID_GPU_KER_CUH_

#include <host_defines.h>
#include <cstdlib>
#include "Vector/map_vector_sparse.hpp"
#include "DataBlock.cuh"

template<typename AggregateT, unsigned int p>
using BlockTypeOf = typename std::remove_reference<typename boost::fusion::result_of::at_c<typename AggregateT::type, p>::type>::type;

template<typename AggregateT, unsigned int p>
using ScalarTypeOf = typename std::remove_reference<typename boost::fusion::result_of::at_c<typename AggregateT::type, p>::type>::type::scalarType;

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
class SparseGridGpu_ker
{
private:
    openfpm::vector_sparse_gpu_ker<AggregateBlockT, indexT, layout_base> blockMap;

public:
    typedef AggregateBlockT AggregateType;

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

    template<unsigned int p>
    inline __device__ auto insertBlock(unsigned int blockId) -> BlockTypeOf<AggregateBlockT, p>*;

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
    block.setElementDevice(offset);
    return block[offset];
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ auto SparseGridGpu_ker<AggregateBlockT, indexT, layout_base>
::insertBlock(unsigned int blockId) -> BlockTypeOf<AggregateBlockT, p>*
{
    auto &block = blockMap.template insert<p>(blockId);
    block.existBitMask = 0;
    return &block;
}

//template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
//inline __device__ void init()
//{
//    blockMap.init();
//}

//template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
//inline __device__ void flush_block_insert()
//{
//    blockMap.flush_block_insert();
//}

#endif /* SPARSE_GRID_GPU_KER_CUH_ */