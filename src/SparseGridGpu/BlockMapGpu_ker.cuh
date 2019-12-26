#ifndef BLOCK_MAP_GPU_KER_CUH_
#define BLOCK_MAP_GPU_KER_CUH_

#include "util/cuda_util.hpp"
#include <cstdlib>
#include "Vector/map_vector_sparse.hpp"
#include "DataBlock.cuh"
#include "TemplateUtils/encap_shmem.hpp"

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
#if defined(__NVCC__) || defined(__HIPCC__)
        aggregate = other.aggregate;
#else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
    }

    InsertBlockWrapper<AggregateT, pMask> &operator=(const InsertBlockWrapper<AggregateT, pMask> &other)
    {
#if defined(__NVCC__) || defined(__HIPCC__)
        aggregate = other.aggregate;
        return *this;
#else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
    }

    template <unsigned int p>
    inline auto get() -> decltype(aggregate.template get<p>())
    {
#if defined(__NVCC__) || defined(__HIPCC__)
        return aggregate.template get<p>();
#else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
    }

    inline auto getMask() -> decltype(aggregate.template get<pMask>())
    {
#if defined(__NVCC__) || defined(__HIPCC__)
        return aggregate.template get<pMask>();
#else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
    }
};

template<typename AggregateBlockT=aggregate<DataBlock<float, 64>>, typename indexT=int, template<typename> class layout_base=memory_traits_inte>
class BlockMapGpu_ker
{
protected:
    openfpm::vector_sparse_gpu_ker<AggregateBlockT, indexT, layout_base> blockMap;
    const static unsigned char EXIST_BIT = 0;

public:
    static const unsigned int pMask = AggregateBlockT::max_prop_real - 1;
    typedef AggregateBlockT AggregateType;
    typedef InsertBlockWrapper<AggregateBlockT, pMask> InsertBlockWrapperType;

public:
    template<typename BitMaskT>
    inline static __device__ __host__ bool getBit(const BitMaskT &bitMask, unsigned char pos)
    {
        return (bitMask>>pos)&1U;
    }

    template<typename BitMaskT>
    inline static __device__ __host__ void setBit(BitMaskT &bitMask, unsigned char pos)
    {
        bitMask = bitMask | (1U<<pos);
    }

    template<typename BitMaskT>
    inline static __device__ __host__ void unsetBit(BitMaskT &bitMask, unsigned char pos)
    {
        bitMask = bitMask & ~(1U<<pos);
    }

public:
    BlockMapGpu_ker(openfpm::vector_sparse_gpu_ker<AggregateBlockT, indexT, layout_base> blockMap)
            : blockMap(blockMap) {};

    template<unsigned int p>
    inline __device__ auto get(unsigned int linId) const -> ScalarTypeOf<AggregateBlockT, p>;

    template<unsigned int p>
    inline __device__ auto get(unsigned int blockId, unsigned int offset) const -> ScalarTypeOf<AggregateBlockT, p>;

    inline __device__ auto getBlock(unsigned int blockId) -> decltype(blockMap.get(0));

    template<unsigned int p>
    inline __device__ ScalarTypeOf<AggregateBlockT, p> & getReference(unsigned int linId);

    template<unsigned int p>
    inline __device__ ScalarTypeOf<AggregateBlockT, p> & getReference(unsigned int blockId, unsigned int offset);

    template<unsigned int p>
    inline __device__ auto insert(unsigned int linId) -> ScalarTypeOf<AggregateBlockT, p>&
    {
    #if defined(__NVCC__) || defined(__HIPCC__)
        typedef BlockTypeOf<AggregateBlockT, p> BlockT;
        unsigned int blockId = linId / BlockT::size;
        unsigned int offset = linId % BlockT::size;
        return insert<p>(blockId, offset);
    #else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
    #endif //
    }

    template<unsigned int p>
    inline __device__ auto insert(unsigned int blockId, unsigned int offset) -> ScalarTypeOf<AggregateBlockT, p>&;

    template<unsigned int nChunksPerBlocks = 1>
    inline __device__ auto insertBlock(indexT blockId, unsigned int stride = 8192) -> decltype(blockMap.insert(0))
	{
    	int offset = threadIdx.x / stride;
//    	__shared__ int mem[nChunksPerBlocks][encap_shmem<sizeof(blockMap.insert(0))>::nthr];
    	__shared__ int mem_[nChunksPerBlocks];

    	decltype(blockMap.insert(0)) ec_(blockMap.private_get_data(),0);

	#if defined(__NVCC__) || defined(__HIPCC__)
    	if (threadIdx.x % stride == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    	{
    		auto ec = blockMap.insert(blockId);

    		mem_[offset] = ec.private_get_k();

    		// copy to shared to broadcast on all thread
    		//new (mem[offset]) decltype(ec)(ec.private_get_data(),ec.private_get_k());
    	}

    	__syncthreads();;

    	ec_.private_set_k(mem_[offset]);

		return ec_/* *(decltype(blockMap.insert(0)) *)mem[offset]*/;
		#else //
		    std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
		#endif //
	}

    inline __device__ openfpm::vector_sparse_gpu_ker<AggregateBlockT, indexT, layout_base> & getblockMap()
    {
    	return blockMap;
    }

    inline static __device__ unsigned int getBlockId(unsigned int linId)
    {
#if defined(__NVCC__) || defined(__HIPCC__)
        return linId / BlockTypeOf<AggregateBlockT, 0>::size;
#else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
    }

    inline static __device__ unsigned int getOffset(unsigned int linId)
    {
#if defined(__NVCC__) || defined(__HIPCC__)
        return linId % BlockTypeOf<AggregateBlockT, 0>::size;
#else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
    }

    inline __device__ void init()
    {
#if defined(__NVCC__) || defined(__HIPCC__)
        blockMap.init();
#else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
    }

    inline __device__ void flush_block_insert()
    {
#if defined(__NVCC__) || defined(__HIPCC__)
        blockMap.flush_block_insert();
#else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
    }

    template<typename BitMaskT>
    inline static __device__ bool exist(const BitMaskT &bitMask)
    {
        return getBit(bitMask, EXIST_BIT);
    }

    template<typename BitMaskT>
    inline static __device__ void setExist(BitMaskT &bitMask)
    {
        setBit(bitMask, EXIST_BIT);
    }

    template<typename BitMaskT>
    inline static __device__ void unsetExist(BitMaskT &bitMask)
    {
        unsetBit(bitMask, EXIST_BIT);
    }

    inline __device__ ScalarTypeOf<AggregateBlockT, pMask> getMask(unsigned int linId) const
    {
        return get<pMask>(linId);
    }

    inline __device__ void remove(unsigned int linId)
    {
    #if defined(__NVCC__) || defined(__HIPCC__)
        typedef BlockTypeOf<AggregateBlockT, pMask> BlockT;
        unsigned int blockId = linId / BlockT::size;
        unsigned int offset = linId % BlockT::size;
        remove(blockId, offset);
    #else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
    #endif //
    }

    inline __device__ void remove(unsigned int blockId, unsigned int offset)
    {
    #if defined(__NVCC__) || defined(__HIPCC__)
    //    const auto & aggregate = blockMap.get(blockId);
    //    const auto & block = aggregate.template get<p>();
    //    const auto & mask = aggregate.template get<pMask>();
    //    // Now check if the element actually exists
    //    return exist(mask[offset])
    //                ? block[offset]
    //                : blockMap.template getBackground<p>()[offset];
    ////    return blockMap.template get<p>(blockId)[offset];

        const auto sid = blockMap.get_sparse(blockId);
        blockMap.template get<pMask>(sid)[offset] = 0;

    #else //
        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
    #endif //
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

			pc = blockMap.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("blockMap overflow : ") + "\n" + pc.match_str;
				return pc;
			}

			return pc;
		}

#endif

};

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ auto BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>
::get(unsigned int linId) const -> ScalarTypeOf<AggregateBlockT, p>
{
#if defined(__NVCC__) || defined(__HIPCC__)
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    unsigned int blockId = linId / BlockT::size;
    unsigned int offset = linId % BlockT::size;
    return get<p>(blockId, offset);
#else //
    std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ auto BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>
::get(unsigned int blockId, unsigned int offset) const -> ScalarTypeOf<AggregateBlockT, p>
{
#if defined(__NVCC__) || defined(__HIPCC__)

    const auto sid = blockMap.get_sparse(blockId);
    const auto & block = blockMap.template get_ele<p>(sid.id)[offset];
    const auto mask = blockMap.template get_ele<pMask>(sid.id)[offset];
    // Now check if the element actually exists
    return exist(mask)
                ? block
                : blockMap.template getBackground<p>()[offset];

#else //
    std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
}


template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
inline __device__ auto BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>
::getBlock(unsigned int blockId) -> decltype(blockMap.get(0))
{
#if defined(__NVCC__) || defined(__HIPCC__)
    return blockMap.get(blockId);
#else //
    std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ ScalarTypeOf<AggregateBlockT, p> & BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>
::getReference(const unsigned int linId)
{
    // Only call this if you are TOTALLY SURE the element exists! Otherwise KABOOOOOM! :D
#if defined(__NVCC__) || defined(__HIPCC__)
    typedef BlockTypeOf<AggregateBlockT, p> BlockT;
    unsigned int blockId = linId / BlockT::size;
    unsigned int offset = linId % BlockT::size;
    return getReference<p>(blockId, offset);
#else //
    std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
}

template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ ScalarTypeOf<AggregateBlockT, p> & BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>
::getReference(const unsigned int blockId, const unsigned int offset)
{
    // Only call this if you are TOTALLY SURE the element exists! Otherwise KABOOOOOM! :D
#if defined(__NVCC__) || defined(__HIPCC__)
    return blockMap.template get<p>(blockId)[offset];
#else //
    std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
}


template<typename AggregateBlockT, typename indexT, template<typename> class layout_base>
template<unsigned int p>
inline __device__ auto BlockMapGpu_ker<AggregateBlockT, indexT, layout_base>
::insert(unsigned int blockId, unsigned int offset) -> ScalarTypeOf<AggregateBlockT, p>&
{
#if defined(__NVCC__) || defined(__HIPCC__)
    auto aggregate = blockMap.insert(blockId);
    auto &block = aggregate.template get<p>();
    auto &mask = aggregate.template get<pMask>();
    setExist(mask[offset]);
    return block[offset];
#else //
    std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
#endif //
}


#endif /* BLOCK_MAP_GPU_KER_CUH_ */
