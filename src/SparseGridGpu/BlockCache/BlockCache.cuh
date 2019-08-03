//
// Created by tommaso on 27/06/19.
//

#ifndef OPENFPM_PDATA_BLOCKCACHE_CUH
#define OPENFPM_PDATA_BLOCKCACHE_CUH

#include "util/cuda_util.hpp"

namespace BlockCacheUtils
{
    template <typename T>
    inline __device__ __host__ unsigned int SetToZeroIfFalse<true, T>(T value)
    {
        return value;
    };

    template <typename T>
    inline __device__ __host__ unsigned int SetToZeroIfFalse<false, T>(T value)
    {
        return 0;
    };
}

/**
 * BlockCache is an abstraction built on the concept of loading a block into shared
 * memory before using it in a stencil operation.
 * The idea is to provide a way to transparently address shared and global memory via coordinates,
 * caching the block data into shared memory but also allowing addressing non-cached data directly in global.
 */
template <typename SparseGridT, unsigned int chunksPerBlock, bool loadGhostInSharedMemory, unsigned int ... props>
struct BlockCache
{
    static void
};

#endif //OPENFPM_PDATA_BLOCKCACHE_CUH
