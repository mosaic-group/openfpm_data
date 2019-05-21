//
// Created by tommaso on 14/05/19.
//

#ifndef OPENFPM_PDATA_DATABLOCK_CUH
#define OPENFPM_PDATA_DATABLOCK_CUH

#include <cstdlib>
#include <host_defines.h>
#include <cstring>

static const unsigned int DATA_BLOCK_SIZE = 64;

template<typename ScalarT>
struct DataBlock
{
    typedef ScalarT scalarType;

    static const unsigned int size = DATA_BLOCK_SIZE;
    ScalarT block[size];
    size_t existBitMask; // note: size_t is usually 64 bits on 64bit archs
    //todo Should we also add a paddingBitMask here?
//    size_t paddingBitMask;
//    size_t propagateBitMask;

    __device__ __host__ DataBlock() = default;
//    __device__ __host__ DataBlock(){};

    __device__ __host__ DataBlock(const DataBlock &other)
            : existBitMask(other.existBitMask)
    {
        memcpy(block, other.block, size * sizeof(ScalarT));
//        memcpy(&existBitMask, &(other.existBitMask), sizeof(size_t));
        existBitMask = other.existBitMask;
    }

    __device__ __host__ DataBlock operator=(const DataBlock &other)
    {
        memcpy(block, other.block, size * sizeof(ScalarT));
//        memcpy(&existBitMask, &(other.existBitMask), sizeof(size_t));
        existBitMask = other.existBitMask;
        return *this;
    }

    __device__ __host__ DataBlock operator=(float v) // Hack to make things compile
    {
        existBitMask = 0;
        return *this;
    }

    __device__ __host__ ScalarT &operator[](unsigned int i)
    {
        return block[i];
    }

    __device__ __host__ const ScalarT &operator[](unsigned int i) const
    {
        return block[i];
    }

    /**
     * Query the bitmask if the i-th element exists.
     * @param i
     * @return True if the i-th element exists.
     */
    __device__ __host__ bool exist(unsigned int i)
    {
        return (existBitMask >> i) & ((size_t) 1);
    }

    /**
     * Set the i-th element as existing in the bitmask.
     * @param i
     */
    __host__ void setElement(unsigned int i)
    {
        existBitMask = existBitMask | ((size_t) 1) << i;
    }

    /**
     * Set the i-th element as existing in the bitmask. Gpu version
     * @param i
     */
    __device__ void setElementDevice(unsigned int i)
    {
        existBitMask = existBitMask | ((size_t) 1) << i;
//        atomicOr((unsigned long long int *)&existBitMask, (unsigned long long int)((size_t) 1) << i);
    }
};

#endif //OPENFPM_PDATA_DATABLOCK_CUH
