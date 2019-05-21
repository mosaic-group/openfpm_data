//
// Created by tommaso on 14/05/19.
//

#ifndef OPENFPM_PDATA_DATABLOCK_CUH
#define OPENFPM_PDATA_DATABLOCK_CUH

#include <cstdlib>
#include <host_defines.h>
#include <cstring>

//todo: Copy the new DataBlock definition here (the one where bitmask is external)
//todo: Rename static exist/setBit methods to getBit and setBit which can work with generic bitmasks

template<typename ScalarT, unsigned int DataBlockSize=64, typename BitMaskT=unsigned char>
struct DataBlock
{
    typedef ScalarT scalarType;

    static const unsigned int EXISTBIT = 0;

    static const unsigned int size = DataBlockSize;
    ScalarT block[size];

    __device__ __host__ DataBlock() = default;

    __device__ __host__ DataBlock(const DataBlock &other)
    {
#ifdef  __CUDA_ARCH__
        block[threadIdx.x % size] = other.block[threadIdx.x % size];
#else
        memcpy(block, other.block, size * sizeof(ScalarT));
#endif
    }

    __device__ __host__ DataBlock operator=(const DataBlock &other)
    {
#ifdef  __CUDA_ARCH__
        block[threadIdx.x % size] = other.block[threadIdx.x % size];
#else
        memcpy(block, other.block, size * sizeof(ScalarT));
#endif
        return *this;
    }

    __device__ __host__ DataBlock operator=(ScalarT v)
    {
#ifdef  __CUDA_ARCH__
        block[threadIdx.x % size] = v;
#else
        for (unsigned int i = 0; i < size; ++i)
        {
            block[i] = v;
        }
#endif
        return *this;
    }

    __device__ __host__ inline ScalarT &operator[](unsigned int i)
    {
        return block[i];
    }

    __device__ __host__ inline const ScalarT &operator[](unsigned int i) const
    {
        return block[i];
    }

    __device__ __host__ inline static bool exist(BitMaskT &bitMask)
    {
        return (bitMask >> EXISTBIT) & 1U;
    }

    __device__ __host__ inline static void setElement(BitMaskT &bitMask)
    {
        bitMask = bitMask | ( 1U << EXISTBIT );
    }
};

#endif //OPENFPM_PDATA_DATABLOCK_CUH
