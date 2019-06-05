//
// Created by tommaso on 14/05/19.
//

#ifndef OPENFPM_PDATA_DATABLOCK_CUH
#define OPENFPM_PDATA_DATABLOCK_CUH

#include <cstdlib>
#include "util/cuda_util.hpp"
#include <cstring>

//todo: Copy the new DataBlock definition here (the one where bitmask is external)
//todo: Rename static exist/setBit methods to getBit and setBit which can work with generic bitmasks

template<typename ScalarT, unsigned int DataBlockSize=64>
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
#ifdef __NVCC__
        block[threadIdx.x % size] = other.block[threadIdx.x % size];
#endif // __NVCC__
#else // __CUDA_ARCH__
        memcpy(block, other.block, size * sizeof(ScalarT));
#endif // __CUDA_ARCH__
    }

    __device__ __host__ DataBlock operator=(const DataBlock &other)
    {
#ifdef  __CUDA_ARCH__
#ifdef __NVCC__
        block[threadIdx.x % size] = other.block[threadIdx.x % size];
#endif // __NVCC__
#else
        memcpy(block, other.block, size * sizeof(ScalarT));
#endif
        return *this;
    }

    __device__ __host__ DataBlock operator=(ScalarT v)
    {
#ifdef  __CUDA_ARCH__
#ifdef __NVCC__
        block[threadIdx.x % size] = v;
#endif // __NVCC__
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

    template<typename BitMaskT>
    __device__ __host__ inline static bool exist(BitMaskT &bitMask)
    {
        return bitMask==1;
    }

    template<typename BitMaskT>
    __device__ __host__ inline static void setElement(BitMaskT &bitMask)
    {
        bitMask = 1;
    }
};

#endif //OPENFPM_PDATA_DATABLOCK_CUH
