//
// Created by tommaso on 17/06/19.
//

#ifndef OPENFPM_PDATA_BLOCKGEOMETRY_HPP
#define OPENFPM_PDATA_BLOCKGEOMETRY_HPP

#include <boost/mpl/size_t.hpp>
#include <cstring>
#include <Grid/grid_sm.hpp>

template<unsigned int base, unsigned int exponent>
struct IntPow
{
    constexpr static size_t value = base * IntPow<base, exponent - 1>::value;
};

template<unsigned int base>
struct IntPow<base, 0>
{
    constexpr static size_t value = 1;
};

/**
 * This class provides an interface to linearization of coordinates and viceversa when blocks are involved.
 * This can be seen as a lightweight version of grid_sm, with just LinId and InvLinId methods, but
 * tuned for blocked data.
 */
template<unsigned int dim, unsigned int blockEdgeSize>
class BlockGeometry
{
private:
    constexpr static size_t blockSize = IntPow<blockEdgeSize, dim>::value;
    size_t blockSz[dim];
//    size_t sz[dim];

public:
    //todo: Add a constructor from dim3 and uint3 objects

    __host__ __device__ BlockGeometry(const size_t blockDimensions[dim])
    {
        memcpy(blockSz, blockDimensions, dim * sizeof(size_t));
//        for (int d=0; d<dim; ++d)
//        {
//            sz[d] = blockDimensions[d] * blockEdgeSize;
//        }
    }

    __host__ __device__ BlockGeometry(const size_t domainBlockEdgeSize)
    {
        for (int i = 0; i < dim; ++i)
        {
            blockSz[i] = domainBlockEdgeSize;
//            sz[i] = domainBlockEdgeSize * blockEdgeSize;
        }
    }

    template<typename T>
    __host__ __device__ BlockGeometry(const grid_sm<dim, T> blockGrid)
    {
        for (int i = 0; i < dim; ++i)
        {
            blockSz[i] = blockGrid.size(i);
//            sz[i] = blockGrid.size(i) * blockEdgeSize;
        }
    }

    __host__ __device__ BlockGeometry(const BlockGeometry<dim, blockEdgeSize> &other)
    {
        memcpy(blockSz, other.blockSz, dim * sizeof(size_t));
//        memcpy(sz, other.sz, dim * sizeof(size_t));
    }

    __host__ __device__ BlockGeometry &operator=(const BlockGeometry<dim, blockEdgeSize> &other)
    {
        if (&other != this)
        {
            memcpy(blockSz, other.blockSz, dim * sizeof(size_t));
//            memcpy(sz, other.sz, dim * sizeof(size_t));
        }
        return *this;
    }

    template<typename indexT>
    inline __host__ __device__ mem_id LinId(const grid_key_dx<dim, indexT> coord) const
    {
        mem_id blockLinId = coord.get(dim - 1) / blockEdgeSize;
        mem_id localLinId = coord.get(dim - 1) % blockEdgeSize;
        for (int d = dim - 2; d >= 0; --d)
        {
            blockLinId *= blockSz[d];
            localLinId *= blockEdgeSize;
            blockLinId += coord.get(d) / blockEdgeSize;
            localLinId += coord.get(d) % blockEdgeSize;
        }
        return blockLinId * blockSize + localLinId;
    }

    inline __host__ __device__ grid_key_dx<dim> InvLinId(const mem_id linId) const
    {
        mem_id blockLinId = linId / blockSize;
        mem_id localLinId = linId % blockSize;
        grid_key_dx<dim> coord;
        for (int d = 0; d < dim; ++d)
        {
            auto c = blockLinId % blockSz[d];
            c *= blockEdgeSize;
            c += localLinId % blockEdgeSize;
            coord.set_d(d, c);
            blockLinId /= blockSz[d];
            localLinId /= blockEdgeSize;
        }
        return coord;
    }
};

#endif //OPENFPM_PDATA_BLOCKGEOMETRY_HPP
