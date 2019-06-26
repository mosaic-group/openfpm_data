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
    size_t sz[dim];

public:

    BlockGeometry()	{}

    __host__ __device__ BlockGeometry(const size_t blockDimensions[dim])
    {
        memcpy(blockSz, blockDimensions, dim * sizeof(size_t));
        for (int d=0; d<dim; ++d)
        {
            sz[d] = blockDimensions[d] * blockEdgeSize;
        }
    }

    __host__ __device__ BlockGeometry(const size_t domainBlockEdgeSize)
    {
        for (int i = 0; i < dim; ++i)
        {
            blockSz[i] = domainBlockEdgeSize;
            sz[i] = domainBlockEdgeSize * blockEdgeSize;
        }
    }

    template<typename T>
    __host__ __device__ BlockGeometry(const grid_sm<dim, T> blockGrid)
    {
        for (int i = 0; i < dim; ++i)
        {
            blockSz[i] = blockGrid.size(i);
            sz[i] = blockGrid.size(i) * blockEdgeSize;
        }
    }

#ifdef __NVCC__
    //Constructors from dim3 and uint3 objects
    __host__ __device__ BlockGeometry(const dim3 blockDimensions)
    {
    	unsigned int i = 0;
        assert(dim <= 3);
        blockSz[i] = blockDimensions.x;
        sz[i] = blockSz[i] * blockEdgeSize;
        if (dim > 1)
        {
        	++i;
            blockSz[i] = blockDimensions.y;
            sz[i] = blockSz[i] * blockEdgeSize;
            if (dim > 2)
            {
            	++i;
                blockSz[i] = blockDimensions.z;
                sz[i] = blockSz[i] * blockEdgeSize;
            }
        }
    }

/*    __host__ __device__ BlockGeometry(const uint3 blockDimensions)
    {
        assert(dim <= 3);
        blockSz[0] = blockDimensions.x;
        sz[i] = blockSz[0] * blockEdgeSize;
        if (dim > 1)
        {
            blockSz[1] = blockDimensions.y;
            if (dim > 2)
            {
                blockSz[2] = blockDimensions.z;
            }
        }
    }*/

#endif // __NVCC__

    __host__ __device__ BlockGeometry(const BlockGeometry<dim, blockEdgeSize> &other)
    {
        memcpy(blockSz, other.blockSz, dim * sizeof(size_t));
        memcpy(sz, other.sz, dim * sizeof(size_t));
    }

    __host__ __device__ BlockGeometry &operator=(const BlockGeometry<dim, blockEdgeSize> &other)
    {
        if (&other != this)
        {
            memcpy(blockSz, other.blockSz, dim * sizeof(size_t));
            memcpy(sz, other.sz, dim * sizeof(size_t));
        }
        return *this;
    }

    __host__ __device__ const size_t (& getSize() const)[dim]
	{
    	return sz;
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
        return InvLinId(blockLinId, localLinId);
    }

    inline __host__ __device__ grid_key_dx<dim> InvLinId(mem_id blockLinId, mem_id localLinId) const
    {
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

    // Now methods to handle blockGrid coordinates (e.g. to load neighbouring blocks)
    template<typename indexT>
    inline __host__ __device__ mem_id BlockLinId(const grid_key_dx<dim, indexT> blockCoord) const
    {
        mem_id blockLinId = blockCoord.get(dim - 1);
        for (int d = dim - 2; d >= 0; --d)
        {
            blockLinId *= blockSz[d];
            blockLinId += blockCoord.get(d);
        }
        return blockLinId;
    }

    inline __host__ __device__ grid_key_dx<dim> BlockInvLinId(mem_id blockLinId) const
    {
        grid_key_dx<dim> blockCoord;
        for (int d = 0; d < dim; ++d)
        {
            auto c = blockLinId % blockSz[d];
            blockCoord.set_d(d, c);
            blockLinId /= blockSz[d];
        }
        return blockCoord;
    }

    inline size_t size() const
    {
    	size_t sz = 1;

    	for (size_t i = 0 ; i < dim ; i++)
    	{sz *= blockSz[i];}

    	return sz;
    }
};

#endif //OPENFPM_PDATA_BLOCKGEOMETRY_HPP
