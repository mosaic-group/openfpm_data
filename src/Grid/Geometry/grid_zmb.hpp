/*
 * grid_zmb.hpp
 *
 *  Created on: Aug 1, 2019
 *      Author: i-bird
 */

#ifndef GRID_ZMB_HPP_
#define GRID_ZMB_HPP_



#include <boost/mpl/size_t.hpp>
#include <cstring>
#include <Grid/grid_sm.hpp>
#include "SparseGridGpu/TemplateUtils/mathUtils.hpp"
#include "util/zmorton.hpp"

/**
 * This class provides an interface to linearization of coordinates and viceversa when blocks are involved.
 * This can be seen as a lightweight version of grid_sm, with just LinId and InvLinId methods, but
 * tuned for blocked data.
 */
template<unsigned int dim, unsigned int blockEdgeSize, typename indexT>
class grid_zmb : private grid_smb<dim,blockEdgeSize,indexT>
{

public:

    grid_zmb()	{}

    __host__ __device__ grid_zmb(const size_t (& sz)[dim])
    :grid_smb<dim,blockEdgeSize,indexT>(sz)
    {}

    __host__ __device__ grid_zmb(const size_t domainBlockEdgeSize)
    :grid_smb<dim,blockEdgeSize,indexT>(domainBlockEdgeSize)
    {}

    template<typename T>
    __host__ __device__ grid_zmb(const grid_sm<dim, T> blockGrid)
    :grid_smb<dim,blockEdgeSize,indexT>(blockGrid)
    {}

#ifdef __NVCC__
    //Constructors from dim3 and uint3 objects
    __host__ __device__ grid_zmb(const dim3 blockDimensions)
    :grid_smb<dim,blockEdgeSize,indexT>(blockDimensions)
    {}


#endif // __NVCC__

    __host__ __device__ grid_zmb(const grid_zmb<dim, blockEdgeSize, indexT> &other)
    :grid_smb<dim,blockEdgeSize,indexT>(other)
    {}

    __host__ __device__ grid_zmb &operator=(const grid_zmb<dim, blockEdgeSize, indexT> &other)
    {
    	((grid_smb<dim,blockEdgeSize,indexT> *)this)->operator=(other);
        return *this;
    }

    template<typename indexT_>
    inline __host__ __device__ indexT LinId(const grid_key_dx<dim, indexT_> coord) const
    {
    	grid_key_dx<dim> key_b;
    	int localLinId = 0;
    	int sr = 1;

        for (int d = 0 ; d < dim; d++)
        {
            key_b.set_d(d,coord.get(d) / blockEdgeSize);
            localLinId += coord.get(d) % blockEdgeSize * sr;
            sr *= blockEdgeSize;
        }
        return lin_zid(key_b) * this->blockSize + localLinId;
    }

    inline __host__ __device__ grid_key_dx<dim, int> InvLinId(const indexT linId) const
    {
    	indexT linIdB = linId / this->blockSize;
    	int localLinId = linId % this->blockSize;

        return InvLinId(linIdB,localLinId);
    }

    inline __host__ __device__ grid_key_dx<dim, int> InvLinId(indexT blockLinId, indexT localLinId) const
    {
        grid_key_dx<dim,int> k;
        invlin_zid(blockLinId,k);

        for (indexT i = 0 ; i < dim ; i++)
        {
        	k.set_d(i,k.get(i)*blockEdgeSize + localLinId % blockEdgeSize);
        	localLinId /= blockEdgeSize;
        }
        return k;
    }

    // Now methods to handle blockGrid coordinates (e.g. to load neighbouring blocks)
    template<typename indexT_>
    inline __host__ __device__ indexT BlockLinId(const grid_key_dx<dim, indexT_> & blockCoord) const
    {
        return lin_zid(blockCoord);
    }

    inline __host__ __device__ grid_key_dx<dim, int> BlockInvLinId(indexT blockLinId) const
    {
    	grid_key_dx<dim,int> k;
    	invlin_zid(blockLinId,k);

        return k;
    }

    __host__ __device__ const indexT (& getSize() const)[dim]
	{
    	return grid_smb<dim,blockEdgeSize,indexT>::getSize();
	}


    // Now methods to handle blockGrid coordinates (e.g. to load neighbouring blocks)
    template<typename indexT_>
    inline __host__ __device__ grid_key_dx<dim,indexT> getGlobalCoord(const grid_key_dx<dim, indexT_> & blockCoord, unsigned int offset) const
    {
    	return grid_smb<dim,blockEdgeSize,indexT>::getGlobalCoord(blockCoord,offset);
    }

    inline indexT getBlockSize() const
    {
    	return grid_smb<dim,blockEdgeSize,indexT>::getBlockSize();
    }

    __host__ __device__ const indexT & size(int i) const
	{
    	return grid_smb<dim,blockEdgeSize,indexT>::size(i);
	}

    __host__ __device__ inline void swap(grid_zmb<dim, blockEdgeSize, indexT> &other)
    {
        grid_smb<dim,blockEdgeSize,indexT>::swap(other);
    }

    __host__ __device__ inline indexT getBlockEgdeSize() const
    {
        return blockEdgeSize;
    }

    __host__ __device__ inline indexT size() const
    {
        return grid_smb<dim,blockEdgeSize,indexT>::size();
    }
};


#endif /* GRID_ZMB_HPP_ */
