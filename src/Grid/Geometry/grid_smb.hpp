//
// Created by tommaso on 17/06/19.
//

#ifndef OPENFPM_PDATA_BLOCKGEOMETRY_HPP
#define OPENFPM_PDATA_BLOCKGEOMETRY_HPP

#include <boost/mpl/size_t.hpp>
#include <cstring>
#include <Grid/grid_sm.hpp>
#include "SparseGridGpu/TemplateUtils/mathUtils.hpp"

/**
 * This class provides an interface to linearization of coordinates and viceversa when blocks are involved.
 * This can be seen as a lightweight version of grid_sm, with just LinId and InvLinId methods, but
 * tuned for blocked data.
 */
template<unsigned int dim, unsigned int blockEdgeSize, typename indexT = long int>
class grid_smb
{
private:

    indexT blockSz[dim];
    indexT sz[dim];

protected:

    constexpr static indexT blockSize = IntPow<blockEdgeSize, dim>::value;

public:

    grid_smb()	{}

    __host__ __device__ grid_smb(const size_t (& sz)[dim])
    {
        for (int d=0; d<dim; ++d)
        {
            this->sz[d] = sz[d];
            blockSz[d] = sz[d] / blockEdgeSize + ((sz[d] % blockEdgeSize) != 0);
        }
    }

    __host__ __device__ grid_smb(const indexT (& sz)[dim])
    {
        for (int d=0; d<dim; ++d)
        {
            this->sz[d] = sz[d];
            blockSz[d] = sz[d] / blockEdgeSize + ((sz[d] % blockEdgeSize) != 0);
        }
    }

    __host__ __device__ grid_smb(const size_t domainBlockEdgeSize)
    {
        for (int i = 0; i < dim; ++i)
        {
            blockSz[i] = (indexT)domainBlockEdgeSize;
            sz[i] = (indexT)domainBlockEdgeSize * blockEdgeSize;
        }
    }

    template<typename T>
    __host__ __device__ grid_smb(const grid_sm<dim, T> blockGrid)
    {
        for (int i = 0; i < dim; ++i)
        {
            blockSz[i] = blockGrid.size(i);
            sz[i] = blockGrid.size(i) * blockEdgeSize;
        }
    }

#ifdef __NVCC__
    //Constructors from dim3 and uint3 objects
    __host__ __device__ grid_smb(const dim3 blockDimensions)
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


#endif // __NVCC__

    __host__ __device__ grid_smb(const grid_smb<dim, blockEdgeSize> &other)
    {
    	for (indexT i = 0 ; i < dim ; i++)
    	{
            blockSz[i] = other.blockSz[i];
            sz[i] = other.sz[i];
    	}
    }

    __host__ __device__ grid_smb &operator=(const grid_smb<dim, blockEdgeSize> &other)
    {
    	for (indexT i = 0 ; i < dim ; i++)
    	{
            blockSz[i] = other.blockSz[i];
            sz[i] = other.sz[i];
    	}
        return *this;
    }

    __host__ __device__ const indexT (& getSize() const)[dim]
	{
    	return sz;
	}

    /*! \brief Linearize the coordinate index
     *
     * The linearization is given by getting the block indexes and the local coordinate indexes
     *
     * Linearize the block index (blockLinId), linearize the local index (localLinId) and return
     * blockLinId * blockSize + localLinId
     *
     * \param coord coordinates
     *
     * \return linearized index
     *
     */
    template<typename indexT_>
    inline __host__ __device__ indexT LinId(const grid_key_dx<dim, indexT_> coord) const
    {
        //todo: Check (in debug mode only) that the coordinates passed here are valid and not overflowing dimensions (???)
        indexT blockLinId = coord.get(dim - 1) / blockEdgeSize;
        indexT localLinId = coord.get(dim - 1) % blockEdgeSize;
        for (int d = dim - 2; d >= 0; --d)
        {
            blockLinId *= blockSz[d];
            localLinId *= blockEdgeSize;
            blockLinId += coord.get(d) / blockEdgeSize;
            localLinId += coord.get(d) % blockEdgeSize;
        }
        return blockLinId * blockSize + localLinId;
    }

    inline __host__ __device__ grid_key_dx<dim, int> InvLinId(const indexT linId) const
    {
        indexT blockLinId = linId / blockSize;
        indexT localLinId = linId % blockSize;
        return InvLinId(blockLinId, localLinId);
    }

    /*! brief Given a local offset it provide the the position in coordinates
     *
     * \param linId point
     *
     * \return the point in coordinates
     *
     */
    inline __host__ __device__ grid_key_dx<dim, int> LocalInvLinId(unsigned int localLinId) const
    {
        grid_key_dx<dim, int> coord;
        for (int d = 0; d < dim; ++d)
        {
            auto c = localLinId % blockEdgeSize;
            coord.set_d(d, c);
            localLinId /= blockEdgeSize;
        }
        return coord;
    }

    /*! \brief Invert from the linearized block id + local id to the position of the point in coordinates
     *
     * From the point coordinated x,y,z you can get the block coordinated block_x,block_y,block_z
     *                and the local coordinates inside the chunk loc_x,loc_y,loc_z
     *
     * linearizing block coordinated you get blockId and linearizing the local coordinated you get localLinId
     *
     * each point in a sparse grid is identified by the formula blockId*blockSize + localLinId
     *
     * This function invert such formula.
     *
     * \param blockLinId is blockId
     *
     * \param localLinId is localLinId in the formula
     *
     *
     * \return the spatial coordinates of the point
     *
     */
    inline __host__ __device__ grid_key_dx<dim, int> InvLinId(indexT blockLinId, indexT localLinId) const
    {
        grid_key_dx<dim, int> coord;
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
    template<typename indexT_>
    inline __host__ __device__ indexT BlockLinId(const grid_key_dx<dim, indexT_> & blockCoord) const
    {
        indexT blockLinId = blockCoord.get(dim - 1);
        if (blockLinId >= blockSz[dim-1])
        {return -1;}

        for (int d = dim - 2; d >= 0; --d)
        {
            blockLinId *= blockSz[d];
            indexT cur = blockCoord.get(d);
            if (cur >= blockSz[d])
            {
                return -1;
            }
            blockLinId += cur;
        }
        return blockLinId;
    }

    // Now methods to handle blockGrid coordinates (e.g. to load neighbouring blocks)
    template<typename indexT_>
    inline __host__ __device__ grid_key_dx<dim,indexT> getGlobalCoord(const grid_key_dx<dim, indexT_> & blockCoord, unsigned int offset) const
    {
    	grid_key_dx<dim,indexT> k;

    	for (unsigned int i = 0 ; i < dim ; i++)
    	{
    		k.set_d(i,blockCoord.get(i)*blockEdgeSize + offset%blockEdgeSize);
    		offset /= blockEdgeSize;
    	}
    	return k;
    }

    inline __host__ __device__ grid_key_dx<dim, int> BlockInvLinId(indexT blockLinId) const
    {
        grid_key_dx<dim, int> blockCoord;
        for (int d = 0; d < dim; ++d)
        {
            auto c = blockLinId % blockSz[d];
            blockCoord.set_d(d, c);
            blockLinId /= blockSz[d];
        }
        return blockCoord;
    }

    inline indexT size() const
    {
    	indexT sz = 1;

    	for (indexT i = 0 ; i < dim ; i++)
    	{sz *= blockSz[i];}

    	return sz;
    }

    __host__ __device__ inline indexT getBlockSize() const
    {
    	return blockSize;
    }
};

#endif //OPENFPM_PDATA_BLOCKGEOMETRY_HPP
