//
// Created by tommaso on 19/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH

#include <SparseGridGpu/BlockMapGpu.hpp>
#include <SparseGridGpu/TemplateUtils/mathUtils.hpp>
#include "util/cuda_util.hpp"

// Kernels for SparseGridGpu

namespace SparseGridGpuKernels
{
    // This kernel is to be called with 1D parameters (?)
    template <unsigned int dim,
            unsigned int stencilSupportRadius,
            unsigned int pMask,
            typename IndexBufT,
            typename DataBufT,
            typename SparseGridT,
            typename nn_blocksT>
    __global__ void tagBoundaries(IndexBufT indexBuffer, DataBufT dataBuffer, SparseGridT sparseGrid,nn_blocksT nbT)
    {
        //todo: #ifdef __NVCC__
        constexpr unsigned int pIndex = 0;

        typedef typename IndexBufT::value_type IndexAggregateT;
        typedef BlockTypeOf<IndexAggregateT, pIndex> IndexT;

        typedef typename DataBufT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;
        typedef ScalarTypeOf<AggregateT, pMask> MaskT;
        constexpr unsigned int blockSize = MaskBlockT::size;

        // NOTE: here we do 1 chunk per block! (we want to be sure to fit local memory constraints
        // since we will be loading also neighbouring elements!) (beware curse of dimensionality...)
        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x;

        constexpr unsigned int enlargedBlockSize = IntPow<
                sparseGrid.getBlockEdgeSize() + 2 * stencilSupportRadius, dim>::value;
        __shared__ MaskT enlargedBlock[enlargedBlockSize];

        if (dataBlockPos >= indexBuffer.size())
        {
            return;
        }

        const auto dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);
        auto dataBlock = dataBuffer.get(dataBlockPos); // Avoid binary searches as much as possible

        if (dataBlockId < 0)
        {
        	printf("Negative Datablock \n");
        	return;
        }

        openfpm::sparse_index<unsigned int> sdataBlockPos;
        sdataBlockPos.id = dataBlockPos;
        sparseGrid.loadGhostBlock<pMask>(dataBlock,sdataBlockPos,enlargedBlock);

        __syncthreads();

        //Here code for tagging the boundary
        if (offset < blockSize)
        {
            const auto coord = sparseGrid.getCoordInEnlargedBlock(offset);
            const auto linId = sparseGrid.getLinIdInEnlargedBlock(offset);

            MaskT cur = enlargedBlock[linId];
            if (sparseGrid.exist(cur))
            {
                bool isPadding = false;
                for (int d=0; d<dim; ++d)
                {
                    auto nPlusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, 1);
                    auto nMinusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, -1);
                    MaskT neighbourPlus = enlargedBlock[nPlusId];
                    MaskT neighbourMinus = enlargedBlock[nMinusId];
                    isPadding = isPadding || (!sparseGrid.exist(neighbourPlus));
                    isPadding = isPadding || (!sparseGrid.exist(neighbourMinus));
                    if (isPadding) break;
                }
                if (isPadding)
                {
                    sparseGrid.setPadding(enlargedBlock[linId]);
                }
                else
                {
                    sparseGrid.unsetPadding(enlargedBlock[linId]);
                }
            }
        }
        // Write block back to global memory
        __syncthreads();
        sparseGrid.storeBlock<pMask>(dataBlock, enlargedBlock);
    }

    /*! \brief find the neighborhood of each chunk
     *
     * \param indexBuffer Chunk indec buffer
     * \param dataBuffer Output array of the neighborhood chunks
     * \param sparseGrid
     *
     */
    template <unsigned int dim,
            unsigned int nNN,
            typename IndexBufT,
            typename SparseGridT,
            typename nn_blocksT>
    __global__ void findNeighbours(IndexBufT indexBuffer, SparseGridT sparseGrid, nn_blocksT nn_blocks)
    {
        //todo: #ifdef __NVCC__
        constexpr unsigned int pIndex = 0;

        typedef typename IndexBufT::value_type IndexAggregateT;
        typedef BlockTypeOf<IndexAggregateT , pIndex> IndexT;

        const unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

        const unsigned int dataBlockPos = pos / nNN;
        const unsigned int offset = pos % nNN;

        if (dataBlockPos >= indexBuffer.size())
        {return;}

        const auto dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);

        auto neighbourPos = sparseGrid.getNeighboursPos(dataBlockId, offset);
        nn_blocks.template get<0>(dataBlockPos*nNN + offset) = neighbourPos;
    }

    template <unsigned int dim,
            unsigned int pMask,
            typename stencil,
            typename IndexBufT,
            typename DataBufT,
            typename SparseGridT,
            typename... Args>
    __global__ void applyStencilInPlace(
            IndexBufT indexBuffer,
            DataBufT dataBuffer,
            SparseGridT sparseGrid,
            Args... args)
    {
        constexpr unsigned int pIndex = 0;

        typedef typename IndexBufT::value_type IndexAggregateT;
        typedef BlockTypeOf<IndexAggregateT , pIndex> IndexT;

        typedef typename DataBufT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;
        typedef ScalarTypeOf<AggregateT, pMask> MaskT;
        constexpr unsigned int blockSize = MaskBlockT::size;

        // NOTE: here we do 1 chunk per block! (we want to be sure to fit local memory constraints
        // since we will be loading also neighbouring elements!) (beware curse of dimensionality...)
        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x;

        if (dataBlockPos >= indexBuffer.size())
        {
            return;
        }

        auto dataBlockLoad = dataBuffer.get(dataBlockPos); // Avoid binary searches as much as possible

        // todo: Add management of RED-BLACK stencil application! :)
        const unsigned int dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);
        // Read local mask to register
        const auto curMask = dataBlockLoad.template get<pMask>()[offset];
        grid_key_dx<dim, int> pointCoord = sparseGrid.getCoord(dataBlockId * blockSize + offset);

        bool applyStencilHere = true;

        if ( (!sparseGrid.exist(curMask)) || sparseGrid.isPadding(curMask) || offset > blockSize )
        {
            //            return; // We want to apply only on existing AND non-padding elements
            applyStencilHere = false;
        }

        openfpm::sparse_index<unsigned int> sdataBlockPos;
        sdataBlockPos.id = dataBlockPos;

        stencil::stencil(
                sparseGrid, dataBlockId, sdataBlockPos , offset, pointCoord, dataBlockLoad, dataBlockLoad,
                applyStencilHere, args...);
    }

    template <unsigned int dim,
            unsigned int pMask,
            typename stencil,
            typename IndexBufT,
            typename DataBufT,
            typename SparseGridT,
            typename... Args>
    __global__ void applyStencilInsert(
            IndexBufT indexBuffer,
            DataBufT dataBuffer,
            SparseGridT sparseGrid,
            Args... args)
    {
        //todo: #ifdef __NVCC__
        constexpr unsigned int pIndex = 0;

        typedef typename IndexBufT::value_type IndexAggregateT;
        typedef BlockTypeOf<IndexAggregateT , pIndex> IndexT;

        typedef typename DataBufT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;
        typedef ScalarTypeOf<AggregateT, pMask> MaskT;
        constexpr unsigned int blockSize = MaskBlockT::size;

        typedef decltype(sparseGrid.insertBlock(0U)) EncapT;

        // NOTE: here we do 1 chunk per block! (we want to be sure to fit local memory constraints
        // since we will be loading also neighbouring elements!) (beware curse of dimensionality...)
//        int pos = blockIdx.x * blockDim.x + threadIdx.x;
//        const unsigned int dataBlockPos = pos / blockSize;
//        const unsigned int offset = pos % blockSize;
        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x;

        if (dataBlockPos >= indexBuffer.size())
        {
            return;
        }

        sparseGrid.init();
        __syncthreads();

        auto dataBlockLoad = dataBuffer.get(dataBlockPos); // Avoid binary searches as much as possible
        const auto dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);

        auto dataBlockStore = sparseGrid.insertBlockNew(dataBlockId);

        // Read local mask to register
        const auto curMask = dataBlockLoad.template get<pMask>()[offset];

        grid_key_dx<dim, int> pointCoord = sparseGrid.getCoord(dataBlockId * blockSize + offset);


        bool applyStencilHere = true;

        if ( (!sparseGrid.exist(curMask)) || sparseGrid.isPadding(curMask) || offset > blockSize )
        {
            applyStencilHere = false;
        }

        openfpm::sparse_index<unsigned int> sdataBlockId;
        sdataBlockId.id = dataBlockPos;

        stencil::stencil(
                sparseGrid, dataBlockId, sdataBlockId, offset, pointCoord, dataBlockLoad, dataBlockStore,
                applyStencilHere, args...);
        sparseGrid.setExist(dataBlockStore.template get<pMask>()[offset]);

        __syncthreads();
        sparseGrid.flush_block_insert();
    }
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
