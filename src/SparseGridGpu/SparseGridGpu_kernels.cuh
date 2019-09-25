//
// Created by tommaso on 19/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH

#include <SparseGridGpu/BlockMapGpu.hpp>
#include <SparseGridGpu/TemplateUtils/mathUtils.hpp>
#include "util/cuda_util.hpp"

#ifndef SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE
#define SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE
#endif

#ifndef SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE_NO_SHARED
#define SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE_NO_SHARED
#endif

// Kernels for SparseGridGpu

namespace SparseGridGpuKernels
{
    // This kernel is to be called with 1D parameters (?)
    template <unsigned int dim,
            unsigned int stencilSupportRadius,
            unsigned int pMask,
            typename NN_type,
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

        const long long dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);
        auto dataBlock = dataBuffer.get(dataBlockPos); // Avoid binary searches as much as possible

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
                bool isPadding = NN_type::isPadding(sparseGrid,coord,enlargedBlock);
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
            typename nNN_type,
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

        const unsigned int dataBlockPos = pos / nNN_type::nNN;
        const unsigned int offset = pos % nNN_type::nNN;

        if (dataBlockPos >= indexBuffer.size())
        {return;}

        const auto dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);

        auto neighbourPos = sparseGrid.template getNeighboursPos<nNN_type>(dataBlockId, offset);

        nn_blocks.template get<0>(dataBlockPos*nNN_type::nNN + offset) = neighbourPos;
    }

    template <unsigned int dim,
            unsigned int pMask,
            typename stencil,
            typename IndexBufT,
            typename DataBufT,
            typename SparseGridT,
            typename... Args>
    __host__ void applyStencilInPlaceHost(
            IndexBufT &indexBuffer,
            DataBufT &dataBuffer,
            SparseGridT & sparseGrid,
            Args... args)
    {
        constexpr unsigned int pIndex = 0;

        typedef typename IndexBufT::value_type IndexAggregateT;
        typedef BlockTypeOf<IndexAggregateT , pIndex> IndexT;

        typedef typename DataBufT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;
        typedef ScalarTypeOf<AggregateT, pMask> MaskT;
        constexpr unsigned int blockSize = MaskBlockT::size;
        const auto bufferSize = indexBuffer.size();

        for (size_t blockId=0; blockId<bufferSize; ++blockId)
        {
            const unsigned int dataBlockPos = blockId;
            auto dataBlockLoad = dataBuffer.get(dataBlockPos); // Avoid binary searches as much as possible
            const unsigned int dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);

            for (size_t elementId=0; elementId<blockSize; ++elementId)
            {
                const unsigned int offset = elementId;

                // Read local mask to register
                const auto curMask = dataBlockLoad.template get<pMask>()[offset];
                grid_key_dx<dim, int> pointCoord = sparseGrid.getCoord(dataBlockId * blockSize + offset);

                bool applyStencilHere = true;

                if ((!sparseGrid.exist(curMask)) || sparseGrid.isPadding(curMask) || offset > blockSize)
                {
                    //            return; // We want to apply only on existing AND non-padding elements
                    applyStencilHere = false;
                }

                openfpm::sparse_index<unsigned int> sdataBlockPos;
                sdataBlockPos.id = dataBlockPos;

                stencil::stencilHost(
                        sparseGrid, dataBlockId, sdataBlockPos, offset, pointCoord, dataBlockLoad, dataBlockLoad,
                        applyStencilHere, args...);
            }
        }
    }

    template <unsigned int dim,
            unsigned int pMask,
            typename stencil,
            typename IndexBufT,
            typename DataBufT,
            typename SparseGridT,
            typename... Args>
    __global__ void
    SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE_NO_SHARED
    applyStencilInPlaceNoShared(
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

        int p = blockIdx.x * blockDim.x + threadIdx.x;

        auto & pntBuff = sparseGrid.getPointBuffer();

        if (p >= pntBuff.size())
        {
            return;
        }

        auto id = pntBuff.template get<0>(p);

        const unsigned int dataBlockPos = id / blockSize;
        const unsigned int offset = id % blockSize;

        auto dataBlockLoad = dataBuffer.get(dataBlockPos);

        const unsigned int dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);
        grid_key_dx<dim, int> pointCoord = sparseGrid.getCoord(dataBlockId * blockSize + offset);

        bool applyStencilHere = false;

        if (offset < blockSize)
        {
            // Read local mask to register
            const auto curMask = dataBlockLoad.template get<pMask>()[offset];
            applyStencilHere = sparseGrid.exist(curMask) && (!sparseGrid.isPadding(curMask));
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
    __global__ void
    SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE
    applyStencilInPlace(
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
        grid_key_dx<dim, int> pointCoord = sparseGrid.getCoord(dataBlockId * blockSize + offset);

        bool applyStencilHere = false;

        if (offset < blockSize)
        {
            // Read local mask to register
            const auto curMask = dataBlockLoad.template get<pMask>()[offset];
            applyStencilHere = sparseGrid.exist(curMask) && (!sparseGrid.isPadding(curMask));
        }

        openfpm::sparse_index<unsigned int> sdataBlockPos;
        sdataBlockPos.id = dataBlockPos;

        stencil::stencil(
                sparseGrid, dataBlockId, sdataBlockPos , offset, pointCoord, dataBlockLoad, dataBlockLoad,
                applyStencilHere, args...);
    }

    template<unsigned int pMask,
    		 typename dataBuffType,
    		 typename scanType,
    		 typename outType>
    __global__ void fill_e_points(dataBuffType dataBuf, scanType scanBuf, outType output)
    {
        typedef typename dataBuffType::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;
        constexpr unsigned int blockSize = MaskBlockT::size;

        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x % blockSize;

        __shared__ int ato_cnt;

        if (threadIdx.x == 0)
        {ato_cnt = 0;}

        __syncthreads();

        if (dataBlockPos >= scanBuf.size() - 1)
        {
            return;
        }

        int predicate = dataBuf.template get<pMask>(dataBlockPos)[offset] & 0x1;

        int id = atomicAdd(&ato_cnt,predicate);

        __syncthreads();

        if (predicate == true)
        {
        	output.template get<0>(id + scanBuf.template get<0>(dataBlockPos)) = offset + dataBlockPos * blockSize;
        }
    }

    template<unsigned int pMask,
    		 typename dataBufferType,
    		 typename outType>
    __global__ void calc_exist_points(dataBufferType dataBuf, outType output)
    {
    	typedef typename dataBufferType::value_type AggregateT;
    	typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;
    	constexpr unsigned int blockSize = MaskBlockT::size;

        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x % blockSize;

        __shared__ int ato_cnt;

        if (threadIdx.x == 0)
        {ato_cnt = 0;}

        __syncthreads();

        if (dataBlockPos >= output.size())
        {
            return;
        }

        int predicate = dataBuf.template get<pMask>(dataBlockPos)[offset] & 0x1;

        atomicAdd(&ato_cnt,predicate);

        __syncthreads();

        output.template get<0>(dataBlockPos) = ato_cnt;
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
        const unsigned int offset = threadIdx.x % blockSize;

        if (dataBlockPos >= indexBuffer.size())
        {
            return;
        }

        sparseGrid.init();
        __syncthreads();

        auto dataBlockLoad = dataBuffer.get(dataBlockPos); // Avoid binary searches as much as possible
        const auto dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);

        auto dataBlockStore = sparseGrid.insertBlock(dataBlockId);

        grid_key_dx<dim, int> pointCoord = sparseGrid.getCoord(dataBlockId * blockSize + offset);

        bool applyStencilHere = false;

        // Read local mask to register
        const auto curMask = dataBlockLoad.template get<pMask>()[offset];
        applyStencilHere = sparseGrid.exist(curMask) && (!sparseGrid.isPadding(curMask));
        if (applyStencilHere)
        {
            // Mark the current element in the new block as existing
            sparseGrid.setExist(dataBlockStore.template get<pMask>()[offset]);
        }

        openfpm::sparse_index<unsigned int> sdataBlockId;
        sdataBlockId.id = dataBlockPos;

        stencil::stencil(
                sparseGrid, dataBlockId, sdataBlockId, offset, pointCoord, dataBlockLoad, dataBlockStore,
                applyStencilHere, args...);

        __syncthreads();
        sparseGrid.flush_block_insert();
    }

    // Apply in-place operator on boundary
    template <unsigned int dim,
            unsigned int pMask,
            typename stencil,
            typename IndexBufT,
            typename DataBufT,
            typename SparseGridT,
            typename... Args>
    __global__ void applyBoundaryStencilInPlace(
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
        const unsigned int dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);
        grid_key_dx<dim, int> pointCoord = sparseGrid.getCoord(dataBlockId * blockSize + offset);

        bool applyStencilHere = false;

        if (offset < blockSize)
        {
            // Read local mask to register
            const auto curMask = dataBlockLoad.template get<pMask>()[offset];
            applyStencilHere = sparseGrid.isPadding(curMask) && sparseGrid.exist(curMask);
        }

        openfpm::sparse_index<unsigned int> sdataBlockPos;
        sdataBlockPos.id = dataBlockPos;

        stencil::stencil(
                sparseGrid, dataBlockId, sdataBlockPos , offset, pointCoord, dataBlockLoad, dataBlockLoad,
                applyStencilHere, args...);
    }
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
