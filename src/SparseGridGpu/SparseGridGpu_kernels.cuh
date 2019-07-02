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
            typename SparseGridT>
    __global__ void tagBoundaries(IndexBufT indexBuffer, DataBufT dataBuffer, SparseGridT sparseGrid)
    {
        //todo: #ifdef __NVCC__
        constexpr unsigned int pIndex = 0;
        
        typedef typename IndexBufT::value_type IndexAggregateT;
        typedef BlockTypeOf<IndexAggregateT , pIndex> IndexT;

        typedef typename DataBufT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;
        typedef ScalarTypeOf<AggregateT, pMask> MaskT;
        constexpr unsigned int blockSize = MaskBlockT::size;

        // NOTE: here we do 1 chunk per block! (we want to be sure to fit local memory constraints
        // since we will be loading also neighbouring elements!) (beware curse of dimensionality...)
//        int pos = blockIdx.x * blockDim.x + threadIdx.x;
//        const unsigned int dataBlockPos = pos / blockSize;
//        const unsigned int offset = pos % blockSize;
        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x;

        constexpr unsigned int enlargedBlockSize = IntPow<sparseGrid.getBlockEdgeSize() + 2*stencilSupportRadius, dim>::value;
        __shared__ MaskT enlargedBlock[enlargedBlockSize];

        if (dataBlockPos >= indexBuffer.size())
        {
            return;
        }
        const unsigned int dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);
        grid_key_dx<dim> dataBlockCoord = sparseGrid.getBlockCoord(dataBlockId);
        auto dataBlock = dataBuffer.get(dataBlockPos); // Avoid binary searches as much as possible
        // Read block and ghost layer from global to shared memory
        sparseGrid.loadBlock<pMask>(dataBlock, enlargedBlock);
        sparseGrid.loadGhost<pMask>(dataBlockCoord, enlargedBlock);
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
        // todo:
        //#else // __NVCC__
        //        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
        //#endif // __NVCC__
    }
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
