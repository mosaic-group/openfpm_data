//
// Created by tommaso on 19/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH

#include <SparseGridGpu/BlockMapGpu.hpp>
#include "util/cuda_util.hpp"

// Kernels for SparseGridGpu

namespace SparseGridGpuKernels
{
    // This kernel is to be called with 1D parameters (?)
    template <unsigned int dim,
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

        int pos = blockIdx.x * blockDim.x + threadIdx.x;
        // NOTE: here we do 1 chunk per block! (we want to be sure to fit local memory constraints
        // since we will be loading also neighbouring elements!) (beware curse of dimensionality...)
        const unsigned int dataBlockPos = pos / blockSize;
        const unsigned int offset = pos % blockSize;

        const unsigned int enlargedBlockSize = sparseGrid.getEnlargedBlockSize();
        __shared__ MaskT enlargedBlock[enlargedBlockSize];

        if (dataBlockPos < indexBuffer.size())
        {
            const unsigned int dataBlockId = indexBuffer.template get<pIndex>(dataBlockPos);
            grid_key_dx<dim> dataBlockCoord = sparseGrid.getBlockCoord(dataBlockId);
            // Read block and ghost layer from global to shared memory
            sparseGrid.loadBlock(dataBlockCoord, enlargedBlock);
            sparseGrid.loadGhost(dataBlockCoord, enlargedBlock);
            //todo: Here code for tagging the boundary
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
//                        MaskT neighbourPlus =
                        //todo: finish this
                    }
                }
            }
            // Write block back to global memory
            sparseGrid.storeBlock(dataBlockCoord, enlargedBlock);
        }
        // todo:
        //#else // __NVCC__
        //        std::cout << __FILE__ << ":" << __LINE__ << " error: you are supposed to compile this file with nvcc, if you want to use it with gpu" << std::endl;
        //#endif // __NVCC__
    }
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
