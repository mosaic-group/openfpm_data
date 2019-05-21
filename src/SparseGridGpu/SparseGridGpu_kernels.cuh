//
// Created by tommaso on 16/05/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH

//#ifdef __NVCC__

namespace SparseGridGpuKernels
{
    template<unsigned int p, typename InsertBufferT, typename ScalarT>
    __global__ void initializeInsertBuffer(InsertBufferT insertBuffer, ScalarT backgroundValue)
    {
        typedef typename InsertBufferT::value_type AggregateT;
        typedef BlockTypeOf<AggregateT, p> BlockT;

        int pos = blockIdx.x * blockDim.x + threadIdx.x;
        unsigned int dataBlockId = pos / BlockT::size;
        unsigned int offset = pos % BlockT::size;

        insertBuffer.template get<p>(dataBlockId).existBitMask = 0;
        insertBuffer.template get<p>(dataBlockId)[offset] = backgroundValue;
    }
}
//#endif //__NVCC__

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
