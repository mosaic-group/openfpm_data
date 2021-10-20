//
// Created by tommaso on 19/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
#define OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH

#include <SparseGridGpu/BlockMapGpu.hpp>
#include <SparseGridGpu/TemplateUtils/mathUtils.hpp>
#include "util/cuda_util.hpp"
#include "SparseGrid/cp_block.hpp"

#ifndef SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE
#define SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE
#endif

#ifndef SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE_NO_SHARED
#define SPARSEGRIDGPU_LAUNCH_BOUND_APPLY_STENCIL_IN_PLACE_NO_SHARED
#endif

enum mask_sparse
{
	NOT_EXIST = 0,
	EXIST = 1,
	PADDING = 2,
	EXIST_AND_PADDING = 3
};

// Kernels for SparseGridGpu
namespace SparseGridGpuKernels
{
	template<typename SparseGridGpuType, typename pointers_type, typename headers_type>
	__global__ void unpack_headers(pointers_type pointers, headers_type headers, int * result, unsigned int sz_pack, int n_slot)
	{
		int t = threadIdx.x;

		if (t > pointers.size()) {return;}

		unsigned char * data_pack = (unsigned char *)pointers.template get<0>(t);

		while (data_pack < pointers.template get<1>(t) )
		{
			int ih = pointers.template get<2>(t);
			if (n_slot > ih)
			{
                if (sizeof(typename SparseGridGpuType::indexT_) == 8)
                {headers.template get<0>(t*n_slot + ih) = *(size_t *)data_pack;}
                else
                {
                    unsigned int dp1 = *(unsigned int *)data_pack;
                    unsigned int dp2 = *(unsigned int *)&(data_pack[4]);
                    headers.template get<0>(t*n_slot + ih) = (size_t)dp1 + (((size_t)dp2) << 32);
                }
				data_pack += sizeof(size_t);
				data_pack += SparseGridGpuType::unpack_headers(headers,data_pack,t*n_slot + ih,sz_pack);
				pointers.template get<2>(t) += 1;
			}
			else
			{
				// report error
				result[0] = 1;
				return;
			}			
		}
	}


	template<unsigned int dim>
	struct stencil_cross_func_impl
	{
		template<typename ScalarT, typename coordType, typename SparseGridT, unsigned int enlargedBlockSize, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil(ScalarT & res, ScalarT & cur, coordType & coord ,
				            ScalarT (& enlargedBlock)[enlargedBlockSize],
				            lambda_func f,
				            SparseGridT & sparseGrid, ArgsT ... args)
		{
			cross_stencil<dim,ScalarT> cs;

            for (int d = 0; d < dim; ++d)
            {
                auto nPlusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, 1);
                auto nMinusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, d, -1);
                ScalarT neighbourPlus = enlargedBlock[nPlusId];
                ScalarT neighbourMinus = enlargedBlock[nMinusId];

                cs.xm[d] = neighbourMinus;
                cs.xp[d] = neighbourPlus;
            }

            res = f(cur,cs, args ...);
		}
	};

	template<unsigned int dim>
	struct stencil_conv_func_impl
	{
		template<typename ScalarT, typename coordType, typename CpBlockType, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil(ScalarT & res, coordType & coord ,
				            CpBlockType & cpb,
				            lambda_func f,
				            ArgsT ... args)
		{
			printf("Convolution operation on GPU: Dimension not implemented \n");
		}

		template<typename ScalarT, typename coordType, typename CpBlockType, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil2(ScalarT & res1, ScalarT & res2, coordType & coord ,
				            CpBlockType & cpb1,
				            CpBlockType & cpb2,
				            lambda_func f,
				            ArgsT ... args)
		{
			printf("Convolution operation on GPU: Dimension not implemented \n");
		}
	};

	template<>
	struct stencil_conv_func_impl<3>
	{
		template<typename ScalarT, typename coordType, typename CpBlockType, typename DataBlockWrapperT, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil_block(ScalarT & res, coordType & coord ,
				            CpBlockType & cpb,
							DataBlockWrapperT & DataBlockLoad,
							int offset,
				            lambda_func f,
				            ArgsT ... args)
		{
			res = f(cpb,DataBlockLoad,offset,coord[0],coord[1],coord[2]);
		}

		template<typename ScalarT, typename coordType, typename CpBlockType, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil(ScalarT & res, coordType & coord ,
				            CpBlockType & cpb,
				            lambda_func f,
				            ArgsT ... args)
		{
			res = f(cpb,coord[0],coord[1],coord[2]);
		}

		template<typename ScalarT, typename coordType, typename CpBlockType, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil2(ScalarT & res1, ScalarT & res2, coordType & coord ,
				            CpBlockType & cpb1,
				            CpBlockType & cpb2,
				            lambda_func f,
				            ArgsT ... args)
		{
			f(res1,res2,cpb1,cpb2,coord[0],coord[1],coord[2]);
		}

		template<typename ScalarT, typename coordType, typename CpBlockType, typename DataBlockWrapperT, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil2_block(ScalarT & res1, ScalarT & res2, coordType & coord ,
				            CpBlockType & cpb1,
				            CpBlockType & cpb2,
							DataBlockWrapperT & DataBlockLoad,
							int offset,
				            lambda_func f,
				            ArgsT ... args)
		{
			f(res1,res2,cpb1,cpb2,DataBlockLoad,offset,coord[0],coord[1],coord[2]);
		}
	};

	template<>
	struct stencil_conv_func_impl<2>
	{
		template<typename ScalarT, typename coordType, typename CpBlockType, typename DataBlockWrapperT, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil_block(ScalarT & res, coordType & coord,
				            CpBlockType & cpb,
							DataBlockWrapperT & DataBlockLoad,
							int offset,
				            lambda_func f,
				            ArgsT ... args)
		{
			res = f(cpb,DataBlockLoad,offset,coord[0],coord[1]);
		}

		template<typename ScalarT, typename coordType, typename CpBlockType, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil(ScalarT & res, coordType & coord ,
				            CpBlockType & cpb,
				            lambda_func f,
				            ArgsT ... args)
		{
			res = f(cpb,coord[0],coord[1]);
		}

		template<typename ScalarT, typename coordType, typename CpBlockType, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil2(ScalarT & res1, ScalarT & res2, coordType & coord ,
				            CpBlockType & cpb1,
				            CpBlockType & cpb2,
				            lambda_func f,
				            ArgsT ... args)
		{
			f(res1,res2,cpb1,cpb2,coord[0],coord[1]);
		}

		template<typename ScalarT, typename coordType, typename CpBlockType, typename DataBlockWrapperT, typename lambda_func, typename ... ArgsT>
		__device__ static inline void stencil2_block(ScalarT & res1, ScalarT & res2, coordType & coord ,
				            CpBlockType & cpb1,
				            CpBlockType & cpb2,
							DataBlockWrapperT & DataBlockLoad,
							int offset,
				            lambda_func f,
				            ArgsT ... args)
		{
			f(res1,res2,cpb1,cpb2,DataBlockLoad,offset,coord[0],coord[1]);
		}
	};

	template<unsigned int dim, unsigned int n_loop, unsigned int p_src, unsigned int p_dst, unsigned int stencil_size>
	struct stencil_cross_func_conv
	{
		typedef NNStar<dim> stencil_type;

		static constexpr unsigned int supportRadius = stencil_size;

		template<typename SparseGridT, typename DataBlockWrapperT, typename lambda_func, typename ... ArgT>
		static inline __device__ void stencil(
				SparseGridT & sparseGrid,
				const unsigned int dataBlockId,
				openfpm::sparse_index<unsigned int> dataBlockIdPos,
				unsigned int offset,
				grid_key_dx<dim, int> & pointCoord,
				DataBlockWrapperT & dataBlockLoad,
				DataBlockWrapperT & dataBlockStore,
				unsigned char curMask,
				lambda_func f,
				ArgT ... args)
		{
	        typedef typename SparseGridT::AggregateBlockType AggregateT;
	        typedef ScalarTypeOf<AggregateT, p_src> ScalarT;

	        constexpr unsigned int enlargedBlockSize = IntPow<
	                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

	        __shared__ ScalarT enlargedBlock[enlargedBlockSize];

	        for (int i = 0; i < n_loop ; i++)
	        {
	        	if (i*IntPow<SparseGridT::getBlockEdgeSize(), dim>::value + threadIdx.x < enlargedBlockSize)
	        	{
	        		enlargedBlock[i*IntPow<SparseGridT::getBlockEdgeSize(), dim>::value + threadIdx.x] = sparseGrid.getblockMap().template getBackground<p_src>()[0];
	        	}
	        }

	        __syncthreads();

	        typedef typename vmpl_create_constant<dim,SparseGridT::blockEdgeSize_>::type block_sizes;
	        typedef typename vmpl_sum_constant<2*stencil_size,block_sizes>::type vmpl_sizes;

	        cp_block<ScalarT,stencil_size,vmpl_sizes,dim> cpb(enlargedBlock);

	        sparseGrid.template loadGhostBlock<p_src>(dataBlockLoad, dataBlockIdPos, enlargedBlock);

	        __syncthreads();

	        ScalarT res = 0;

	        if ((curMask & mask_sparse::EXIST) && !(curMask & mask_sparse::PADDING))
	        {
	        	int coord[dim];

				unsigned int linIdTmp = offset;
				for (unsigned int d = 0; d < dim; ++d)
				{
					coord[d] = linIdTmp % SparseGridT::blockEdgeSize_;
					linIdTmp /= SparseGridT::blockEdgeSize_;
				}

	            stencil_conv_func_impl<dim>::stencil(res,coord,cpb,f,args...);

	            dataBlockStore.template get<p_dst>()[offset] = res;
	        }
		}

	    template <typename SparseGridT, typename CtxT>
	    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
	    {
	        // No flush
	    }
	};


	template<unsigned int dim, unsigned int n_loop, unsigned int p_src, unsigned int p_dst, unsigned int stencil_size>
	struct stencil_cross_func_conv_block_read
	{
		typedef NNStar<dim> stencil_type;

		static constexpr unsigned int supportRadius = stencil_size;

		template<typename SparseGridT, typename DataBlockWrapperT, typename lambda_func, typename ... ArgT>
		static inline __device__ void stencil(
				SparseGridT & sparseGrid,
				const unsigned int dataBlockId,
				openfpm::sparse_index<unsigned int> dataBlockIdPos,
				unsigned int offset,
				grid_key_dx<dim, int> & pointCoord,
				DataBlockWrapperT & dataBlockLoad,
				DataBlockWrapperT & dataBlockStore,
				unsigned char curMask,
				lambda_func f,
				ArgT ... args)
		{
	        typedef typename SparseGridT::AggregateBlockType AggregateT;
	        typedef ScalarTypeOf<AggregateT, p_src> ScalarT;

	        constexpr unsigned int enlargedBlockSize = IntPow<
	                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

	        __shared__ ScalarT enlargedBlock[enlargedBlockSize];

	        for (int i = 0; i < n_loop ; i++)
	        {
	        	if (i*IntPow<SparseGridT::getBlockEdgeSize(), dim>::value + threadIdx.x < enlargedBlockSize)
	        	{
	        		enlargedBlock[i*IntPow<SparseGridT::getBlockEdgeSize(), dim>::value + threadIdx.x] = sparseGrid.getblockMap().template getBackground<p_src>()[0];
	        	}
	        }

	        __syncthreads();

	        typedef typename vmpl_create_constant<dim,SparseGridT::blockEdgeSize_>::type block_sizes;
	        typedef typename vmpl_sum_constant<2*stencil_size,block_sizes>::type vmpl_sizes;

	        cp_block<ScalarT,stencil_size,vmpl_sizes,dim> cpb(enlargedBlock);

	        sparseGrid.template loadGhostBlock<p_src>(dataBlockLoad, dataBlockIdPos, enlargedBlock);

	        __syncthreads();

	        ScalarT res = 0;

	        if ((curMask & mask_sparse::EXIST) && !(curMask & mask_sparse::PADDING))
	        {
	        	int coord[dim];

				unsigned int linIdTmp = offset;
				for (unsigned int d = 0; d < dim; ++d)
				{
					coord[d] = linIdTmp % SparseGridT::blockEdgeSize_;
					linIdTmp /= SparseGridT::blockEdgeSize_;
				}

	            stencil_conv_func_impl<dim>::stencil_block(res,coord,cpb,dataBlockLoad,offset,f,args...);

	            dataBlockStore.template get<p_dst>()[offset] = res;
	        }
		}

	    template <typename SparseGridT, typename CtxT>
	    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
	    {
	        // No flush
	    }
	};

	template<unsigned int dim, unsigned int n_loop, unsigned int p_src1, unsigned int p_src2, unsigned int p_dst1, unsigned int p_dst2, unsigned int stencil_size>
	struct stencil_func_conv2_b
	{
		typedef NNStar<dim> stencil_type;

		static constexpr unsigned int supportRadius = stencil_size;

		template<typename SparseGridT, typename DataBlockWrapperT, typename lambda_func, typename ... ArgT>
		static inline __device__ void stencil(
				SparseGridT & sparseGrid,
				const unsigned int dataBlockId,
				openfpm::sparse_index<unsigned int> dataBlockIdPos,
				unsigned int offset,
				grid_key_dx<dim, int> & pointCoord,
				DataBlockWrapperT & dataBlockLoad,
				DataBlockWrapperT & dataBlockStore,
				unsigned char curMask,
				lambda_func f,
				ArgT ... args)
		{
	        typedef typename SparseGridT::AggregateBlockType AggregateT;
	        typedef ScalarTypeOf<AggregateT, p_src1> ScalarT1;
	        typedef ScalarTypeOf<AggregateT, p_src1> ScalarT2;

	        constexpr unsigned int enlargedBlockSize = IntPow<
	                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

	        __shared__ ScalarT1 enlargedBlock1[enlargedBlockSize];
	        __shared__ ScalarT2 enlargedBlock2[enlargedBlockSize];

	        // fill with background

	        typedef typename vmpl_create_constant<dim,SparseGridT::blockEdgeSize_>::type block_sizes;
	        typedef typename vmpl_sum_constant<2*stencil_size,block_sizes>::type vmpl_sizes;

	        cp_block<ScalarT1,stencil_size,vmpl_sizes,dim> cpb1(enlargedBlock1);
	        cp_block<ScalarT2,stencil_size,vmpl_sizes,dim> cpb2(enlargedBlock2);

	        sparseGrid.template loadGhostBlock<p_src1>(dataBlockLoad, dataBlockIdPos, enlargedBlock1);
	        sparseGrid.template loadGhostBlock<p_src2>(dataBlockLoad, dataBlockIdPos, enlargedBlock2);

	        __syncthreads();

	        ScalarT1 res1 = 0;
	        ScalarT2 res2 = 0;

	        if ((curMask & mask_sparse::EXIST) && !(curMask & mask_sparse::PADDING))
	        {
	        	int coord[dim];

				unsigned int linIdTmp = offset;
				for (unsigned int d = 0; d < dim; ++d)
				{
					coord[d] = linIdTmp % SparseGridT::blockEdgeSize_;
					linIdTmp /= SparseGridT::blockEdgeSize_;
				}

	            stencil_conv_func_impl<dim>::stencil2_block(res1,res2,coord,cpb1,cpb2,dataBlockLoad,offset,f,args...);

	            dataBlockStore.template get<p_dst1>()[offset] = res1;
	            dataBlockStore.template get<p_dst2>()[offset] = res2;
	        }
		}

	    template <typename SparseGridT, typename CtxT>
	    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
	    {
	        // No flush
	    }
	};

	template<unsigned int dim, unsigned int n_loop, unsigned int p_src1, unsigned int p_src2, unsigned int p_dst1, unsigned int p_dst2, unsigned int stencil_size>
	struct stencil_func_conv2
	{
		typedef NNStar<dim> stencil_type;

		static constexpr unsigned int supportRadius = stencil_size;

		template<typename SparseGridT, typename DataBlockWrapperT, typename lambda_func, typename ... ArgT>
		static inline __device__ void stencil(
				SparseGridT & sparseGrid,
				const unsigned int dataBlockId,
				openfpm::sparse_index<unsigned int> dataBlockIdPos,
				unsigned int offset,
				grid_key_dx<dim, int> & pointCoord,
				DataBlockWrapperT & dataBlockLoad,
				DataBlockWrapperT & dataBlockStore,
				unsigned char curMask,
				lambda_func f,
				ArgT ... args)
		{
	        typedef typename SparseGridT::AggregateBlockType AggregateT;
	        typedef ScalarTypeOf<AggregateT, p_src1> ScalarT1;
	        typedef ScalarTypeOf<AggregateT, p_src1> ScalarT2;

	        constexpr unsigned int enlargedBlockSize = IntPow<
	                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

	        __shared__ ScalarT1 enlargedBlock1[enlargedBlockSize];
	        __shared__ ScalarT2 enlargedBlock2[enlargedBlockSize];

	        // fill with background

	        typedef typename vmpl_create_constant<dim,SparseGridT::blockEdgeSize_>::type block_sizes;
	        typedef typename vmpl_sum_constant<2*stencil_size,block_sizes>::type vmpl_sizes;

	        cp_block<ScalarT1,stencil_size,vmpl_sizes,dim> cpb1(enlargedBlock1);
	        cp_block<ScalarT2,stencil_size,vmpl_sizes,dim> cpb2(enlargedBlock2);

	        sparseGrid.template loadGhostBlock<p_src1>(dataBlockLoad, dataBlockIdPos, enlargedBlock1);
	        sparseGrid.template loadGhostBlock<p_src2>(dataBlockLoad, dataBlockIdPos, enlargedBlock2);

	        __syncthreads();

	        ScalarT1 res1 = 0;
	        ScalarT2 res2 = 0;

	        if ((curMask & mask_sparse::EXIST) && !(curMask & mask_sparse::PADDING))
	        {
	        	int coord[dim];

				unsigned int linIdTmp = offset;
				for (unsigned int d = 0; d < dim; ++d)
				{
					coord[d] = linIdTmp % SparseGridT::blockEdgeSize_;
					linIdTmp /= SparseGridT::blockEdgeSize_;
				}

	            stencil_conv_func_impl<dim>::stencil2(res1,res2,coord,cpb1,cpb2,f,args...);

	            dataBlockStore.template get<p_dst1>()[offset] = res1;
	            dataBlockStore.template get<p_dst2>()[offset] = res2;
	        }
		}

	    template <typename SparseGridT, typename CtxT>
	    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
	    {
	        // No flush
	    }
	};

	template<unsigned int dim, unsigned int p_src, unsigned int p_dst, unsigned int stencil_size>
	struct stencil_cross_func
	{
		typedef NNStar<dim> stencil_type;

		static constexpr unsigned int supportRadius = stencil_size;

		template<typename SparseGridT, typename DataBlockWrapperT, typename lambda_func, typename ... ArgT>
		static inline __device__ void stencil(
				SparseGridT & sparseGrid,
				const unsigned int dataBlockId,
				openfpm::sparse_index<unsigned int> dataBlockIdPos,
				unsigned int offset,
				grid_key_dx<dim, int> & pointCoord,
				DataBlockWrapperT & dataBlockLoad,
				DataBlockWrapperT & dataBlockStore,
				unsigned char curMask,
				lambda_func f,
				ArgT ... args)
		{
	        typedef typename SparseGridT::AggregateBlockType AggregateT;
	        typedef ScalarTypeOf<AggregateT, p_src> ScalarT;

	        constexpr unsigned int enlargedBlockSize = IntPow<
	                SparseGridT::getBlockEdgeSize() + 2 * supportRadius, dim>::value;

	        __shared__ ScalarT enlargedBlock[enlargedBlockSize];

	        sparseGrid.template loadGhostBlock<p_src>(dataBlockLoad, dataBlockIdPos, enlargedBlock);

	        __syncthreads();

	        decltype(sparseGrid.getLinIdInEnlargedBlock(0)) linId = 0;
	        ScalarT res = 0;

	        if ((curMask & mask_sparse::EXIST) && !(curMask & mask_sparse::PADDING))
	        {
	            const auto coord = sparseGrid.getCoordInEnlargedBlock(offset);
	            // const auto linId = sparseGrid.getLinIdInEnlargedBlock(offset);
	            linId = sparseGrid.getLinIdInEnlargedBlock(offset);
	            ScalarT cur = enlargedBlock[linId];

	            stencil_cross_func_impl<dim>::stencil(res,cur,coord,enlargedBlock,f,sparseGrid,args ...);


	        }
	        __syncthreads();
	        if ((curMask & mask_sparse::EXIST) && !(curMask & mask_sparse::PADDING))
	        {
	            enlargedBlock[linId] = res;
	        }
	        __syncthreads();
	        sparseGrid.template storeBlock<p_dst>(dataBlockStore, enlargedBlock);
		}

	    template <typename SparseGridT, typename CtxT>
	    static inline void __host__ flush(SparseGridT & sparseGrid, CtxT & ctx)
	    {
	        // No flush
	    }
	};

    // This kernel is to be called with 1D parameters (?)
    template <unsigned int dim,
            unsigned int stencilSupportRadius,
            unsigned int pMask,
            typename NN_type,
            typename checker_type,
            typename IndexBufT,
            typename DataBufT,
            typename SparseGridT,
            typename nn_blocksT>
    __global__ void tagBoundaries(IndexBufT indexBuffer, DataBufT dataBuffer, SparseGridT sparseGrid,nn_blocksT nbT, checker_type chk)
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
        sparseGrid.template loadGhostBlock<pMask>(dataBlock,sdataBlockPos,enlargedBlock);

        __syncthreads();

        bool check = chk.check(sparseGrid,dataBlockId,offset);

        //Here code for tagging the boundary
        if (offset < blockSize && check == true)
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
        sparseGrid.template storeBlock<pMask>(dataBlock, enlargedBlock);
    }

    /*! \brief construct the link between 2 sparse grid
     *
     *
     */
    template<unsigned int dim, unsigned int pMask, unsigned int chunk_size , typename SparseGridType, typename outputType>
    __global__ void link_construct(SparseGridType grid_up, SparseGridType grid_cu, outputType out)
    {
        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x;

        auto & indexBuffer = grid_cu.getIndexBuffer();
        auto & dataBuffer = grid_cu.getDataBuffer();

        // if the point is a padding
        if (dataBuffer.template get <pMask>(dataBlockPos)[offset] & 0x2)
        {
        	auto id = indexBuffer.template get<0>(dataBlockPos);
        	grid_key_dx<dim,int> pos = grid_cu.getCoord(id*chunk_size + offset);

        	printf("HERE %d %d \n",pos.get(0),pos.get(1));

        	for (int i = 0 ; i < dim ; i++)
        	{pos.set_d(i,pos.get(i) / 2);}

        	if (grid_up.template get<pMask>(pos) == 0x1)
        	{
        		atomicAdd(&out.template get<0>(dataBlockPos),1);
        	}
        }
    }

    /*! \brief count the padding particles
     *
     *
     */
    template<unsigned int dim, unsigned int pMask, unsigned int chunk_size , typename SparseGridType, typename outputType, typename BoxType>
    __global__ void count_paddings(SparseGridType grid_cu, outputType out, BoxType box)
    {
        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x;

        auto & indexBuffer = grid_cu.getIndexBuffer();
        auto & dataBuffer = grid_cu.getDataBuffer();

        auto id = indexBuffer.template get<0>(dataBlockPos);

        // check if the point is inside the box
        auto coord = grid_cu.getCoord(id,offset);

        bool active = box.isInsideKey(coord);

        if (active == false)
        {return;}

        // if the point is a padding
        if (dataBuffer.template get <pMask>(dataBlockPos)[offset] & 0x2)
        {
        	atomicAdd(&out.template get<0>(dataBlockPos),1);
        }
    }

    /*! \brief count the padding particles
     *
     *
     */
    template<unsigned int pMask, typename SparseGridType, typename ScanType, typename outputType, typename BoxType>
    __global__ void collect_paddings(SparseGridType grid_cu, ScanType stp, outputType out, BoxType box)
    {
        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x;

        __shared__ int counter;
        counter = 0;
        __syncthreads();

        auto & indexBuffer = grid_cu.getIndexBuffer();
        auto & dataBuffer = grid_cu.getDataBuffer();

        auto id = indexBuffer.template get<0>(dataBlockPos);

        // check if the point is inside the box
        auto coord = grid_cu.getCoord(id,offset);

        bool active = box.isInsideKey(coord);

        if (active == false)
        {return;}

        int pad_offset = stp.template get<0>(dataBlockPos);

        // if the point is a padding
        if (dataBuffer.template get <pMask>(dataBlockPos)[offset] & 0x2)
        {
        	int cnt = atomicAdd(&counter,1);

        	out.template get<0>(pad_offset + cnt) = dataBlockPos;
        	out.template get<1>(pad_offset + cnt) = offset;
        }
    }

    /*! \brief construct the link between 2 sparse grid
     *
     *
     */
    template<unsigned int dim, unsigned int pMask, unsigned int chunk_size,
    		 typename padPointType , typename SparseGridType,
    		 typename outputType>
    __global__ void link_construct_dw_count(padPointType padPoints, SparseGridType grid_dw, SparseGridType grid_cu, outputType out, Point<dim,int> p_dw)
    {
    	const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

    	if (p >= padPoints.size())	{return;}

        const unsigned int dataBlockPos = padPoints.template get<0>(p);
        const unsigned int offset = padPoints.template get<1>(p);

        auto & indexBuffer = grid_cu.getIndexBuffer();
        auto & dataBuffer = grid_cu.getDataBuffer();

		auto id = indexBuffer.template get<0>(dataBlockPos);
		grid_key_dx<dim,int> pos = grid_cu.getCoord(id*chunk_size + offset);

		for (int i = 0 ; i < dim ; i++)
		{pos.set_d(i,pos.get(i) * 2 + p_dw.get(i) );}

		for (int j = 0 ; j < 2*dim ; j++)
		{
			grid_key_dx<dim,int> kc;
			for (int k = 0 ; k < dim ; k++)
			{
				kc.set_d(k,pos.get(k) + ((j >> k) & 0x1) );
			}

			if (grid_dw.template get<pMask>(kc) & 0x1)
			{
				int a = atomicAdd(&out.template get<0>(p),1);
			}
		}
    }

    /*! \brief construct the link between 2 sparse grid
     *
     *
     */
    template<unsigned int dim, unsigned int pMask, unsigned int chunk_size,
    		 typename padPointType , typename SparseGridType,
    		 typename outputType>
    __global__ void link_construct_up_count(padPointType padPoints, SparseGridType grid_up, SparseGridType grid_cu, outputType out, Point<dim,int> p_up)
    {
    	const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

    	if (p >= padPoints.size())	{return;}

        const unsigned int dataBlockPos = padPoints.template get<0>(p);
        const unsigned int offset = padPoints.template get<1>(p);

        auto & indexBuffer = grid_cu.getIndexBuffer();
        auto & dataBuffer = grid_cu.getDataBuffer();

		auto id = indexBuffer.template get<0>(dataBlockPos);
		grid_key_dx<dim,int> pos = grid_cu.getCoord(id*chunk_size + offset);

		for (int i = 0 ; i < dim ; i++)
		{pos.set_d(i,(pos.get(i) - p_up.get(i)) / 2);}

		if (grid_up.template get<pMask>(pos) & 0x1)
		{
			int a = atomicAdd(&out.template get<0>(p),1);
		}
    }

    /*! \brief construct the link between 2 sparse grid
     *
     *
     */
    template<unsigned int dim, unsigned int pMask, unsigned int chunk_size,
    typename padPointType , typename SparseGridType, typename scanType, typename outputType>
    __global__ void link_construct_insert_dw(padPointType padPoints, SparseGridType grid_dw, SparseGridType grid_cu, scanType scan, outputType out, Point<dim,int> p_dw)
    {
    	const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

    	if (p >= padPoints.size())	{return;}

        const unsigned int dataBlockPos = padPoints.template get<0>(p);
        const unsigned int offset = padPoints.template get<1>(p);

        auto & indexBuffer = grid_cu.getIndexBuffer();
        auto & dataBuffer = grid_cu.getDataBuffer();

        auto & dataBuffer_dw = grid_dw.getDataBuffer();

		auto id = indexBuffer.template get<0>(dataBlockPos);
		grid_key_dx<dim,int> pos = grid_cu.getCoord(id*chunk_size + offset);

		for (int i = 0 ; i < dim ; i++)
		{pos.set_d(i,pos.get(i) * 2 + p_dw.get(i) );}

		unsigned int dataBlockPos_dw;
		unsigned int offset_dw;

		int link_offset = scan.template get<0>(p);

		int c = 0;
		for (int j = 0 ; j < 2*dim ; j++)
		{
			grid_key_dx<dim,int> kc;
			for (int k = 0 ; k < dim ; k++)
			{
				kc.set_d(k,pos.get(k) + ((j >> k) & 0x1) );
			}

			grid_dw.get_sparse(kc,dataBlockPos_dw,offset_dw);

			if (dataBuffer_dw.template get<pMask>(dataBlockPos_dw)[offset_dw] & 0x1)
			{
				out.template get<0>(link_offset + c) = dataBlockPos_dw;
				out.template get<1>(link_offset + c) = offset_dw;

				c++;
			}
		}
    }

    /*! \brief construct the link between 2 sparse grid
     *
     *
     */
    template<unsigned int dim, unsigned int pMask, unsigned int chunk_size,
    typename padPointType , typename SparseGridType, typename scanType, typename outputType>
    __global__ void link_construct_insert_up(padPointType padPoints, SparseGridType grid_up, SparseGridType grid_cu, scanType scan, outputType out, Point<dim,int> p_up)
    {
    	const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

    	if (p >= padPoints.size())	{return;}

        const unsigned int dataBlockPos = padPoints.template get<0>(p);
        const unsigned int offset = padPoints.template get<1>(p);

        auto & indexBuffer = grid_cu.getIndexBuffer();
        auto & dataBuffer = grid_cu.getDataBuffer();

        auto & dataBuffer_dw = grid_up.getDataBuffer();

		auto id = indexBuffer.template get<0>(dataBlockPos);
		grid_key_dx<dim,int> pos = grid_cu.getCoord(id*chunk_size + offset);

		for (int i = 0 ; i < dim ; i++)
		{pos.set_d(i,(pos.get(i) - p_up.get(i)) / 2);}

		unsigned int dataBlockPos_dw;
		unsigned int offset_dw;

		int link_offset = scan.template get<0>(p);

		grid_up.get_sparse(pos,dataBlockPos_dw,offset_dw);

		if (dataBuffer_dw.template get<pMask>(dataBlockPos_dw)[offset_dw] & 0x1)
		{
			out.template get<0>(link_offset) = dataBlockPos_dw;
			out.template get<1>(link_offset) = offset_dw;
		}
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
    __global__ void
    applyStencilInPlaceNoShared(
    		Box<dim,int> bx,
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

        unsigned char curMask;

        if (offset < blockSize)
        {
            // Read local mask to register
            curMask = dataBlockLoad.template get<pMask>()[offset];
            if (bx.isInsideKey(pointCoord) == false)	{curMask = 0;}
        }

        openfpm::sparse_index<unsigned int> sdataBlockPos;
        sdataBlockPos.id = dataBlockPos;

        stencil::stencil(
                sparseGrid, dataBlockId, sdataBlockPos , offset, pointCoord, dataBlockLoad, dataBlockLoad,
                curMask, args...);
    }

    template <unsigned int dim,
            unsigned int pMask,
            typename stencil,
            typename IndexBufT,
            typename DataBufT,
            typename SparseGridT,
            typename... Args>
    __global__ void
    applyStencilInPlace(
    		Box<dim,int> bx,
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

        unsigned char curMask;

        if (offset < blockSize)
        {
            // Read local mask to register
            curMask = dataBlockLoad.template get<pMask>()[offset];
			for (int i = 0 ; i < dim ; i++)
			{curMask &= (pointCoord.get(i) < bx.getLow(i) || pointCoord.get(i) > bx.getHigh(i))?0:0xFF;}
        }

        openfpm::sparse_index<unsigned int> sdataBlockPos;
        sdataBlockPos.id = dataBlockPos;

        stencil::stencil(
                sparseGrid, dataBlockId, sdataBlockPos , offset, pointCoord, dataBlockLoad, dataBlockLoad,
                curMask, args...);
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

    template<unsigned int dim,
    		 unsigned int pMask,
    		 unsigned int blockEdgeSize,
    		 typename outType,
    		 typename boxesVector_type,
    		 typename grid_smb_type,
    		 typename SparseGrid_type>
    __global__ void calc_remove_points_chunks_boxes(SparseGrid_type sparseGrid,
    											 boxesVector_type boxes,
    											 grid_smb_type grd,
    											 outType output)
    {
        const unsigned int dataBlockPos = blockIdx.x * blockDim.x + threadIdx.x;

        if (dataBlockPos >= sparseGrid.numBlocks())
        {return;}

        auto id = sparseGrid.getIndex(dataBlockPos);
        grid_key_dx<dim,int> pnt = grd.InvLinId(id*grd.getBlockSize());

        Box<dim,unsigned int> b;

        for (int i = 0 ; i < dim ; i++)
        {
        	b.setLow(i,pnt.get(i));
        	b.setHigh(i,pnt.get(i) + blockEdgeSize - 1);
        }

        // this block intersect a remove box section so mark the chunk

        output.template get<1>(dataBlockPos) = 0;
		for (int k = 0 ; k < boxes.size() ; k++ )
		{
			Box<dim,unsigned int> btest = boxes.get(k);

			Box<dim,unsigned int> bout;

			if (btest.Intersect(b,bout) == true)
			{
				output.template get<1>(dataBlockPos) = 1;
			}
		}
    }

    template<typename outType,
    		 typename activeCnkType>
    __global__ void collect_rem_chunks(activeCnkType act,
    								   outType output)
    {
        const unsigned int dataBlockPos = blockIdx.x * blockDim.x + threadIdx.x;

        if (dataBlockPos >= act.size()-1)
        {return;}

        auto id = act.template get<1>(dataBlockPos);
        auto id_p1 = act.template get<1>(dataBlockPos+1);

        if (id != id_p1)
        {
        	output.template get<0>(id) = dataBlockPos;
        }
    }

    template<unsigned int dim, unsigned int pMask,
    						    typename SparseGridType,
    						    typename grid_smb_type,
    						    typename activeCntType,
    						    typename boxesType>
    __global__ void remove_points(SparseGridType sparseGrid,
    							  grid_smb_type grd,
    							  activeCntType active_blocks,
    							  boxesType boxes)
    {
    	constexpr unsigned int blockSize = sparseGrid.getBlockSize();

        const unsigned int dataBlockPos = active_blocks.template get<0>(blockIdx.x);
        const unsigned int offset = threadIdx.x % blockSize;

        if (dataBlockPos >= sparseGrid.numBlocks())
        {return;}

        int predicate = sparseGrid.getMask(dataBlockPos,offset) & 0x1;
        // calculate point coord;
        auto id = sparseGrid.getIndex(dataBlockPos);
        grid_key_dx<dim,int> pnt = grd.InvLinId(id*grd.getBlockSize() + offset);
        Point<dim,unsigned int> p;

        for (int i = 0 ; i < dim ; i++)
        {p.get(i) = pnt.get(i);}

        // if this block intersect any box

        if (predicate == true)
        {
			for (int k = 0 ; k < boxes.size() ; k++ )
			{
				Box<dim,unsigned int> box = boxes.get(k);

				if (box.isInside(p) == true)
				{
					sparseGrid.getMask(dataBlockPos,offset) = 0;
				}
			}
        }
    }

    template<unsigned int dim,
    		 unsigned int pMask,
    		 unsigned int numCnt,
             typename indexT,
    		 typename dataBufferType,
    		 typename outType,
    		 typename boxesVector_type,
    		 typename grid_smb_type,
    		 typename indexBuffer_type>
    __global__ void calc_exist_points_with_boxes(indexBuffer_type indexBuffer,
    											 boxesVector_type boxes,
    											 grid_smb_type grd,
    											 dataBufferType dataBuf,
    											 outType output,
    											 unsigned int stride_size)
    {
    	typedef typename dataBufferType::value_type AggregateT;
    	typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;
    	constexpr unsigned int blockSize = MaskBlockT::size;

        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x % blockSize;

        __shared__ int ato_cnt[numCnt];

        if (threadIdx.x < numCnt)
        {ato_cnt[threadIdx.x] = 0;}

        __syncthreads();

#ifdef SE_CLASS1

        if (numCnt > blockDim.x)
        {printf("Error calc_exist_points_with_boxes assertion failed numCnt >= blockDim.x  %d %d \n",numCnt,(int)blockDim.x);}

#endif

        if (dataBlockPos >= output.size())
        {return;}

        int predicate = dataBuf.template get<pMask>(dataBlockPos)[offset] & 0x1;
        // calculate point coord;
        indexT id = indexBuffer.template get<0>(dataBlockPos);
        grid_key_dx<dim,int> pnt = grd.InvLinId(id*grd.getBlockSize() + offset);
        Point<dim,int> p;

        for (int i = 0 ; i < dim ; i++)
        {p.get(i) = pnt.get(i);}

        // if this block intersect any box

        if (predicate == true)
        {
			for (int k = 0 ; k < boxes.size() ; k++ )
			{
				Box<dim,int> box = boxes.get(k);

				if (box.isInside(p) == true)
				{
					atomicAdd(&ato_cnt[k],1);
				}
			}
        }

        __syncthreads();

        if (threadIdx.x < boxes.size())
        {
        		output.template get<0>(dataBlockPos+threadIdx.x*stride_size) = ato_cnt[threadIdx.x];
        		output.template get<1>(dataBlockPos+threadIdx.x*stride_size) = (ato_cnt[threadIdx.x] != 0);
    	}
    }

    /*! \brief
     *
     * \param indexBuffer is the buffer of indexes
     * \param boxes is the set of boxes
     * \param grd contain information about the grid
     * \param dataBuff contain the information buffer
     * \param pack_output contain the index in the output buffer for the point to serialize
     * \param scan for each data-block it contain the number of points that fall inside the set of boxes
     *
     */
    template<unsigned int dim,
    		 unsigned int pMask,
    		 unsigned int numCnt,
             typename indexT,
    		 typename sparseGridType,
    		 typename packBufferType,
    		 typename scanType,
    		 typename scanItType,
    		 typename outputType,
    		 typename boxesVector_type,
    		 typename grid_smb_type>
    __global__ void get_exist_points_with_boxes(sparseGridType sparseGrid,
    											 boxesVector_type boxes,
    											 grid_smb_type grd,
    											 packBufferType pack_output,
    											 scanType scan,
    											 scanItType scan_it,
    											 outputType output)
    {
    	constexpr unsigned int blockSize = sparseGrid.getBlockSize();

        const unsigned int dataBlockPos = blockIdx.x;
        const unsigned int offset = threadIdx.x % blockSize;

        __shared__ int ato_cnt[numCnt];

        if (threadIdx.x < numCnt)
        {ato_cnt[threadIdx.x] = 0;}

        __syncthreads();

#ifdef SE_CLASS1

        if (numCnt > blockDim.x)
        {printf("Error get_exist_points_with_boxes assertion failed numCnt >= blockDim.x  %d %d \n",numCnt,(int)blockDim.x);}

#endif

        int predicate = sparseGrid.getMask(dataBlockPos,offset) & 0x1;
        // calculate point coord;
        indexT id = sparseGrid.getIndex(dataBlockPos);
        grid_key_dx<dim,int> pnt = grd.InvLinId(id*grd.getBlockSize() + offset);
        Point<dim,int> p_;

        for (int i = 0 ; i < dim ; i++)
        {p_.get(i) = pnt.get(i);}

        // if this block intersect any box

        if (predicate == true)
        {
			for (int k = 0 ; k < boxes.size() ; k++ )
			{
				Box<dim,int> box = boxes.get(k);

				if (box.isInside(p_) == true)
				{
					// We have an atomic counter for every packing box
					int p = atomicAdd(&ato_cnt[k] , 1);

					// we have a scan for every box
					const unsigned int dataBlockPosPack = scan.template get<1>(dataBlockPos + k*(sparseGrid.numBlocks() + 1));
					unsigned int sit = scan.template get<0>(dataBlockPos + k*(sparseGrid.numBlocks() + 1));
					int scan_id = scan.template get<0>(dataBlockPos + k*(sparseGrid.numBlocks() + 1)) + scan_it.template get<0>(k);
					output.template get<0>(scan_id + p) = (offset + dataBlockPos * blockSize) * numCnt + k;
					pack_output.template get<0>(scan_id + p) = p + sit;
				}
			}
        }
    }




    template<unsigned int dim,
    		 unsigned int blockSize,
    		 unsigned int blockEdgeSize,
    		 unsigned int nshifts,
    		 typename indexT,
    		 typename linearizer,
    		 typename shiftTypeVector,
    		 typename outputType>
    __global__ void convert_chunk_ids(indexT * ids,
    										int n_cnk,
    		                                linearizer gridGeoPack,
    		                                grid_key_dx<dim,int> origPack,
    		                                linearizer gridGeo,
    		                                grid_key_dx<dim,int> origUnpack,
    		                                outputType output,
    		                                shiftTypeVector shifts,
    		                                grid_key_dx<dim,int> sz,
    		                                int bs)
    {
    	// points
        const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

        if (p >= n_cnk)
        {return;}

        auto id = ids[p];

        for (int i = 0 ; i < nshifts ; i++)
        {
        	grid_key_dx<dim,int> pos = gridGeoPack.InvLinId(id,0) - origPack + origUnpack;

        	for (int j = 0 ; j < dim ; j++)
        	{
        		pos.set_d(j,pos.get(j) + shifts.template get<0>(i)[j]*blockEdgeSize);
        		if (pos.get(j) < 0)
        		{
        			pos.set_d(j,0);
        		}
        		if (pos.get(j) >= sz.get(j))
        		{
        			pos.set_d(j,pos.get(j) - blockEdgeSize);
        		}
        	}

        	auto plin = gridGeo.LinId(pos);

        	output.template get<0>(p*nshifts + i + bs) = plin / blockSize;
        }
    }

    template<typename vectorType>
    __global__ void set_one(vectorType vt)
    {
    	// points
        const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

        if (p >= vt.size())
        {return;}

        vt.template get<0>(p) = 1;
    }

    template<unsigned int pSegment,typename newMapType, typename mergeMapType,
             typename dataMapType, typename segmentOffsetType,
             typename outMapType>
    __global__ void construct_new_chunk_map(newMapType new_map, dataMapType dataMap,
    		                                mergeMapType merge_id, outMapType outMap,
    		                                segmentOffsetType segments_data, int start_p)
    {
    	// points
        const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

        if (p >= segments_data.size()-1)
        {return;}

        unsigned int st = segments_data.template get<pSegment>(p);

        int segmentSize = segments_data.template get<pSegment>(p + 1)
                          - segments_data.template get<pSegment>(p);

        for (int j = 0 ; j < segmentSize ; j++)
        {
        	int dm = dataMap.template get<0>(st+j);
        	new_map.template get<0>(dm) = outMap.template get<0>(p);
        }
    }

    template<unsigned int pMask, typename AggregateT, typename blockConvertType, typename newMapType, typename dataType_ptrs, typename dataType, unsigned int ... prp>
    __global__ void copy_packed_data_to_chunks(unsigned int * scan,
    										   unsigned short int * offsets,
    										   blockConvertType blc,
    										   newMapType new_map,
    										   dataType_ptrs data_ptrs,
    										   dataType data_buff,
    										   unsigned int n_cnk,
    										   unsigned int n_shf,
    										   unsigned int n_pnt,
    										   unsigned int i,
    										   unsigned int n_accu_cnk)
    {
    	// points
        const unsigned int p = blockIdx.x;

        if (p >= n_cnk)
        {return;}

        int scan_pp = scan[p];
        int n_block_pnt = scan[p+1] - scan_pp;

        if (threadIdx.x < n_block_pnt)
        {
        	unsigned short int off = offsets[scan[p] + threadIdx.x];

        	int conv = blc.template get<0>(i)[off];

        	unsigned short int off_c = conv & 0xFFFF;
        	unsigned short int shf_c = conv >> 16;

        	unsigned int pos_c = new_map.template get<0>(n_shf*p + shf_c + n_accu_cnk);

    		sparsegridgpu_unpack_impl<AggregateT, dataType ,prp ...>
    														spi(pos_c,off_c,data_buff,scan_pp + threadIdx.x,data_ptrs,n_pnt);

    		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(spi);

    		data_buff.template get<pMask>(pos_c)[off_c] |= 0x1;
        }
    }


    template<typename scanPointerType, typename scanType>
    __global__ void last_scan_point(scanPointerType scan_ptr, scanType scan,unsigned int stride, unsigned int n_pack)
    {
    	const unsigned int k = blockIdx.x * blockDim.x + threadIdx.x;

    	if (k >= n_pack)	{return;}

		unsigned int ppos = scan.template get<0>((k+1)*stride-1);
		unsigned int pos = scan.template get<1>((k+1)*stride-1);

		((unsigned int *)scan_ptr.ptr[k])[pos] = ppos;
    }

    template<unsigned int pMask,
             typename AggregateT,
    		 unsigned int n_it,
    		 unsigned int n_prp,
    		 typename indexT,
    		 typename pntBuff_type,
    		 typename pointOffset_type,
    		 typename SparseGridT,
    		 typename scan_type,
    		 unsigned int blockSize,
    		 unsigned int ... prp>
    __global__ void pack_data(pntBuff_type pntBuff,
    						  SparseGridT sparseGrid,
    						  scan_type scan,
    						  pointOffset_type point_offsets,
    						  arr_ptr<n_it> index_ptr,
    						  arr_ptr<n_it> scan_ptr,
    						  arr_arr_ptr<n_it,n_prp> * data_ptr,
    						  arr_ptr<n_it> offset_ptr,
    						  arr_ptr<n_it> mask_ptr,
    						  static_array<n_it,unsigned int> sar)
    {
        const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

        if (p >= pntBuff.size())
        {return;}

        const unsigned int pb = pntBuff.template get<0>(p);
        const unsigned int p_offset = point_offsets.template get<0>(p);

        const unsigned int k = pb % n_it;
        const unsigned int id = pb / n_it;

        const unsigned int dataBlockPos = id / blockSize;
        const unsigned int offset = id % blockSize;

		unsigned int ppos = scan.template get<0>(dataBlockPos + k*(sparseGrid.numBlocks() + 1));
		const unsigned int dataBlockPosPack = scan.template get<1>(dataBlockPos + k*(sparseGrid.numBlocks() + 1));

		sparsegridgpu_pack_impl<AggregateT, SparseGridT ,prp ...>
														spi(dataBlockPos,offset,sparseGrid,p_offset,data_ptr->ptr[k],sar.sa[k]);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(spi);

		((unsigned int *)scan_ptr.ptr[k])[dataBlockPosPack] = ppos;

		((indexT *)index_ptr.ptr[k])[dataBlockPosPack] = sparseGrid.getIndex(dataBlockPos);
		((short int *)offset_ptr.ptr[k])[p_offset] = offset;
		((unsigned char *)mask_ptr.ptr[k])[p_offset] = sparseGrid.getMask(dataBlockPos,offset);
    }
}

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_KERNELS_CUH
