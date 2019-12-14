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

    template<unsigned int dim,
    		 unsigned int pMask,
    		 unsigned int blockEdgeSize,
    		 typename dataBufferType,
    		 typename outType,
    		 typename boxesVector_type,
    		 typename grid_smb_type,
    		 typename indexBuffer_type>
    __global__ void calc_remove_points_chunks_boxes(indexBuffer_type indexBuffer,
    											 boxesVector_type boxes,
    											 grid_smb_type grd,
    											 dataBufferType dataBuf,
    											 outType output)
    {
    	typedef typename dataBufferType::value_type AggregateT;
    	typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;

        const unsigned int dataBlockPos = blockIdx.x * blockDim.x + threadIdx.x;

        if (dataBlockPos >= indexBuffer.size())
        {return;}

        auto id = indexBuffer.template get<0>(dataBlockPos);
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
    							typename dataBufferType,
    						    typename indexBufferType,
    						    typename grid_smb_type,
    						    typename activeCntType,
    						    typename boxesType>
    __global__ void remove_points(indexBufferType indexBuffer,
    							  grid_smb_type grd,
    							  dataBufferType dataBuffer,
    							  activeCntType active_blocks,
    							  boxesType boxes)
    {
    	typedef typename dataBufferType::value_type AggregateT;
    	typedef BlockTypeOf<AggregateT, pMask> MaskBlockT;
    	constexpr unsigned int blockSize = MaskBlockT::size;

        const unsigned int dataBlockPos = active_blocks.template get<0>(blockIdx.x);
        const unsigned int offset = threadIdx.x % blockSize;

        if (dataBlockPos >= dataBuffer.size()-1)
        {return;}

        int predicate = dataBuffer.template get<pMask>(dataBlockPos)[offset] & 0x1;
        // calculate point coord;
        auto id = indexBuffer.template get<0>(dataBlockPos);
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
					dataBuffer.template get<pMask>(dataBlockPos)[offset] = 0;
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

        if (numCnt >= blockDim.x)
        {printf("Error calc_exist_points_with_boxes assertion failed numCnt >= blockDim.x  %d %d \n",numCnt,blockDim.x);}

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
    		 typename dataBufferType,
    		 typename packBufferType,
    		 typename scanType,
    		 typename scanItType,
    		 typename outputType,
    		 typename boxesVector_type,
    		 typename grid_smb_type,
    		 typename indexBuffer_type>
    __global__ void get_exist_points_with_boxes(indexBuffer_type indexBuffer,
    											 boxesVector_type boxes,
    											 grid_smb_type grd,
    											 dataBufferType dataBuf,
    											 packBufferType pack_output,
    											 scanType scan,
    											 scanItType scan_it,
    											 outputType output)
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

        if (numCnt >= blockDim.x)
        {printf("Error get_exist_points_with_boxes assertion failed numCnt >= blockDim.x  %d %d \n",numCnt,blockDim.x);}

#endif

        int predicate = dataBuf.template get<pMask>(dataBlockPos)[offset] & 0x1;
        // calculate point coord;
        indexT id = indexBuffer.template get<0>(dataBlockPos);
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
					const unsigned int dataBlockPosPack = scan.template get<1>(dataBlockPos + k*(indexBuffer.size() + 1));
					unsigned int sit = scan.template get<0>(dataBlockPos + k*(indexBuffer.size() + 1));
					int scan_id = scan.template get<0>(dataBlockPos + k*(indexBuffer.size() + 1)) + scan_it.template get<0>(k);
					output.template get<0>(scan_id + p) = (offset + dataBlockPos * blockSize) * numCnt + k;
					pack_output.template get<0>(scan_id + p) = p + sit;
				}
			}
        }
    }

    template<unsigned int pMask, typename add_data_type>
    __global__ void resetMask(add_data_type add_data, unsigned int start)
    {
    	// points
        const unsigned int bid = blockIdx.x + start;
        const unsigned int tid = threadIdx.x;

        if (bid >= add_data.size())
        {return;}

    	add_data.template get<pMask>(bid)[tid] = 0;
    }

    template<unsigned int blockSize,
    		 unsigned int pMask,
    		 typename AggregateT,
    		 typename indexT,
    		 typename point_buffer,
    		 typename chunk_arr_type,
    		 typename convertion_type,
    		 typename add_index_type,
    		 typename data_ptrs_type,
    		 typename add_data_type,
    		 unsigned int ... prp>
    __global__ void fill_add_buffer(indexT * ptr_id,
    								point_buffer pts,
    								chunk_arr_type chunk_arr,
    								convertion_type conv,
    								add_index_type add_index,
    								data_ptrs_type data_ptrs,
    								add_data_type add_data,
    								unsigned int start)
    {
    	// points
        const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

        if (p >= pts.size())
        {return;}

        auto dataBlockPos = conv.template get<0>(p);
        auto ord = chunk_arr.template get<1>(p);
        auto plin = chunk_arr.template get<0>(p);

        auto dataBlockId = plin / blockSize;
        short int offset = plin % blockSize;

        add_index.template get<0>(dataBlockPos + start) = dataBlockId;

		sparsegridgpu_unpack_impl<AggregateT, add_data_type ,prp ...>
														spi(dataBlockPos + start,offset,add_data,ord,data_ptrs,pts.size());

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(spi);

		add_data.template get<pMask>(dataBlockPos + start)[offset] = 1;
    }

    template<unsigned int blockSize, typename vector_type, typename output_type>
    __global__ void mark_unpack_chunks(vector_type vd, output_type output)
    {
    	// points
        const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

        if (p >= vd.size())
        {return;}

        auto cp = vd.template get<0>(p) / blockSize;
        decltype(cp) cp_p1;

        if (p == vd.size()-1)	{cp_p1 = -1;}
        else {cp_p1 = vd.template get<0>(p+1) / blockSize;}

        output.template get<0>(p) = cp_p1 != cp;
    }

    template<unsigned int dim,
    		 unsigned int blockSize,
    		 typename indexT,
    		 typename linearizer,
    		 typename segType,
    		 typename outputType>
    __global__ void convert_chunk_alignment(indexT * ids, short int * offsets, unsigned int * scan,
    										unsigned int n_seg,
    										segType segments,
    		                                linearizer gridGeoPack,
    		                                grid_key_dx<dim,int> origPack,
    		                                linearizer gridGeo,
    		                                grid_key_dx<dim,int> origUnpack,
    		                                outputType output)
    {
    	// points
        const unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;

        if (p >= output.size())
        {return;}

        // get the chunk index
        int cid = segments.template get<0>(p);

        auto id = ids[cid];

		short int offset = offsets[p];
		grid_key_dx<dim,int> pos = gridGeoPack.InvLinId(id,offset) - origPack + origUnpack;

		auto plin = gridGeo.LinId(pos);

		output.template get<0>(p) = plin;
		output.template get<1>(p) = p;
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

    template<typename scanPointerType, typename scanType>
    __global__ void last_scan_point(scanPointerType scan_ptr, scanType scan,unsigned int stride, unsigned int n_pack)
    {
    	const unsigned int k = blockIdx.x * blockDim.x + threadIdx.x;

    	if (k >= n_pack)	{return;}

		unsigned int ppos = scan.template get<0>((k+1)*stride-1);
		unsigned int pos = scan.template get<1>((k+1)*stride-1);

		((unsigned int *)scan_ptr.ptr[k])[pos] = ppos;
    }

    template<typename AggregateT,
    		 unsigned int n_it,
    		 unsigned int n_prp,
    		 typename indexT,
    		 typename pntBuff_type,
    		 typename pointOffset_type,
    		 typename indexBuffer_type,
    		 typename dataBuffer_type,
    		 typename scan_type,
    		 unsigned int blockSize,
    		 unsigned int ... prp>
    __global__ void pack_data(pntBuff_type pntBuff,
    						  dataBuffer_type dataBuff,
    						  indexBuffer_type indexBuff,
    						  scan_type scan,
    						  pointOffset_type point_offsets,
    						  arr_ptr<n_it> index_ptr,
    						  arr_ptr<n_it> scan_ptr,
    						  arr_arr_ptr<n_it,n_prp> * data_ptr,
    						  arr_ptr<n_it> offset_ptr,
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

		unsigned int ppos = scan.template get<0>(dataBlockPos + k*(indexBuff.size() + 1));
		const unsigned int dataBlockPosPack = scan.template get<1>(dataBlockPos + k*(indexBuff.size() + 1));

		sparsegridgpu_pack_impl<AggregateT, dataBuffer_type ,prp ...>
														spi(dataBlockPos,offset,dataBuff,p_offset,data_ptr->ptr[k],sar.sa[k]);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(spi);

		((unsigned int *)scan_ptr.ptr[k])[dataBlockPosPack] = ppos;

		((indexT *)index_ptr.ptr[k])[dataBlockPosPack] = indexBuff.template get<0>(dataBlockPos);
		((short int *)offset_ptr.ptr[k])[p_offset] = offset;
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
