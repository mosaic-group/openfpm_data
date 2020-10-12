//
// Created by tommaso on 6/06/19.
//

#ifndef OPENFPM_PDATA_SPARSEGRIDGPU_HPP
#define OPENFPM_PDATA_SPARSEGRIDGPU_HPP

constexpr int BLOCK_SIZE_STENCIL = 128;

#include <cstdlib>
#include <SparseGridGpu/BlockMapGpu.hpp>
#include <Grid/iterators/grid_skin_iterator.hpp>
#include <SparseGridGpu/Geometry/grid_smb.hpp>
#include "SparseGridGpu_ker.cuh"
#include "SparseGridGpu_kernels.cuh"
#include "Iterators/SparseGridGpu_iterator_sub.hpp"
#include "Geometry/grid_zmb.hpp"
#include "util/stat/common_statistics.hpp"
#include "Iterators/SparseGridGpu_iterator.hpp"
#include "util/cuda/moderngpu/kernel_load_balance.hxx"
#include "Space/SpaceBox.hpp"

#if defined(OPENFPM_DATA_ENABLE_IO_MODULE) || defined(PERFORMANCE_TEST)
#include "VTKWriter/VTKWriter.hpp"
#endif

#ifdef OPENFPM_PDATA
#include "VCluster/VCluster.hpp"
#endif

constexpr int NO_ITERATOR_INIT = 0;

// todo: Move all the following utils into some proper file inside TemplateUtils

enum tag_boundaries
{
	NO_CALCULATE_EXISTING_POINTS,
	CALCULATE_EXISTING_POINTS
};

template<unsigned int dim>
struct default_edge
{
	typedef boost::mpl::int_<2> type;
};


template<>
struct default_edge<1>
{
	typedef boost::mpl::int_<64> type;
};

template<>
struct default_edge<2>
{
	typedef boost::mpl::int_<8> type;
};

template<>
struct default_edge<3>
{
	typedef boost::mpl::int_<4> type;
};

template<typename T>
struct type_identity
{
	typedef T type;
};

template<typename T, unsigned int dim, unsigned int blockEdgeSize>
struct process_data_block
{
	typedef type_identity<DataBlock<T,IntPow<blockEdgeSize,dim>::value>> type;
};

template<typename T, unsigned int dim, unsigned int blockEdgeSize, unsigned int N1>
struct process_data_block<T[N1],dim,blockEdgeSize>
{
	typedef type_identity<DataBlock<T,IntPow<blockEdgeSize,dim>::value>[N1]> type;
};

template<unsigned int dim, unsigned int blockEdgeSize, typename ... aggr_list>
struct aggregate_transform_datablock_impl
{
	typedef aggregate<typename process_data_block<aggr_list,dim,blockEdgeSize>::type::type ...> type;
};

template<unsigned int dim, unsigned int blockEdgeSize, typename aggr>
struct aggregate_convert
{
};

template<unsigned int dim, unsigned int blockEdgeSize, typename ... types>
struct aggregate_convert<dim,blockEdgeSize,aggregate<types ...>>
{
	typedef typename aggregate_transform_datablock_impl<dim,blockEdgeSize,types ...>::type type;
};

template<typename aggr>
struct aggregate_add
{
};

template<typename ... types>
struct aggregate_add<aggregate<types ...>>
{
	typedef aggregate<types ..., unsigned char> type;
};

/////////////

enum StencilMode
{
    STENCIL_MODE_INPLACE = 1,
    STENCIL_MODE_INPLACE_NO_SHARED = 3
};

/*! \brief get the type of the block
 *
 *
 */
template<typename SGridGpu, unsigned int prp, unsigned int stencil_size>
struct GetCpBlockType
{
	typedef cp_block<typename boost::mpl::at<typename SGridGpu::device_grid_type::value_type::type,boost::mpl::int_<prp>>::type,
	         stencil_size,
	         typename vmpl_sum_constant<2*stencil_size,typename vmpl_create_constant<SGridGpu::dims,SGridGpu::device_grid_type::blockEdgeSize_>::type >::type,
	         SGridGpu::device_grid_type::dims> type;
};

#include "encap_num.hpp"

/*! \brief get the type of the insertBlock
 *
 *
 */
template<typename SGridGpu>
struct GetAddBlockType
{
	typedef enc_num<typename SGridGpu::device_grid_type::insert_encap> type;
};

/*! \brief Check if is padding
 *
 *
 */
template<unsigned int dim>
struct NNfull_is_padding_impl
{
	__device__ inline static bool is_padding()
	{
		printf("NNfull_is_padding_impl with dim: %d not implemented yet \n",dim);

		return false;
	}
};

/*! \brief Check if is padding
 *
 *
 */
template<>
struct NNfull_is_padding_impl<3>
{
	template<typename sparseGrid_type, typename coord_type, typename Mask_type,unsigned int eb_size>
	__device__ inline static bool is_padding(sparseGrid_type & sparseGrid, coord_type & coord, Mask_type (& enlargedBlock)[eb_size])
	{
		bool isPadding_ = false;
		for (int i = 0 ; i < 3 ; i++)
		{
			for (int j = 0 ; j < 3 ; j++)
			{
				for (int k = 0 ; k < 3 ; k++)
				{
					grid_key_dx<3,int> key;

					key.set_d(0,i-1);
					key.set_d(1,j-1);
					key.set_d(2,k-1);

					auto nPlusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, key);
					typename std::remove_all_extents<Mask_type>::type neighbourPlus = enlargedBlock[nPlusId];
					isPadding_ = isPadding_ || (!sparseGrid.exist(neighbourPlus));
					if (isPadding_) break;
				}
			}
		}
		return isPadding_;
	}
};


/*! \brief Check if is padding
 *
 *
 */
template<>
struct NNfull_is_padding_impl<2>
{
	template<typename sparseGrid_type, typename coord_type, typename Mask_type,unsigned int eb_size>
	__device__ inline static bool is_padding(sparseGrid_type & sparseGrid, coord_type & coord, Mask_type (& enlargedBlock)[eb_size])
	{
		bool isPadding_ = false;
		for (int i = 0 ; i < 3 ; i++)
		{
			for (int j = 0 ; j < 3 ; j++)
			{
				grid_key_dx<2,int> key;

				key.set_d(0,i-1);
				key.set_d(1,j-1);

				auto nPlusId = sparseGrid.getNeighbourLinIdInEnlargedBlock(coord, key);
				typename std::remove_all_extents<Mask_type>::type neighbourPlus = enlargedBlock[nPlusId];
				isPadding_ = isPadding_ || (!sparseGrid.exist(neighbourPlus));
				if (isPadding_) break;
			}
		}
		return isPadding_;
	}
};

template<unsigned int dim>
struct NNFull
{
	static const int nNN = IntPow<3, dim>::value;

	template<typename indexT, typename blockCoord_type, typename blockMap_type, typename SparseGrid_type>
	__device__ static inline indexT getNNpos(blockCoord_type & blockCoord,
								  blockMap_type & blockMap,
								  SparseGrid_type & sparseGrid,
								  const unsigned int offset)
	{
        // point to the background element
        int neighbourPos = blockMap.size();
        if (offset < nNN && offset != nNN / 2)
        {
        	int cnt = offset;
        	for (int i = 0 ; i < dim ; i++)
        	{
        		int dPos = cnt % 3;
        		cnt /= 3;
        		blockCoord.set_d(i, blockCoord.get(i) + dPos - 1);
        	}

            neighbourPos = blockMap.get_sparse(sparseGrid.getBlockLinId(blockCoord)).id;
        }
        return neighbourPos;
	}

	template<typename indexT, unsigned int blockEdgeSize, typename coordType>
	__host__ static inline indexT getNNskin(coordType & coord, int stencilSupportRadius)
	{
		// linearize the coord

		indexT neighbourNum = 0;

		indexT accu = 1;
		for(int i = 0 ; i < dim ; i++)
		{
			int c = static_cast<int>(coord.get(i)) - static_cast<int>(stencilSupportRadius);
			if (c < 0)
			{
				neighbourNum += 0;
			}
            else if (c >= blockEdgeSize)
            {
                neighbourNum += 2*accu;
            }
            else
            {
            	neighbourNum += accu;
            }
			accu *= 3;
		}

        return neighbourNum;
	}


	template<typename sparseGrid_type, typename coord_type, typename Mask_type,unsigned int eb_size>
	__device__ static inline bool isPadding(sparseGrid_type & sparseGrid, coord_type & coord, Mask_type (& enlargedBlock)[eb_size])
	{
		return NNfull_is_padding_impl<3>::template is_padding(sparseGrid,coord,enlargedBlock);
	}

	/*! \brief given a coordinate writtel in local coordinate for a given it return the neighborhood chunk position and the offset
	 *        in the neighborhood chunk
	 *
	 * \param coord local coordinated
	 * \param NN_index index of the neighborhood chunk
	 * \param offset_nn offset in local coordinates
	 *
	 * \return true if it is inside false otherwise
	 *
	 */
	template<unsigned int blockEdgeSize, typename indexT2>
	__device__ static inline bool getNNindex_offset(grid_key_dx<dim,indexT2> & coord, unsigned int & NN_index, unsigned int & offset_nn)
	{
		bool out = false;
		NN_index = 0;
		offset_nn = 0;

		int cnt = 1;
		int cnt_off = 1;
		for (unsigned int i = 0 ; i < dim ; i++)
		{
			int p = 1 - ((int)(coord.get(i) < 0)) + ((int)(coord.get(i) >= (int)blockEdgeSize));

			NN_index += p*cnt;

			offset_nn += (coord.get(i) + (1 - p)*(int)blockEdgeSize)*cnt_off;

			cnt *= 3;
			cnt_off *= blockEdgeSize;

			out |= (p != 1);
		}

		return out;
	}
};



template<unsigned int nNN_, unsigned int nLoop_>
struct ct_par
{
	static const unsigned int nNN = nNN_;
	static const unsigned int nLoop = nLoop_;
};

template<typename copy_type>
struct copy_prop_to_vector_block_impl
{
	template<typename T, typename dst_type, typename src_type>
	static inline void copy(src_type & src, dst_type & dst, unsigned int bPos)
	{
		dst.template get<T::value>() = src.template get<T::value>()[bPos];
	}
};

template<typename copy_type,unsigned int N1>
struct copy_prop_to_vector_block_impl<copy_type[N1]>
{
	template<typename T, typename dst_type, typename src_type>
	static inline void copy(src_type & src, dst_type & dst, unsigned int bPos)
	{
		for (int i = 0 ; i < N1 ; i++)
		{
			dst.template get<T::value>()[i] = src.template get<T::value>()[i][bPos];
		}
	}
};

template<typename Tsrc,typename Tdst>
class copy_prop_to_vector_block
{
	//! source
	Tsrc src;

	//! destination
	Tdst dst;

	size_t pos;

	unsigned int bPos;

public:

	copy_prop_to_vector_block(Tsrc src, Tdst dst,size_t pos, size_t bPos)
	:src(src),dst(dst),pos(pos),bPos(bPos)
	{}

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename std::remove_reference<decltype(dst.template get<T::value>())>::type copy_rtype;

		copy_prop_to_vector_block_impl<copy_rtype>::template copy<T>(src,dst,bPos);

		//meta_copy<copy_rtype>::meta_copy_(src.template get<T::value>()[bPos],dst.template get<T::value>());
	}

};



template<typename AggregateT, unsigned int n_it, unsigned int ... prp>
class data_ptr_fill
{
	typedef typename to_boost_vmpl<prp...>::type vprp;

	//! data pointers
	void * base_ptr;

	mutable size_t tot = 0;

	size_t i;

	size_t sz = 0;

	arr_arr_ptr<n_it,sizeof...(prp)> & arrs;

public:

	data_ptr_fill(void * base_ptr,size_t i,  arr_arr_ptr<n_it,sizeof...(prp)> & arrs, size_t sz)
	:base_ptr(base_ptr),i(i),sz(sz),arrs(arrs)
	{}

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<vprp,T>::type prp_cp;

		// Remove the reference from the type to copy
		typedef typename boost::mpl::at<typename AggregateT::type,prp_cp>::type pack_type;

		arrs.ptr[i][T::value] = (void *)((((unsigned char *)base_ptr) + tot));

		tot += sz * sizeof(pack_type);
	}

};

template<typename SparseGridType>
struct sparse_grid_section
{
	SparseGridType * grd;

	// Source box
	Box<SparseGridType::dims,size_t> src;

	// destination Box
	Box<SparseGridType::dims,size_t> dst;

	sparse_grid_section(SparseGridType & grd,const Box<SparseGridType::dims,size_t> & src, const Box<SparseGridType::dims,size_t> & dst)
	:grd(&grd),src(src),dst(dst)
	{}
};


template<unsigned int dim,
		 typename AggregateT,
		 unsigned int blockEdgeSize = default_edge<dim>::type::value,
		 unsigned int threadBlockSize = 128,
		 typename indexT=long int,
		 template<typename> class layout_base=memory_traits_inte,
		 typename linearizer = grid_smb<dim, blockEdgeSize>>
class SparseGridGpu : public BlockMapGpu<
        typename aggregate_convert<dim,blockEdgeSize,AggregateT>::type,
        threadBlockSize, indexT, layout_base>
{
public:

	 static constexpr unsigned int dims = dim;

private:

	typedef BlockMapGpu<
	        typename aggregate_convert<dim,blockEdgeSize,AggregateT>::type,
	        threadBlockSize, indexT, layout_base> BMG;

    const static unsigned char PADDING_BIT = 1;
	typedef typename aggregate_convert<dim,blockEdgeSize,AggregateT>::type AggregateBlockT;
    linearizer gridGeometry;
    grid_sm<dim, int> extendedBlockGeometry;
    grid_sm<dim, int> gridSize;
    unsigned int stencilSupportRadius;
    unsigned int ghostLayerSize;
    int req_index;
    int req_index_swp;
    int req_index_swp_r;

    AggregateT bck;

    typedef SparseGridGpu<dim,AggregateT,blockEdgeSize,threadBlockSize,indexT,layout_base,linearizer> self;

    // Queue of remove sections
    openfpm::vector_gpu<Box<dim,unsigned int>> rem_sects;

    // Queue of copy sections
    openfpm::vector<sparse_grid_section<self>> copySect;

    CudaMemory mem;

    // pointers for unpack
    openfpm::vector<void *> index_ptrs;
    openfpm::vector<void *> index_ptrs_swp;
    openfpm::vector<void *> index_ptrs_swp_r;
    openfpm::vector<void *> scan_ptrs;
    openfpm::vector<void *> scan_ptrs_swp;
    openfpm::vector<void *> scan_ptrs_swp_r;
    openfpm::vector<void *> data_ptrs;
    openfpm::vector<void *> data_ptrs_swp;
    openfpm::vector<void *> data_ptrs_swp_r;
    openfpm::vector<void *> offset_ptrs;
    openfpm::vector<void *> offset_ptrs_swp;
    openfpm::vector<void *> offset_ptrs_swp_r;
    openfpm::vector<void *> mask_ptrs;
    openfpm::vector<void *> mask_ptrs_swp;
    openfpm::vector<void *> mask_ptrs_swp_r;

    // pointers for copyRemove
    openfpm::vector<void *> offset_ptrs_cp;
    openfpm::vector<void *> offset_ptrs_cp_swp;
    openfpm::vector<void *> offset_ptrs_cp_swp_r;
    openfpm::vector<void *> scan_ptrs_cp;
    openfpm::vector<void *> scan_ptrs_cp_swp;
    openfpm::vector<void *> scan_ptrs_cp_swp_r;
    openfpm::vector<void *> data_base_ptr_cp;
    openfpm::vector<void *> data_base_ptr_cp_swp;
    openfpm::vector<void *> data_base_ptr_cp_swp_r;
    openfpm::vector<int> n_cnk_cp;
    openfpm::vector<int> n_cnk_cp_swp;
    openfpm::vector<int> n_cnk_cp_swp_r;
    openfpm::vector<int> n_pnt_cp;
    openfpm::vector<int> n_pnt_cp_swp;
    openfpm::vector<int> n_pnt_cp_swp_r;
    openfpm::vector<int> n_shifts_cp;
    openfpm::vector<int> n_shift_cp_swp;
    openfpm::vector<int> n_shift_cp_swp_r;
	typedef typename aggregate_convert<dim,blockEdgeSize,aggregate<int>>::type convertAggr;

	// Map to convert blocks from missaligned chunks
	openfpm::vector_gpu<convertAggr> convert_blk;
	openfpm::vector_gpu<convertAggr> convert_blk_swp;
	openfpm::vector_gpu<convertAggr> convert_blk_swp_r;
	openfpm::vector<Box<dim,size_t>> box_cp;
	openfpm::vector<Box<dim,size_t>> box_cp_swp;
	openfpm::vector<Box<dim,size_t>> box_cp_swp_r;

    //! set of existing points
    //! the formats is id/blockSize = data block poosition id % blockSize = offset
    openfpm::vector_gpu<aggregate<indexT>> e_points;
    openfpm::vector_gpu<aggregate<indexT>> e_points_swp;
    openfpm::vector_gpu<aggregate<indexT>> e_points_swp_r;

    //! Helper array to pack points
    openfpm::vector_gpu<aggregate<unsigned int>> pack_output;
    openfpm::vector_gpu<aggregate<unsigned int>> pack_output_swp;
    openfpm::vector_gpu<aggregate<unsigned int>> pack_output_swp_r;

    //! For stencil in a block-wise computation we have to load blocks + ghosts area. The ghost area live in neighborhood blocks
    //! For example the left ghost margin live in the right part of the left located neighborhood block, the right margin live in the
    //! left part of the of the right located neighborhood block, the top ...
    //! The first index indicate the index of the point in the block + ghost area, the second index indicate the correspondent neighborhood
    //! index (in a star like 0 mean negative x 1 positive x, 1 mean negative y and so on)
    openfpm::vector_gpu<aggregate<short int,short int>> ghostLayerToThreadsMapping;

    openfpm::vector_gpu<aggregate<indexT>> nn_blocks;

    //! temporal
    mutable openfpm::vector_gpu<aggregate<indexT,unsigned int>> tmp;
    mutable openfpm::vector_gpu<aggregate<indexT,unsigned int>> tmp_swp;
    mutable openfpm::vector_gpu<aggregate<indexT,unsigned int>> tmp_swp_r;

    //! temporal 2
    mutable openfpm::vector_gpu<aggregate<indexT>> tmp2;

    //! temporal 3
    mutable openfpm::vector_gpu<aggregate<indexT>> tmp3;

    //! contain the scan of the point for each iterator
    mutable openfpm::vector_gpu<aggregate<indexT>> scan_it;

    //! Map between the (Last) added chunks and their position in chunks data
    mutable openfpm::vector_gpu<aggregate<int>> new_map;
    mutable openfpm::vector_gpu<aggregate<int>> new_map_swp;
    mutable openfpm::vector_gpu<aggregate<int>> new_map_swp_r;

    //! the set of all sub-set to pack
    mutable openfpm::vector_gpu<Box<dim,int>> pack_subs;
    mutable openfpm::vector_gpu<Box<dim,int>> pack_subs_swp;
    mutable openfpm::vector_gpu<Box<dim,int>> pack_subs_swp_r;

    //! Size of the index vector packed. These varaible are used to understand if the option
    //! KEEP_GEOMETRY can be used keep geometry option infact require that when we record the
    //! packing variables the number of chunks (and chunks indexes) does not change
    mutable int index_size_swp = -1;
    mutable int index_size_swp_r = -1;

    //! links of the padding points with real points of a coarse sparsegrid
    openfpm::vector_gpu<aggregate<size_t>> links_up;

    //! scan offsets of the links down
    openfpm::vector_gpu<aggregate<unsigned int>> link_dw_scan;

    //! links of the padding points with real points of a finer sparsegrid
    openfpm::vector_gpu<aggregate<int,short int>> link_dw;

    //! scan offsets of the links down
    openfpm::vector_gpu<aggregate<unsigned int>> link_up_scan;

    //! links of the padding points with real points of a finer sparsegrid
    openfpm::vector_gpu<aggregate<int,short int>> link_up;

    //! Memory to remove copy finalize
    ExtPreAlloc<CudaMemory> * prAlloc_prp;

    bool findNN = false;

    inline void swap_internal_remote()
    {
		n_cnk_cp_swp_r.swap(n_cnk_cp);
		n_pnt_cp_swp_r.swap(n_pnt_cp);
		n_shift_cp_swp_r.swap(n_shifts_cp);
		convert_blk_swp_r.swap(convert_blk);
		box_cp_swp_r.swap(box_cp);
		new_map_swp_r.swap(new_map);
    }

    inline void swap_internal_local()
    {
		offset_ptrs_cp_swp.swap(offset_ptrs_cp);
		scan_ptrs_cp_swp.swap(scan_ptrs_cp);
		data_base_ptr_cp_swp.swap(data_base_ptr_cp);
		n_cnk_cp_swp.swap(n_cnk_cp);
		n_pnt_cp_swp.swap(n_pnt_cp);
		n_shift_cp_swp.swap(n_shifts_cp);
		convert_blk_swp.swap(convert_blk);
		box_cp_swp.swap(box_cp);
		new_map_swp.swap(new_map);
    }

    inline void swap_local_pack()
    {
		index_ptrs_swp.swap(index_ptrs);
		scan_ptrs_swp.swap(scan_ptrs);
		data_ptrs_swp.swap(data_ptrs);
		offset_ptrs_swp.swap(offset_ptrs);
		mask_ptrs_swp.swap(mask_ptrs);

		e_points_swp.swap(e_points);
		pack_output_swp.swap(pack_output);
		tmp_swp.swap(tmp);

		pack_subs_swp.swap(pack_subs);
		index_size_swp = private_get_index_array().size();
    }

    inline void swap_remote_pack()
    {
		index_ptrs_swp_r.swap(index_ptrs);
		scan_ptrs_swp_r.swap(scan_ptrs);
		data_ptrs_swp_r.swap(data_ptrs);
		offset_ptrs_swp_r.swap(offset_ptrs);
		mask_ptrs_swp_r.swap(mask_ptrs);

		e_points_swp_r.swap(e_points);
		pack_output_swp_r.swap(pack_output);
		tmp_swp_r.swap(tmp);

		pack_subs_swp_r.swap(pack_subs);
		//req_index_swp_r = req_index;
		index_size_swp_r = private_get_index_array().size();
    }

protected:
    static constexpr unsigned int blockSize = BlockTypeOf<AggregateBlockT, 0>::size;
    typedef AggregateBlockT AggregateInternalT;

public:

	//! it define that this data-structure is a grid
	typedef int yes_i_am_grid;

    static constexpr unsigned int blockEdgeSize_ = blockEdgeSize;

    typedef linearizer grid_info;

    typedef linearizer linearizer_type;

    template<typename Tfunc> using layout_mfunc = memory_traits_inte<Tfunc>;

    typedef sparse_grid_gpu_index<self> base_key;

    typedef decltype(std::declval<BMG>().toKernel().insertBlock(0)) insert_encap;

    /*! \brief return the size of the grid
     *
     * \return Return the size of the grid
     *
     */
    inline size_t size() const
    {
        return this->countExistingElements();
    }

    /*! \brief This is a meta-function return which type of sub iterator a grid produce
     *
     * \return the type of the sub-grid iterator
     *
     */
    template <typename stencil = no_stencil>
    static SparseGridGpu_iterator_sub<dim,self> type_of_subiterator()
    {
        return SparseGridGpu_iterator_sub<dim,self>();
    }

    /*! \brief This is a meta-function return which type of iterator a grid produce
     *
     * \return the type of the sub-grid iterator
     *
     */
    static SparseGridGpu_iterator<dim,self> type_of_iterator()
    {
        return SparseGridGpu_iterator<dim,self>(std::declval<self>());
    }

    template<typename dim3T>
    inline static int dim3SizeToInt(dim3T d)
    {
        return d.x * d.y * d.z;
    }

    inline static int dim3SizeToInt(size_t d)
    {
        return d;
    }

    inline static int dim3SizeToInt(unsigned int d)
    {
        return d;
    }

    template<typename ... v_reduce>
    void flush(mgpu::ofp_context_t &context, flush_type opt = FLUSH_ON_HOST)
    {
        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
                ::template flush<v_reduce ...>(context, opt);

        findNN = false;
    }



    void saveUnpackVariableIfNotKeepGeometry(int opt, bool is_unpack_remote)
    {
    	if (is_unpack_remote == true)
    	{swap_internal_remote();}

    	if (is_unpack_remote == false)
    	{swap_internal_local();}
    }

    void RestoreUnpackVariableIfKeepGeometry(int opt, bool is_unpack_remote)
    {
		if (opt & KEEP_GEOMETRY && is_unpack_remote == true)
		{swap_internal_remote();}

		if (opt & KEEP_GEOMETRY && is_unpack_remote == false)
		{swap_internal_local();}
    }


    void savePackVariableIfNotKeepGeometry(int opt, bool is_pack_remote)
    {
		if (is_pack_remote == false)
		{
			swap_local_pack();
			req_index_swp = req_index;
		}

		if (is_pack_remote == true)
		{
			swap_remote_pack();
			req_index_swp_r = req_index;
		}
    }

    void RestorePackVariableIfKeepGeometry(int opt, bool is_pack_remote)
    {
		if (opt & KEEP_GEOMETRY && is_pack_remote == false)
		{
			swap_local_pack();
			req_index = req_index_swp;
		}

		if (opt & KEEP_GEOMETRY && is_pack_remote == true)
		{
			swap_remote_pack();
			req_index = req_index_swp_r;
		}
    }

    template<unsigned int n_it>
    void calculatePackingPointsFromBoxes(int opt,size_t tot_pnt)
    {
		if (!(opt & KEEP_GEOMETRY))
		{
	    	auto & indexBuffer = private_get_index_array();
	    	auto & dataBuffer = private_get_data_array();

			e_points.resize(tot_pnt);
			pack_output.resize(tot_pnt);

			ite_gpu<1> ite;

			ite.wthr.x = indexBuffer.size();
			ite.wthr.y = 1;
			ite.wthr.z = 1;
			ite.thr.x = getBlockSize();
			ite.thr.y = 1;
			ite.thr.z = 1;

			// Launch a kernel that count the number of element on each chunks
			CUDA_LAUNCH((SparseGridGpuKernels::get_exist_points_with_boxes<dim,
																		BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
																		n_it,
																		indexT>),
					 ite,
					 indexBuffer.toKernel(),
					 pack_subs.toKernel(),
					 gridGeometry,
					 dataBuffer.toKernel(),
					 pack_output.toKernel(),
					 tmp.toKernel(),
					 scan_it.toKernel(),
					 e_points.toKernel());
		}
    }

private:

    void computeSizeOfGhostLayer()
    {
        unsigned int term1 = 1;
        for (int i = 0; i < dim; ++i)
        {
            term1 *= blockEdgeSize + 2 * stencilSupportRadius;
        }
        unsigned int term2 = 1;
        for (int i = 0; i < dim; ++i)
        {
            term2 *= blockEdgeSize;
        }
        ghostLayerSize = term1 - term2;
    }

    void allocateGhostLayerMapping()
    {
    	ghostLayerToThreadsMapping.resize(ghostLayerSize);
    }

    template<typename stencil_type>
    void computeGhostLayerMapping()
    {
        size_t dimensions[dim],
                origin[dim],
                innerDomainBegin[dim], innerDomainEnd[dim],
                outerBoxBegin[dim], outerBoxEnd[dim],
                bc[dim];
        for (int i = 0; i < dim; ++i)
        {
            dimensions[i] = blockEdgeSize + 2 * stencilSupportRadius;
            origin[i] = 0;
            innerDomainBegin[i] = stencilSupportRadius - 1;
            innerDomainEnd[i] = dimensions[i] - stencilSupportRadius;
            outerBoxBegin[i] = origin[i];
            outerBoxEnd[i] = dimensions[i];
            bc[i] = NON_PERIODIC;
        }
        grid_sm<dim, void> enlargedGrid;
        enlargedGrid.setDimensions(dimensions);
        Box<dim, size_t> outerBox(outerBoxBegin, outerBoxEnd);
        Box<dim, size_t> innerBox(innerDomainBegin, innerDomainEnd);

        grid_skin_iterator_bc<dim> gsi(enlargedGrid, innerBox, outerBox, bc);

        unsigned int i = 0;
        while (gsi.isNext())
        {
            auto coord = gsi.get();
            assert(i < ghostLayerSize);
            mem_id linId = enlargedGrid.LinId(coord);
            // Mapping
            ghostLayerToThreadsMapping.template get<gt>(i) = linId;
            // Now compute the neighbour position to use
            ghostLayerToThreadsMapping.template get<nt>(i) = stencil_type::template getNNskin<indexT,blockEdgeSize>(coord,stencilSupportRadius);
            //
            ++i;
            ++gsi;
        }
        assert(i == ghostLayerSize);

        ghostLayerToThreadsMapping.template hostToDevice<gt,nt>();
    }

    void initialize(const size_t (& res)[dim])
    {
    	gridGeometry = linearizer(res);

        computeSizeOfGhostLayer();
        allocateGhostLayerMapping();
        computeGhostLayerMapping<NNStar<dim>>();

        size_t extBlockDims[dim];
        for (int d=0; d<dim; ++d)
        {
            extBlockDims[d] = blockEdgeSize + 2*stencilSupportRadius;
        }
        extendedBlockGeometry.setDimensions(extBlockDims);
        gridSize.setDimensions(res);
    }


    template <typename stencil, typename... Args>
    void applyStencilInPlace(const Box<dim,int> & box, StencilMode & mode,Args... args)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;

        // NOTE: Here we want to work only on one data chunk per block!
        constexpr unsigned int chunksPerBlock = 1;
        const unsigned int localThreadBlockSize = dataChunkSize * chunksPerBlock;
        const unsigned int threadGridSize = numScalars % localThreadBlockSize == 0
                                            ? numScalars / localThreadBlockSize
                                            : 1 + numScalars / localThreadBlockSize;

        constexpr unsigned int nLoop = UIntDivCeil<(IntPow<blockEdgeSize + 2, dim>::value - IntPow<blockEdgeSize, dim>::value), (blockSize * chunksPerBlock)>::value; // todo: This works only for stencilSupportSize==1

        CUDA_LAUNCH_DIM3((SparseGridGpuKernels::applyStencilInPlace
                <dim,
                BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
                stencil>),
                threadGridSize, localThreadBlockSize,
                        box,
                        indexBuffer.toKernel(),
                        dataBuffer.toKernel(),
                        this->template toKernelNN<stencil::stencil_type::nNN, nLoop>(),
                        args...);
    }

    template <typename stencil, typename... Args>
    void applyStencilInPlaceNoShared(const Box<dim,int> & box, StencilMode & mode,Args... args)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;

        auto ite = e_points.getGPUIterator(BLOCK_SIZE_STENCIL);

        CUDA_LAUNCH((SparseGridGpuKernels::applyStencilInPlaceNoShared
                <dim,
                BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
                stencil>),
                ite,
                        box,
                        indexBuffer.toKernel(),
                        dataBuffer.toKernel(),
                        this->template toKernelNN<stencil::stencil_type::nNN, 0>(),
                        args...);
    }

    template<typename ids_type>
    void fill_chunks_boxes(openfpm::vector<SpaceBox<dim,double>> & chunks_box, ids_type & chunk_ids, Point<dim,double> & spacing, Point<dim,double> & offset)
    {
    	for (int i = 0 ; i < chunk_ids.size() ; i++)
    	{
    		SpaceBox<dim,double> box;

    		auto c_pos = gridGeometry.InvLinId(chunk_ids.template get<0>(i)*blockSize);

    		for (int j = 0 ; j < dim ; j++)
    		{
    			box.setLow(j,c_pos.get(j) * spacing[j] - 0.5*spacing[j] + offset.get(j)*spacing[j]);
    			box.setHigh(j,(c_pos.get(j) + blockEdgeSize)*spacing[j] - 0.5*spacing[j] + offset.get(j)*spacing[j]);
    		}

    		chunks_box.add(box);
    	}
    }

    template<typename MemType, unsigned int ... prp>
    void preUnpack(ExtPreAlloc<MemType> * prAlloc_prp, mgpu::ofp_context_t & ctx, int opt)
    {
		if ((opt & rem_copy_opt::KEEP_GEOMETRY) == false)
		{
			// Convert the packed chunk ids

			prAlloc_prp->reset();
			Unpack_stat ups;

			for (size_t i = 0 ; i < copySect.size() ; i++)
			{
				auto sub_it = this->getIterator(copySect.get(i).dst.getKP1(),copySect.get(i).dst.getKP2(),NO_ITERATOR_INIT);

				copySect.get(i).grd->template addAndConvertPackedChunkToTmp<prp ...>(*prAlloc_prp,sub_it,ups,ctx);
			}
		}
    }


	template<unsigned int ... prp>
	void removeCopyToFinalize_phase1(mgpu::ofp_context_t & ctx, int opt)
	{
		if ((opt & rem_copy_opt::KEEP_GEOMETRY) == false)
		{removePoints(ctx);}
	}

	template<unsigned int ... prp>
	void removeCopyToFinalize_phase2(mgpu::ofp_context_t & ctx, int opt)
	{
		// Pack information
		Pack_stat sts;

		if ((opt & rem_copy_opt::KEEP_GEOMETRY) == false)
		{
			this->packReset();

			size_t req = 0;
			// First we do counting of point to copy (as source)

			for (size_t i = 0 ; i < copySect.size() ; i++)
			{
				auto sub_it = this->getIterator(copySect.get(i).src.getKP1(),copySect.get(i).src.getKP2(),NO_ITERATOR_INIT);

				this->packRequest(sub_it,req);
			}

			this->template packCalculate<prp...>(req,ctx);

			mem.resize(req);

			// Create an object of preallocated memory for properties
			prAlloc_prp = new ExtPreAlloc<CudaMemory>(req,mem);
			prAlloc_prp->incRef();

			for (size_t i = 0 ; i < copySect.size() ; i++)
			{
				auto sub_it = this->getIterator(copySect.get(i).src.getKP1(),copySect.get(i).src.getKP2(),NO_ITERATOR_INIT);

				this->pack<prp ...>(*prAlloc_prp,sub_it,sts);
			}
		}
		else
		{
			size_t req = mem.size();

			// Create an object of preallocated memory for properties
			prAlloc_prp = new ExtPreAlloc<CudaMemory>(req,mem);
			prAlloc_prp->incRef();
		}

		this->template packFinalize<prp ...>(*prAlloc_prp,sts,opt,false);

		preUnpack<CudaMemory,prp ...>(prAlloc_prp,ctx,opt);

		prAlloc_prp->decRef();
		delete prAlloc_prp;
	}

	template<unsigned int ... prp>
	void removeCopyToFinalize_phase3(mgpu::ofp_context_t & ctx, int opt, bool is_unpack_remote)
	{
		ite_gpu<1> ite;

		if ((opt & rem_copy_opt::KEEP_GEOMETRY) == false)
		{
			if (tmp2.size() == 0)
			{return;}

			// Fill the add buffer given tmp and than flush

			setGPUInsertBuffer(tmp2.size(),1ul);

			auto & add_buff = this->blockMap.private_get_vct_add_index();
			add_buff.swap(tmp2);

			auto & nadd_buff = this->blockMap.private_get_vct_nadd_index();
			ite = nadd_buff.getGPUIterator();
			CUDA_LAUNCH(SparseGridGpuKernels::set_one,ite,nadd_buff.toKernel());

			int sz_b =  this->private_get_index_array().size();

			this->template flush<sLeft_<prp>...>(ctx,flush_type::FLUSH_ON_DEVICE);

			// get the map of the new inserted elements

			auto & m_map = this->getMergeIndexMapVector();
			auto & a_map = this->getMappingVector();
			auto & o_map = this->getSegmentToOutMap();
			auto & segments_data = this->getSegmentToMergeIndexMap();

			new_map.resize(a_map.size());

			// construct new to old map

			ite = segments_data.getGPUIterator();

			if (ite.nblocks() != 0)
			CUDA_LAUNCH(SparseGridGpuKernels::construct_new_chunk_map<1>,ite,new_map.toKernel(),a_map.toKernel(),m_map.toKernel(),o_map.toKernel(),segments_data.toKernel(),sz_b);

			convert_blk.template hostToDevice<0>();
		}
		else
		{
			ite.wthr.x = 1;
			ite.wthr.y = 1;
			ite.wthr.z = 1;

			ite.thr.x = 1;
			ite.thr.y = 1;
			ite.thr.z = 1;
		}

		// Restore
		RestoreUnpackVariableIfKeepGeometry(opt,is_unpack_remote);

		// for each packed chunk

		size_t n_accu_cnk = 0;
		for (size_t i = 0 ; i < n_cnk_cp.size() ; i++)
		{
			arr_arr_ptr<1,sizeof...(prp)> data;
			size_t n_pnt = n_pnt_cp.get(i);

			void * data_base_ptr = data_base_ptr_cp.get(i);
			data_ptr_fill<AggregateT,1,prp...> dpf(data_base_ptr,0,data,n_pnt);
			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prp)>>(dpf);

			ite.wthr.x = n_cnk_cp.get(i);

			// calculate best number of threads
			Box<dim,size_t> ub = box_cp.get(i);

			ite.thr.x = 1;
			for (int j = 0 ; j < dim ; j++)
			{
				size_t l = ub.getHigh(j) - ub.getLow(j) + 1;

				if (l >= blockEdgeSize)
				{ite.thr.x *= blockEdgeSize;}
				else
				{ite.thr.x *= l;}
			}

			// copy to new (1 block for each packed chunk)
			if (ite.nblocks() != 0 && ite.thr.x != 0)
			{
				auto & chunks = private_get_data_array();

				CUDA_LAUNCH((SparseGridGpuKernels::copy_packed_data_to_chunks<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
						 	 	 	 	 	 	 	 	 	 	 	 	 	  AggregateT,decltype(convert_blk.toKernel()),decltype(new_map.toKernel()),
						                                                      decltype(data),decltype(chunks.toKernel()),prp... >),ite,
																				 (unsigned int *)scan_ptrs_cp.get(i),
																				 (unsigned short int *)offset_ptrs_cp.get(i),
																				 convert_blk.toKernel(),
																				 new_map.toKernel(),
																				 data,
																				 chunks.toKernel(),
																				 n_cnk_cp.get(i),
																				 n_shifts_cp.get(i),
																				 n_pnt_cp.get(i),
																				 i,
																				 n_accu_cnk);
			}

			n_accu_cnk += n_cnk_cp.get(i)*n_shifts_cp.get(i);
		}

		// Save
		saveUnpackVariableIfNotKeepGeometry(opt,is_unpack_remote);
	}

    template<unsigned int n_it, unsigned int ... prp>
    void pack_sg_implement(ExtPreAlloc<CudaMemory> & mem,
						   Pack_stat & sts,
						   int opt,
						   bool is_pack_remote)
    {
    	arr_ptr<n_it> index_ptr;
    	arr_arr_ptr<n_it,sizeof...(prp)> data_ptr;
    	arr_ptr<n_it> scan_ptr;
    	arr_ptr<n_it> offset_ptr;
    	arr_ptr<n_it> mask_ptr;
    	static_array<n_it,unsigned int> sar;

    	auto & indexBuffer = private_get_index_array();
    	auto & dataBuffer = private_get_data_array();

		if (req_index != pack_subs.size())
		{std::cerr << __FILE__ << ":" << __LINE__ << " error the packing request number differ from the number of packed objects " << req_index << "  " << pack_subs.size() << std::endl;}

    	size_t tot_pnt = 0;
    	size_t tot_cnk = 0;

    	sparsegridgpu_pack_request<AggregateT,prp ...> spq;
    	boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prp)>>(spq);

    	// Calculate total points

    	for (size_t i = 0 ; i < pack_subs.size() ; i++)
    	{
    		size_t n_pnt = tmp.template get<0>((i+1)*(indexBuffer.size() + 1)-1);
    		sar.sa[i] = n_pnt;
    		tot_pnt += n_pnt;
    	}

    	// CUDA require aligned access, here we suppose 8 byte alligned and we ensure 8 byte aligned after
    	// the cycle
    	for (size_t i = 0 ; i < pack_subs.size() ; i++)
    	{
    		size_t n_cnk = tmp.template get<1>((i+1)*(indexBuffer.size() + 1)-1);

    		// fill index_ptr data_ptr scan_ptr
    		index_ptr.ptr[i] = index_ptrs.get(i);
    		scan_ptr.ptr[i] = scan_ptrs.get(i);

    		// for all properties fill the data pointer

    		data_ptr_fill<AggregateT,n_it,prp...> dpf(data_ptrs.get(i),i,data_ptr,tmp.template get<0>((i+1)*(indexBuffer.size() + 1)-1));
    		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prp)>>(dpf);

    		offset_ptr.ptr[i] = offset_ptrs.get(i);
    		mask_ptr.ptr[i] = mask_ptrs.get(i);

    		tot_cnk += n_cnk;
    	}

    	ite_gpu<1> ite;

    	if (tot_pnt != 0)
    	{
    		calculatePackingPointsFromBoxes<n_it>(opt,tot_pnt);

			ite = e_points.getGPUIterator();

			// Here we copy the array of pointer of properties into a CudaMemory array

			CudaMemory mem;
			mem.allocate(sizeof(data_ptr));

			// copy
			arr_arr_ptr<n_it,sizeof...(prp)> * arr_data = (arr_arr_ptr<n_it,sizeof...(prp)> *)mem.getPointer();

			for(int i = 0 ; i < n_it ; i++)
			{
				for (int j = 0 ; j < sizeof...(prp) ; j++)
				{
					arr_data->ptr[i][j] = data_ptr.ptr[i][j];
				}
			}

			mem.hostToDevice();

			CUDA_LAUNCH((SparseGridGpuKernels::pack_data<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
					               AggregateT,
								   n_it,
								   sizeof...(prp),
								   indexT,
								   decltype(e_points.toKernel()),
								   decltype(pack_output.toKernel()),
								   decltype(indexBuffer.toKernel()),
								   decltype(dataBuffer.toKernel()),
								   decltype(tmp.toKernel()),
								   self::blockSize,
								   prp ...>),
								   ite,
								   e_points.toKernel(),
								   dataBuffer.toKernel(),
								   indexBuffer.toKernel(),
								   tmp.toKernel(),
								   pack_output.toKernel(),
								   index_ptr,
								   scan_ptr,
								   (arr_arr_ptr<n_it,sizeof...(prp)> *)mem.getDevicePointer(),
								   offset_ptr,
								   mask_ptr,
								   sar);
    	}

    	ite.wthr.x = 1;
    	ite.wthr.y = 1;
    	ite.wthr.z = 1;
    	ite.thr.x = pack_subs.size();
    	ite.thr.y = 1;
    	ite.thr.z = 1;

    	if (pack_subs.size() != 0)
		{CUDA_LAUNCH(SparseGridGpuKernels::last_scan_point,ite,scan_ptr,tmp.toKernel(),indexBuffer.size()+1,pack_subs.size());}
    }


	/*! \brief unpack the sub-grid object
	 *
	 * \tparam prp properties to unpack
	 *
	 * \param mem preallocated memory from where to unpack the object
	 * \param sub sub-grid iterator
	 * \param obj object where to unpack
	 *
	 */
	template<unsigned int ... prp, typename S2>
	void addAndConvertPackedChunkToTmp(ExtPreAlloc<S2> & mem,
				SparseGridGpu_iterator_sub<dim,self> & sub_it,
				Unpack_stat & ps,
				mgpu::ofp_context_t &context)
	{
    	sparsegridgpu_pack_request<AggregateT,prp ...> spq;
    	boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prp)>>(spq);

		// First get the number of chunks

		size_t n_cnk;

		// Unpack the number of chunks
		mem.deviceToHost(ps.getOffset(),ps.getOffset() + sizeof(size_t) + 2*dim*sizeof(int));
		Unpacker<size_t,S2>::unpack(mem,n_cnk,ps);

		grid_key_dx<dim,int> origPack_pnt;
		grid_key_dx<dim,int> origPack_cnk;
		size_t sz[dim];

		// Unpack origin of the chunk indexing
		for (int i = 0 ; i < dim ; i++)
		{
			int tmp;
			Unpacker<int,S2>::unpack(mem,tmp,ps);
			origPack_cnk.set_d(i,((int)(tmp / blockEdgeSize))*blockEdgeSize);
			origPack_pnt.set_d(i,tmp);
		}

		for (int i = 0 ; i < dim ; i++)
		{
			int tmp;
			Unpacker<int,S2>::unpack(mem,tmp,ps);
			sz[i] = tmp;
		}

		size_t actual_offset = n_cnk*sizeof(indexT);
		// get the id pointers
		indexT * ids = (indexT *)((unsigned char *)mem.getDevicePointer() + ps.getOffset());
		unsigned int * scan = (unsigned int *)((unsigned char *)mem.getDevicePointer() + ps.getOffset() + n_cnk*sizeof(indexT));

		mem.deviceToHost(ps.getOffset() + actual_offset + n_cnk*sizeof(unsigned int),
						 ps.getOffset() + actual_offset + n_cnk*sizeof(unsigned int) + sizeof(unsigned int));



		// Unpack number of points
		// calculate the number of total points
		size_t n_pnt = *(unsigned int *)((unsigned char *)mem.getPointer() + ps.getOffset() + actual_offset + n_cnk*sizeof(unsigned int));
		actual_offset += align_number(sizeof(indexT),(n_cnk+1)*sizeof(unsigned int));

		void * data_base_ptr = (void *)((unsigned char *)mem.getDevicePointer() + ps.getOffset() + actual_offset );

		actual_offset += align_number(sizeof(indexT),n_pnt*(spq.point_size));

		short int * offsets = (short int *)((unsigned char *)mem.getDevicePointer() + ps.getOffset() + actual_offset);

		offset_ptrs_cp.add(offsets);
		scan_ptrs_cp.add(scan);
		n_cnk_cp.add(n_cnk);
		n_pnt_cp.add(n_pnt);
		data_base_ptr_cp.add(data_base_ptr);

		Box<dim,size_t> bx;

		for (int i = 0 ; i < dim ; i++)
		{
			bx.setLow(i,sub_it.getStart().get(i));
			bx.setHigh(i,sub_it.getStop().get(i));
		}

		box_cp.add(bx);

		actual_offset += align_number(sizeof(indexT),n_pnt*sizeof(short));

		if (n_cnk != 0)
		{
			openfpm::vector_gpu<aggregate<int[dim]>> shifts;

			int n_shift = 1;
			shifts.add();

			for (int i = 0 ; i < dim ; i++)
			{shifts.last().template get<0>()[i] = 0;}

			for (int i = 0 ; i < dim ; i++)
			{
				int op_q = origPack_pnt.get(i) % blockEdgeSize;
				int ou_q = sub_it.getStart().get(i) % blockEdgeSize;
				int quot = abs(ou_q - op_q) % blockEdgeSize;
				int squot = openfpm::math::sgn(ou_q - op_q);
				if (quot != 0)
				{
					n_shift *= 2;

					int sz = shifts.size();
					for (int j = 0 ; j < sz ; j++)
					{
						shifts.add();
						for (int k = 0 ; k < dim ; k++)
						{
							shifts.last().template get<0>()[k] = shifts.template get<0>(j)[k] + ((i == k)?squot:0);
						}
					}
				}
			}

			shifts.template hostToDevice<0>();

			linearizer gridGeoPack(sz);

			int bs = 0;
			size_t sz[1] = {n_cnk};
			grid_sm<1,void> g(sz);
			auto ite = g.getGPUIterator();

			grid_key_dx<dim,int> sz_g;
			grid_key_dx<dim,int> origUnpack_cnk;

			for (int i = 0 ; i < dim ; i++)
			{
				sz_g.set_d(i,gridGeometry.getSize()[i]);
				origUnpack_cnk.set_d(i,(int)(sub_it.getStart().get(i) / blockEdgeSize)*blockEdgeSize);
			}

			bs = tmp2.size();
			tmp2.resize(tmp2.size() + n_cnk * shifts.size());

			n_shifts_cp.add(shifts.size());

			switch (shifts.size())
			{
			case 1:
				// Calculate for each chunk the indexes where they should go + active points
				CUDA_LAUNCH((SparseGridGpuKernels::convert_chunk_ids<dim,blockSize,blockEdgeSize,1,indexT>),ite,ids,
															  n_cnk,
															  gridGeoPack,origPack_cnk,
															  gridGeometry,origUnpack_cnk,
															  tmp2.toKernel(),
															  shifts.toKernel(),
															  sz_g,
															  bs);
				break;
			case 2:
				// Calculate for each chunk the indexes where they should go + active points
				CUDA_LAUNCH((SparseGridGpuKernels::convert_chunk_ids<dim,blockSize,blockEdgeSize,2,indexT>),ite,ids,
															  n_cnk,
															  gridGeoPack,origPack_cnk,
															  gridGeometry,origUnpack_cnk,
															  tmp2.toKernel(),
															  shifts.toKernel(),
															  sz_g,
															  bs);
				break;
			case 4:
				// Calculate for each chunk the indexes where they should go + active points
				CUDA_LAUNCH((SparseGridGpuKernels::convert_chunk_ids<dim,blockSize,blockEdgeSize,4,indexT>),ite,ids,
															  n_cnk,
															  gridGeoPack,origPack_cnk,
															  gridGeometry,origUnpack_cnk,
															  tmp2.toKernel(),
															  shifts.toKernel(),
															  sz_g,
															  bs);
				break;
			case 8:
				// Calculate for each chunk the indexes where they should go + active points
				CUDA_LAUNCH((SparseGridGpuKernels::convert_chunk_ids<dim,blockSize,blockEdgeSize,8,indexT>),ite,ids,
															  n_cnk,
															  gridGeoPack,origPack_cnk,
															  gridGeometry,origUnpack_cnk,
															  tmp2.toKernel(),
															  shifts.toKernel(),
															  sz_g,
															  bs);
				break;
			}

			convertChunkIds(offsets,origPack_pnt,sub_it);
		}
		else
		{
			convert_blk.add();
			n_shifts_cp.add(0);
		}

		actual_offset += align_number(sizeof(indexT),n_pnt*sizeof(unsigned char));

		ps.addOffset(actual_offset);
	}

	/*! \brief convert the offset index from the packed to the add buffer
	 *
	 *
	 */
	template<typename origPackType, typename IteratorType>
	void convertChunkIds(short int * offset, origPackType & origPack, IteratorType & sub_it)
	{
		int quot_diff[dim];
		for (int i = 0 ; i < dim ; i++)
		{
			int op_q = origPack.get(i) % blockEdgeSize;
			int ou_q = sub_it.getStart().get(i) % blockEdgeSize;
			int quot = abs(ou_q - op_q) % blockEdgeSize;
			quot_diff[i] = openfpm::math::sgn(ou_q - op_q)*quot;
		}

		convert_blk.add();

		// Create conversion block

		for (int j = 0 ; j < this->blockSize ; j++)
		{
			int offset = 0;
			int bpos = 0;
			int bp_c = 1;
			int pos = 0;
			int pos_c = 1;

			int x = j;
			for (int i = 0 ; i < dim ; i++)
			{
				int c = x % blockEdgeSize;

				if (quot_diff[i] + c < 0)
				{
					offset += pos_c*(quot_diff[i] + c + blockEdgeSize);
					bpos += bp_c*1;
				}
				else if (quot_diff[i] + c >= blockEdgeSize)
				{
					offset += pos_c*(quot_diff[i] + c - blockEdgeSize);
					bpos += bp_c*1;
				}
				else
				{
					offset += pos_c*(quot_diff[i] + c);
				}

				pos += pos_c*c;
				pos_c *= blockEdgeSize;
				bp_c *= (quot_diff[i] != 0)?2:1;
				x /= blockEdgeSize;
			}

			convert_blk.template get<0>(convert_blk.size()-1)[pos] = (bpos << 16) + offset;
		}
	}

public:

    typedef AggregateT value_type;

    typedef self device_grid_type;

    SparseGridGpu()
	:stencilSupportRadius(1)
    {};

    /*! \brief resize the SparseGrid
     *
     * \param res indicate the resolution in each dimension
     *
     */
    void resize(size_t (& res)[dim])
    {
    	initialize(res);
    }

    /*! \brief Constructor from glock geometry
     *
     *
     */
    SparseGridGpu(const size_t (& res)[dim], unsigned int stencilSupportRadius = 1)
    :stencilSupportRadius(stencilSupportRadius)
    {
    	initialize(res);
    };

    /*! \brief Constructor from glock geometry
     *
     *
     */
    SparseGridGpu(linearizer & gridGeometry, unsigned int stencilSupportRadius = 1)
            : gridGeometry(gridGeometry),
              stencilSupportRadius(stencilSupportRadius)
    {
    	initialize(gridGeometry.getSize());
    };

    SparseGridGpu_ker
            <
                    dim,
                    blockEdgeSize,
                    typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                    ct_par<0,1>,
                    indexT,
                    layout_base,
                    decltype(extendedBlockGeometry),
                    linearizer,
                    AggregateT
            > toKernel()
    {
        SparseGridGpu_ker
                <
                        dim,
                        blockEdgeSize,
                        typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                        ct_par<0,1>,
                        indexT,
                        layout_base,
                        decltype(extendedBlockGeometry),
                        linearizer,
                        AggregateT
                > toKer(
                BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.toKernel(),
                gridGeometry,
                extendedBlockGeometry,
                stencilSupportRadius,
                ghostLayerToThreadsMapping.toKernel(),
                nn_blocks.toKernel(),
                e_points.toKernel(),
                ghostLayerSize,
                bck);
        return toKer;
    }

    template<unsigned int nNN, unsigned int nLoop>
    SparseGridGpu_ker
            <
                    dim,
                    blockEdgeSize,
                    typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                    ct_par<nNN,nLoop>,
                    indexT,
                    layout_base,
                    decltype(extendedBlockGeometry),
                    linearizer,
                    AggregateT
            > toKernelNN()
    {
        SparseGridGpu_ker
                <
                        dim,
                        blockEdgeSize,
                        typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT,
                        ct_par<nNN,nLoop>,
                        indexT,
                        layout_base,
                        decltype(extendedBlockGeometry),
                        linearizer,
                        AggregateT
                > toKer(
                BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.toKernel(),
                gridGeometry,
                extendedBlockGeometry,
                stencilSupportRadius,
                ghostLayerToThreadsMapping.toKernel(),
                nn_blocks.toKernel(),
                e_points.toKernel(),
                ghostLayerSize,
                bck);
        return toKer;
    }

    /*! \brief Return the grid information object
     *
     * \return grid information object
     *
     */
    linearizer & getGrid()
    {
    	return gridGeometry;
    }

    /*! \brief Set the neighborhood type
     *
     * \tparam stencil_type Type of stencil
     *
     */
    template<typename stencil_type>
    void setNNType()
    {
    	computeGhostLayerMapping<stencil_type>();
    }


    constexpr static unsigned int getBlockEdgeSize()
    {
        return blockEdgeSize;
    }

    constexpr unsigned int getBlockSize() const
    {
        return blockSize;
    }

    // Geometry
    template<typename CoordT>
    inline size_t getLinId(CoordT &coord)
    {
        return gridGeometry.LinId(coord);
    }

    inline grid_key_dx<dim, int> getCoord(size_t linId) const
    {
        return gridGeometry.InvLinId(linId);
    }

    inline ite_gpu<dim> getGridGPUIterator(const grid_key_dx<dim, int> & start, const grid_key_dx<dim, int> & stop, size_t n_thr = threadBlockSize)
    {
    	return gridSize.getGPUIterator(start,stop,n_thr);
    }

    /*! \brief Get an element using the point coordinates
     *
     * \tparam p property index
     *
     * \param coord point coordinates
     *
     * \return the element
     *
     */
    template<typename CoordT>
    base_key get_sparse(const grid_key_dx<dim,CoordT> & coord) const
    {
    	base_key k(*this,0,0);

    	const auto & blockMap = this->private_get_blockMap();

    	auto glid = gridGeometry.LinId(coord);

    	auto bid = glid / blockSize;
    	auto lid = glid % blockSize;

    	auto key = blockMap.get_sparse(bid);

    	k.set_cnk_pos_id(key.id);
    	k.set_data_id(lid);

        return k;
    }

    /*! \brief Get an element using the point coordinates
     *
     * \tparam p property index
     *
     * \param coord point coordinates
     *
     * \return the element
     *
     */
    template<unsigned int p, typename CoordT>
    auto get(const grid_key_dx<dim,CoordT> & coord) const -> const ScalarTypeOf<AggregateBlockT, p> &
    {
        return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::template get<p>(gridGeometry.LinId(coord));
    }

    /*! \brief Get an element using sparse_grid_gpu_index (using this index it guarantee that the point exist)
     *
     * \tparam p property index
     *
     * \param element
     *
     * \return the element
     *
     */
    template<unsigned int p>
    auto get(const sparse_grid_gpu_index<self> & coord) const -> const ScalarTypeOf<AggregateBlockT, p> &
    {
        return private_get_data_array().template get<p>(coord.get_cnk_pos_id())[coord.get_data_id()];
    }

    /*! \brief This function check if keep geometry is possible for this grid
     *
     * \return true if skip labelling is possible
     *
     */
    bool isSkipLabellingPossible()
    {
    	return (index_size_swp_r == private_get_index_array().size()) && (index_size_swp == private_get_index_array().size());
    }

    /*! \brief Get an element using sparse_grid_gpu_index (using this index it guarantee that the point exist)
     *
     * \tparam p property index
     *
     * \param element
     *
     * \return the element
     *
     */
    template<unsigned int p>
    auto get(const sparse_grid_gpu_index<self> & coord) -> ScalarTypeOf<AggregateBlockT, p> &
    {
        return private_get_data_array().template get<p>(coord.get_cnk_pos_id())[coord.get_data_id()];
    }

    /*! \brief Return the flag of the point
     *
     * It indicate for example is if the point is a padding point (internaly it return the pMask flag)
     *
     * \return the flag
     *
     */
    unsigned char getFlag(const sparse_grid_gpu_index<self> & coord) const
    {
    	return private_get_data_array().template get<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>(coord.get_cnk_pos_id())[coord.get_data_id()];
    }

    template<unsigned int p, typename CoordT>
    auto insert(const CoordT &coord) -> ScalarTypeOf<AggregateBlockT, p> &
    {
        return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::template insert<p>(gridGeometry.LinId(coord));
    }

	/*! \brief construct link between levels
	 *
	 * \praram grid_up grid level up
	 * \param grid_dw grid level down
	 *
	 */
    void construct_link(self & grid_up, self & grid_dw, mgpu::ofp_context_t &context)
    {
/*        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

    	ite_gpu<1> ite;

    	ite.wthr.x = indexBuffer.size();
    	ite.wthr.y = 1;
    	ite.wthr.z = 1;

    	ite.thr.x = getBlockSize();
    	ite.thr.y = 1;
    	ite.thr.z = 1;

    	openfpm::vector_gpu<aggregate<unsigned int>> output;
    	output.resize(indexBuffer.size() + 1);

    	CUDA_LAUNCH((SparseGridGpuKernels::link_construct<dim,
    								BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
    								blockSize>),ite,grid_up.toKernel(),this->toKernel(),output.toKernel());


    	openfpm::scan((unsigned int *)output.template getDeviceBuffer<0>(),output.size(),(unsigned int *)output.template getDeviceBuffer<0>(),context);

    	output.template deviceToHost<0>(output.size()-1,output.size()-1);

    	unsigned int np_lup = output.template get<0>(output.size()-1);

    	links_up.resize(np_lup);

    	CUDA_LAUNCH((SparseGridGpuKernels::link_construct_insert<dim,
    								BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
    								blockSize>),ite,grid_up.toKernel(),this->toKernel(),output.toKernel(),links_up.toKernel());

*/
    }

    /*! \brief Get the offsets for each point of the links down
     *
     * \return the offsets of the links down
     *
     */
    openfpm::vector_gpu<aggregate<unsigned int>> & getDownLinksOffsets()
    {
    	return link_dw_scan;
    }

    /*! \brief Get the links down for each point
     *
     * \return the links dow for each point
     *
     */
    openfpm::vector_gpu<aggregate<int,short int>> & getDownLinks()
    {
    	return link_dw;
    }

    /*! \brief Get the offsets for each point of the links up
     *
     * \return the offsets of the links up
     *
     */
    openfpm::vector_gpu<aggregate<unsigned int>> & getUpLinksOffsets()
    {
    	return link_up_scan;
    }

    /*! \brief Get the links up for each point
     *
     * \return the links up for each point
     *
     */
    openfpm::vector_gpu<aggregate<int,short int>> & getUpLinks()
    {
    	return link_up;
    }

	/*! \brief construct link on the down level
	 *
	 * \param grid_dw grid level down
	 * \param db domain box
	 * \param p_dw point offset when you go down
	 * \param gpu context
	 *
	 */
    void construct_link_dw(self & grid_dw, const Box<dim,int> & db_, Point<dim,int> p_dw, mgpu::ofp_context_t &context)
    {
    	Box<dim,int> db = db_;

        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        // Count padding points

        // First we count the padding points
    	ite_gpu<1> ite;

    	ite.wthr.x = indexBuffer.size();
    	ite.wthr.y = 1;
    	ite.wthr.z = 1;

    	ite.thr.x = getBlockSize();
    	ite.thr.y = 1;
    	ite.thr.z = 1;

    	openfpm::vector_gpu<aggregate<unsigned int>> output;
    	output.resize(indexBuffer.size()+1);

    	output.fill<0>(0);

    	CUDA_LAUNCH((SparseGridGpuKernels::count_paddings<dim,
    								BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
    								blockSize>),ite,this->toKernel(),output.toKernel(),db);



    	openfpm::scan((unsigned int *)output.template getDeviceBuffer<0>(),output.size(),(unsigned int *)output.template getDeviceBuffer<0>(),context);

    	output.template deviceToHost<0>(output.size()-1,output.size()-1);
        unsigned int padding_points = output.template get<0>(output.size()-1);

        // get the padding points

        openfpm::vector_gpu<aggregate<unsigned int,short int>> pd_points;
        pd_points.resize(padding_points);

        CUDA_LAUNCH((SparseGridGpuKernels::collect_paddings<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>),ite,this->toKernel(),output.toKernel(),pd_points.toKernel(),db);

        // Count number of link down for padding points

        // Calculate ghost

    	link_dw_scan.resize(padding_points+1);
    	link_dw_scan.fill<0>(0);

    	ite = link_dw_scan.getGPUIterator();

    	CUDA_LAUNCH((SparseGridGpuKernels::link_construct_dw_count<dim,
    								BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
    								blockSize>),
    								ite,pd_points.toKernel(),grid_dw.toKernel(),this->toKernel(),link_dw_scan.toKernel(),p_dw);

    	openfpm::scan((unsigned int *)link_dw_scan.template getDeviceBuffer<0>(),link_dw_scan.size(),(unsigned int *)link_dw_scan.template getDeviceBuffer<0>(),context);

    	link_dw_scan.template deviceToHost<0>(link_dw_scan.size()-1,link_dw_scan.size()-1);

    	size_t np_ldw = link_dw_scan.template get<0>(link_dw_scan.size()-1);

    	link_dw.resize(np_ldw);

    	CUDA_LAUNCH((SparseGridGpuKernels::link_construct_insert_dw<dim,
    								BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
    								blockSize>),ite,pd_points.toKernel(),grid_dw.toKernel(),this->toKernel(),link_dw_scan.toKernel(),link_dw.toKernel(),p_dw);

    	link_dw_scan.resize(link_dw_scan.size()-1);
    }

	/*! \brief construct link on the up levels
	 *
	 * \praram grid_up grid level up
	 *
	 */
    void construct_link_up(self & grid_up,  const Box<dim,int> & db_, Point<dim,int> p_up, mgpu::ofp_context_t &context)
    {
    	Box<dim,int> db = db_;

        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        // Count padding points

        // First we count the padding points
    	ite_gpu<1> ite;

    	ite.wthr.x = indexBuffer.size();
    	ite.wthr.y = 1;
    	ite.wthr.z = 1;

    	ite.thr.x = getBlockSize();
    	ite.thr.y = 1;
    	ite.thr.z = 1;

    	openfpm::vector_gpu<aggregate<unsigned int>> output;
    	output.resize(indexBuffer.size()+1);

    	output.fill<0>(0);

    	CUDA_LAUNCH((SparseGridGpuKernels::count_paddings<dim,
    								BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
    								blockSize>),ite,this->toKernel(),output.toKernel(),db);



    	openfpm::scan((unsigned int *)output.template getDeviceBuffer<0>(),output.size(),(unsigned int *)output.template getDeviceBuffer<0>(),context);

    	output.template deviceToHost<0>(output.size()-1,output.size()-1);
        unsigned int padding_points = output.template get<0>(output.size()-1);

        // get the padding points

        openfpm::vector_gpu<aggregate<unsigned int,short int>> pd_points;
        pd_points.resize(padding_points);

        CUDA_LAUNCH((SparseGridGpuKernels::collect_paddings<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>),ite,this->toKernel(),output.toKernel(),pd_points.toKernel(),db);

        // Count number of link down for padding points

        // Calculate ghost

    	link_up_scan.resize(padding_points+1);
    	link_up_scan.fill<0>(0);

    	ite = link_up_scan.getGPUIterator();

    	CUDA_LAUNCH((SparseGridGpuKernels::link_construct_up_count<dim,
    								BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
    								blockSize>),
    								ite,pd_points.toKernel(),grid_up.toKernel(),this->toKernel(),link_up_scan.toKernel(),p_up);

    	openfpm::scan((unsigned int *)link_up_scan.template getDeviceBuffer<0>(),link_up_scan.size(),(unsigned int *)link_up_scan.template getDeviceBuffer<0>(),context);

    	link_up_scan.template deviceToHost<0>(link_up_scan.size()-1,link_up_scan.size()-1);

    	size_t np_lup = link_up_scan.template get<0>(link_up_scan.size()-1);

    	link_up.resize(np_lup);

    	CUDA_LAUNCH((SparseGridGpuKernels::link_construct_insert_up<dim,
    								BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
    								blockSize>),ite,pd_points.toKernel(),grid_up.toKernel(),this->toKernel(),link_up_scan.toKernel(),link_up.toKernel(),p_up);

    	link_up_scan.resize(link_up_scan.size()-1);
    }

    /*! \Brief Before insert any element you have to call this function to initialize the insert buffer
     *
     * \param nBlock number of blocks the insert buffer has
     * \param nSlot maximum number of insertion each block does
     *
     */
    template<typename dim3T>
    void setGPUInsertBuffer(dim3T nBlock, dim3T nSlot)
    {
        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
        ::setGPUInsertBuffer(
                dim3SizeToInt(nBlock),
                dim3SizeToInt(nSlot)
        );
    }

	/*! \brief In case we manually set the added index buffer and the add data buffer we have to call this
	 *         function before flush
	 *
	 *
	 */
	void preFlush()
	{
		BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::preFlush();
	}

    template<typename stencil_type = NNStar<dim>, typename checker_type = No_check>
    void tagBoundaries(mgpu::ofp_context_t &context, checker_type chk = checker_type(), tag_boundaries opt = tag_boundaries::NO_CALCULATE_EXISTING_POINTS)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        const unsigned int dataChunkSize = BlockTypeOf<AggregateBlockT, 0>::size;
        unsigned int numScalars = indexBuffer.size() * dataChunkSize;

        if (numScalars == 0) return;
        if (findNN == false)
        {
        	findNeighbours<stencil_type>();
        	findNN = true;
        }

        // NOTE: Here we want to work only on one data chunk per block!

        unsigned int localThreadBlockSize = dataChunkSize;
        unsigned int threadGridSize = numScalars % dataChunkSize == 0
                                    ? numScalars / dataChunkSize
                                    : 1 + numScalars / dataChunkSize;

        constexpr unsigned int nLoop = UIntDivCeil<(IntPow<blockEdgeSize + 2, dim>::value - IntPow<blockEdgeSize, dim>::value), (blockSize * 1)>::value; // todo: This works only for stencilSupportSize==1
//        constexpr unsigned int nLoop = IntPow<blockEdgeSize + 2, dim>::value; // todo: This works only for stencilSupportSize==1

        if (stencilSupportRadius == 1)
        {
            CUDA_LAUNCH_DIM3((SparseGridGpuKernels::tagBoundaries<
                    dim,
                    1,
                    BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
                    stencil_type,
                    checker_type>),
                    threadGridSize, localThreadBlockSize,indexBuffer.toKernel(), dataBuffer.toKernel(), this->template toKernelNN<stencil_type::nNN, nLoop>(), nn_blocks.toKernel(),chk);
        }
        else if (stencilSupportRadius == 2)
        {
        	CUDA_LAUNCH_DIM3((SparseGridGpuKernels::tagBoundaries<
                    dim,
                    2,
                    BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
                    stencil_type,
                    checker_type>),
                    threadGridSize, localThreadBlockSize,indexBuffer.toKernel(), dataBuffer.toKernel(), this->template toKernelNN<stencil_type::nNN, nLoop>(), nn_blocks.toKernel(),chk);
        }
        else if (stencilSupportRadius == 0)
        {
        	CUDA_LAUNCH_DIM3((SparseGridGpuKernels::tagBoundaries<
                    dim,
                    0,
                    BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
                    stencil_type,
                    checker_type>),
                    threadGridSize, localThreadBlockSize,indexBuffer.toKernel(), dataBuffer.toKernel(), this->template toKernelNN<stencil_type::nNN, nLoop>(), nn_blocks.toKernel(),chk);
        }
        else
        {
            //todo: KABOOOOOOM!!!
            std::cout << __FILE__ << ":" << __LINE__ << " error: stencilSupportRadius supported only up to 2, passed: " << stencilSupportRadius << std::endl;

        }

        if (opt == tag_boundaries::CALCULATE_EXISTING_POINTS)
        {
        	// first we calculate the existing points
        	openfpm::vector_gpu<aggregate<indexT>> block_points;

        	block_points.resize(indexBuffer.size() + 1);
        	block_points.template get<0>(block_points.size()-1) = 0;
        	block_points.template hostToDevice<0>(block_points.size()-1,block_points.size()-1);

        	ite_gpu<1> ite;

        	ite.wthr.x = indexBuffer.size();
        	ite.wthr.y = 1;
        	ite.wthr.z = 1;
        	ite.thr.x = getBlockSize();
        	ite.thr.y = 1;
        	ite.thr.z = 1;

        	CUDA_LAUNCH((SparseGridGpuKernels::calc_exist_points<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>),
        				 ite,
        				 dataBuffer.toKernel(),
        				 block_points.toKernel());

        	// than we scan
        	openfpm::scan((indexT *)block_points.template getDeviceBuffer<0>(),block_points.size(),(indexT *)block_points.template getDeviceBuffer<0>(),context);

        	// Get the total number of points
        	block_points.template deviceToHost<0>(block_points.size()-1,block_points.size()-1);
        	size_t tot = block_points.template get<0>(block_points.size()-1);
        	e_points.resize(tot);

        	// we fill e_points
        	CUDA_LAUNCH((SparseGridGpuKernels::fill_e_points<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>),ite,
        				 dataBuffer.toKernel(),
        				 block_points.toKernel(),
        				 e_points.toKernel())

        }

        cudaDeviceSynchronize();
    }

    template<typename NNtype = NNStar<dim>>
    void findNeighbours()
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();

        const unsigned int numBlocks = indexBuffer.size();
        const unsigned int numScalars = numBlocks * NNtype::nNN;
        nn_blocks.resize(numScalars);

        if (numScalars == 0) return;

        // NOTE: Here we want to work only on one data chunk per block!

        unsigned int localThreadBlockSize = NNtype::nNN;

        unsigned int threadGridSize = numScalars % localThreadBlockSize == 0
                                      ? numScalars / localThreadBlockSize
                                      : 1 + numScalars / localThreadBlockSize;

        CUDA_LAUNCH_DIM3((SparseGridGpuKernels::findNeighbours<dim,NNtype>),
                threadGridSize, localThreadBlockSize,indexBuffer.toKernel(), this->toKernel(),nn_blocks.toKernel());

        findNN = true;
    }

    size_t countExistingElements() const
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        constexpr unsigned int pMask = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask;
        typedef typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT BAggregateT;
        typedef BlockTypeOf<BAggregateT, pMask> MaskBlockT;
        constexpr unsigned int blockSize = MaskBlockT::size;
        const auto bufferSize = indexBuffer.size();

        size_t numExistingElements = 0;

        for (size_t blockId=0; blockId<bufferSize; ++blockId)
        {
            auto dataBlock = dataBuffer.get(blockId); // Avoid binary searches as much as possible
            for (size_t elementId=0; elementId<blockSize; ++elementId)
            {
                const auto curMask = dataBlock.template get<pMask>()[elementId];

                if (this->exist(curMask))
                {
                    ++numExistingElements;
                }
            }
        }

        return numExistingElements;
    }

    size_t countBoundaryElements()
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        constexpr unsigned int pMask = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask;
        typedef typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT BAggregateT;
        typedef BlockTypeOf<BAggregateT, pMask> MaskBlockT;
        constexpr unsigned int blockSize = MaskBlockT::size;
        const auto bufferSize = indexBuffer.size();

        size_t numBoundaryElements = 0;

        for (size_t blockId=0; blockId<bufferSize; ++blockId)
        {
            auto dataBlock = dataBuffer.get(blockId); // Avoid binary searches as much as possible
            for (size_t elementId=0; elementId<blockSize; ++elementId)
            {
                const auto curMask = dataBlock.template get<pMask>()[elementId];

                if (this->exist(curMask) && this->isPadding(curMask))
                {
                    ++numBoundaryElements;
                }
            }
        }

        return numBoundaryElements;
    }

    // Also count mean+stdDev of occupancy of existing blocks
    void measureBlockOccupancyMemory(double &mean, double &deviation)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        constexpr unsigned int pMask = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask;
        typedef typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT BAggregateT;
        typedef BlockTypeOf<BAggregateT, pMask> MaskBlockT;
        constexpr unsigned int blockSize = MaskBlockT::size;
        const auto bufferSize = indexBuffer.size();

        openfpm::vector<double> measures;

        for (size_t blockId=0; blockId<bufferSize; ++blockId)
        {
            auto dataBlock = dataBuffer.get(blockId); // Avoid binary searches as much as possible
            size_t numElementsInBlock = 0;
            for (size_t elementId=0; elementId<blockSize; ++elementId)
            {
                const auto curMask = dataBlock.template get<pMask>()[elementId];

                if (this->exist(curMask))
                {
                    ++numElementsInBlock;
                }
            }
            double blockOccupancy = static_cast<double>(numElementsInBlock)/blockSize;
            measures.add(blockOccupancy);
        }

        standard_deviation(measures, mean, deviation);
    }

    // Also count mean+stdDev of occupancy of existing blocks
    void measureBlockOccupancy(double &mean, double &deviation)
    {
        // Here it is crucial to use "auto &" as the type, as we need to be sure to pass the reference to the actual buffers!
        auto & indexBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
        auto & dataBuffer = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();

        constexpr unsigned int pMask = BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask;
        typedef typename BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::AggregateInternalT BAggregateT;
        typedef BlockTypeOf<BAggregateT, pMask> MaskBlockT;
        constexpr unsigned int blockSize = MaskBlockT::size;
        const auto bufferSize = indexBuffer.size();

        openfpm::vector<double> measures;

        for (size_t blockId=0; blockId<bufferSize; ++blockId)
        {
            auto dataBlock = dataBuffer.get(blockId); // Avoid binary searches as much as possible
            size_t numElementsInBlock = 0;
            for (size_t elementId=0; elementId<blockSize; ++elementId)
            {
                const auto curMask = dataBlock.template get<pMask>()[elementId];

                if (this->exist(curMask) && !this->isPadding(curMask))
                {
                    ++numElementsInBlock;
                }
            }
            double blockOccupancy = static_cast<double>(numElementsInBlock)/blockSize;
            measures.add(blockOccupancy);
        }

        standard_deviation(measures, mean, deviation);
    }

    /*! \brief Apply a convolution using a cross like stencil
     *
     * in 2D for example the stencil is
     *
     \verbatim

         *
       * * *
         *

     \endverbatim
     *
     *
     */
	template<unsigned int prop_src, unsigned int prop_dst, unsigned int stencil_size, typename lambda_f, typename ... ArgsT >
	void conv_cross(grid_key_dx<3> start, grid_key_dx<3> stop , lambda_f func, ArgsT ... args)
	{
		Box<dim,int> box;

		for (int i = 0 ; i < dim ; i++)
		{
			box.setLow(i,start.get(i));
			box.setHigh(i,stop.get(i));
		}

		applyStencils< SparseGridGpuKernels::stencil_cross_func<dim,prop_src,prop_dst,stencil_size> >(box,STENCIL_MODE_INPLACE,func, args ...);
	}


    /*! \brief Apply a free type convolution using blocks
     *
     *
     */
	template<unsigned int prop_src, unsigned int prop_dst, unsigned int stencil_size, typename lambda_f, typename ... ArgsT >
	void conv(grid_key_dx<3> start, grid_key_dx<3> stop , lambda_f func, ArgsT ... args)
	{
		Box<dim,int> box;

		for (int i = 0 ; i < dim ; i++)
		{
			box.setLow(i,start.get(i));
			box.setHigh(i,stop.get(i));
		}

		constexpr unsigned int nLoop = UIntDivCeil<(IntPow<blockEdgeSize + 2, dim>::value), (blockSize)>::value;

		applyStencils< SparseGridGpuKernels::stencil_cross_func_conv<dim,nLoop,prop_src,prop_dst,stencil_size> >(box,STENCIL_MODE_INPLACE,func, args ...);
	}

    /*! \brief Apply a free type convolution using blocks
     *
     *
     */
	template<unsigned int prop_src1, unsigned int prop_src2, unsigned int prop_dst1 , unsigned int prop_dst2, unsigned int stencil_size, typename lambda_f, typename ... ArgsT >
	void conv2(grid_key_dx<dim> start, grid_key_dx<dim> stop , lambda_f func, ArgsT ... args)
	{
		Box<dim,int> box;

		for (int i = 0 ; i < dim ; i++)
		{
			box.setLow(i,start.get(i));
			box.setHigh(i,stop.get(i));
		}

        constexpr unsigned int nLoop = UIntDivCeil<(IntPow<blockEdgeSize + 2, dim>::value), (blockSize)>::value;

		applyStencils< SparseGridGpuKernels::stencil_func_conv2<dim,nLoop,prop_src1,prop_src2,prop_dst1,prop_dst2,stencil_size> >(box,STENCIL_MODE_INPLACE,func, args ...);
	}

	/*! \brief Return a Box with the  range if the SparseGrid
	 *
	 * \return the range of the grid
	 *
	 */
	Box<dim,int> getBox()
	{
		Box<dim,int> b;

		for (int i = 0 ; i < dim ; i++)
		{
			b.setLow(i,0);
			b.setHigh(i,gridGeometry.getSize()[i]);
		}

		return b;
	}

    //todo: Move implems into a functor for compile time choice of stencil mode
    template<typename stencil, typename... Args>
    void applyStencils(const Box<dim,int> & box, StencilMode mode, Args... args)
    {
        if (findNN == false)
        {
        	findNeighbours<typename stencil::stencil_type>();
        	findNN = true;
        }

        // Apply the given stencil on all elements which are not boundary-tagged
        // The idea is to have this function launch a __global__ function (made by us) on all existing blocks
        // then this kernel checks if elements exist && !padding and on those it calls the user-provided
        // __device__ functor. The mode of the stencil application is used basically to choose how we load the block
        // that we pass to the user functor as storeBlock: in case of Insert, we get the block through an insert (and
        // we also call the necessary aux functions); in case of an In-place we just get the block from the data buffer.
        switch (mode)
        {
            case STENCIL_MODE_INPLACE:
                applyStencilInPlace<stencil>(box,mode,args...);
                break;
            case STENCIL_MODE_INPLACE_NO_SHARED:
                applyStencilInPlaceNoShared<stencil>(box,mode,args...);
                break;
        }
    }
    template<typename stencil1, typename stencil2, typename ... otherStencils, typename... Args>
    void applyStencils(Box<dim,int> box, StencilMode mode, Args... args)
    {
        applyStencils<stencil1>(box,mode, args...);
        applyStencils<stencil2, otherStencils ...>(box,mode, args...);
    }

    template<typename BitMaskT>
    inline static bool isPadding(BitMaskT &bitMask)
    {
        return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
        ::getBit(bitMask, PADDING_BIT);
    }

    template<typename BitMaskT>
    inline static void setPadding(BitMaskT &bitMask)
    {
        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
        ::setBit(bitMask, PADDING_BIT);
    }

    template<typename BitMaskT>
    inline static void unsetPadding(BitMaskT &bitMask)
    {
        BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>
        ::unsetBit(bitMask, PADDING_BIT);
    }

    /*! \brief Linearization of  block coordinates
     *
     * \param blockCoord block coordinates
     *
     * \return the linearized index
     *
     */
    template<typename CoordT>
    inline size_t getBlockLinId(const CoordT & blockCoord) const
    {
        return gridGeometry.BlockLinId(blockCoord);
    }

    /*! \brief Insert the point on host side and flush directly
     *
     * First you have to move everything on host with deviceToHost, insertFlush and than move to GPU again
     *
     * \param grid point where to insert
     *
     * \return a reference to the data to fill
     *
     *
     */
    template<unsigned int p, typename CoordT>
    auto insertFlush(const grid_key_dx<dim,CoordT> &coord) -> ScalarTypeOf<AggregateBlockT, p> &
    {
    	// Linearized block_id
    	auto lin = gridGeometry.LinId(coord);
    	indexT block_id = lin / blockSize;
    	indexT local_id = lin % blockSize;

    	typedef BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base> BMG;

    	auto block_data = this->insertBlockFlush(block_id);
    	block_data.template get<BMG::pMask>()[local_id] = 1;

    	return block_data.template get<p>()[local_id];
    }

    template<unsigned int p>
    void print_vct_add_data()
    {
    	typedef BlockMapGpu<
    	        typename aggregate_convert<dim,blockEdgeSize,AggregateT>::type,
    	        threadBlockSize, indexT, layout_base> BMG;

    	auto & bM = BMG::blockMap.private_get_vct_add_data();
    	auto & vI = BMG::blockMap.private_get_vct_add_index();
    	bM.template deviceToHost<p>();
    	vI.template deviceToHost<0>();

    	std::cout << "vct_add_data: " << std::endl;

    	for (size_t i = 0 ; i < bM.size() ; i++)
    	{
    		std::cout << i << "  index: " << vI.template get<0>(i) << "  BlockData: " << std::endl;
    		for (size_t j = 0 ; j < blockSize ; j++)
    		{
    			std::cout << (int)bM.template get<p>(i)[j] << "  ";
    		}

    		std::cout << std::endl;
    	}
    }

    /*! \brief set the background for property p
     *
     * \tparam p property p
     *
     */
    template<unsigned int p>
    void setBackgroundValue(ScalarTypeOf<AggregateBlockT, p> backgroundValue)
    {
        bck.template get<p>() = backgroundValue;

        BMG::template setBackgroundValue<p>(backgroundValue);
    }

    /////////////////////////////////// DISTRIBUTED INTERFACE ///////////////////////

    /*! \brief memory requested to pack this object
     *
     * \param req request
     *
     */
	template<int ... prp> inline
	void packRequest(size_t & req, mgpu::ofp_context_t &context) const
    {
    	ite_gpu<1> ite;

    	auto & indexBuffer = private_get_index_array();
    	auto & dataBuffer = private_get_data_array();

    	ite.wthr.x = indexBuffer.size();
    	ite.wthr.y = 1;
    	ite.wthr.z = 1;
    	ite.thr.x = getBlockSize();
    	ite.thr.y = 1;
    	ite.thr.z = 1;

    	tmp.resize(indexBuffer.size() + 1);

		// Launch a kernel that count the number of element on each chunks
    	CUDA_LAUNCH((SparseGridGpuKernels::calc_exist_points<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>),
    				 ite,
    				 dataBuffer.toKernel(),
    				 tmp.toKernel());

    	openfpm::scan((indexT *)tmp. template getDeviceBuffer<0>(),
    						tmp.size(), (indexT *)tmp. template getDeviceBuffer<0>(), context);

    	tmp.template deviceToHost<0>(tmp.size()-1,tmp.size()-1);

    	sparsegridgpu_pack_request<AggregateT,prp ...> spq;

    	boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                (prp)>>(spq);

		size_t n_pnt = tmp.template get<0>(tmp.size()-1);


		                        // 4 byte each chunks        data                      // we use short to pack the offset
		                        // for each counter
		req = sizeof(indexT) +               // byte required to pack the number
			  sizeof(indexT)*indexBuffer.size() +    // byte required to pack the chunk indexes
			  sizeof(indexT)*tmp.size() +            // byte required to pack the scan of the chunks points
			  n_pnt*(spq.point_size + sizeof(short int) + sizeof(unsigned char));    // byte required to pack data + offset position
    }

	/*! \brief Calculate the size to pack part of this structure
	 *
	 * \warning this function does not return the byte immediately but it buffer the request until you call
	 *          packCalculate
	 *
	 * \warning to reset call packReset()
	 *
	 * \tparam prp set of properties to pack
	 *

	 * \param sub sub-grid to pack
	 * \param req output byte
	 *
	 */
	template<int ... prp> inline
	void packRequest(SparseGridGpu_iterator_sub<dim,self> & sub_it,
					 size_t & req) const
	{
		pack_subs.add();

		for (int i = 0 ; i < dim ; i++)
		{
			pack_subs.template get<0>(pack_subs.size()-1)[i] = sub_it.getStart().get(i);
			pack_subs.template get<1>(pack_subs.size()-1)[i] = sub_it.getStop().get(i);
		}
	}

	/*! \brief Reset the pack calculation
	 *
	 *
	 */
	void packReset()
	{
		pack_subs.clear();

	    index_ptrs.clear();
	    scan_ptrs.clear();
	    data_ptrs.clear();
	    offset_ptrs.clear();
	    mask_ptrs.clear();

		req_index = 0;
	}

	/*! \brief Calculate the size of the information to pack
	 *
	 * \param req output size (it does not reset the counter it accumulate)
	 * \param context gpu contect
	 *
	 */
	template<int ... prp> inline
	void packCalculate(size_t & req, mgpu::ofp_context_t &context)
	{
    	ite_gpu<1> ite;
		pack_subs.template hostToDevice<0,1>();

    	auto & indexBuffer = private_get_index_array();
    	auto & dataBuffer = private_get_data_array();

    	ite.wthr.x = indexBuffer.size();
    	ite.wthr.y = 1;
    	ite.wthr.z = 1;
    	ite.thr.x = getBlockSize();
    	ite.thr.y = 1;
    	ite.thr.z = 1;

    	tmp.resize((indexBuffer.size() + 1)*pack_subs.size());

    	if (indexBuffer.size() != 0)
    	{
			if (pack_subs.size() <= 32)
			{
				// Launch a kernel that count the number of element on each chunks
				CUDA_LAUNCH((SparseGridGpuKernels::calc_exist_points_with_boxes<dim,
																			BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
																			32,
																			indexT>),
						 ite,
						 indexBuffer.toKernel(),
						 pack_subs.toKernel(),
						 gridGeometry,
						 dataBuffer.toKernel(),
						 tmp.toKernel(),
						 indexBuffer.size() + 1);
			}
			else if (pack_subs.size() <= 64)
			{
				// Launch a kernel that count the number of element on each chunks
				CUDA_LAUNCH((SparseGridGpuKernels::calc_exist_points_with_boxes<dim,
																			BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
																			64,
																			indexT>),
						 ite,
						 indexBuffer.toKernel(),
						 pack_subs.toKernel(),
						 gridGeometry,
						 dataBuffer.toKernel(),
						 tmp.toKernel(),
						 indexBuffer.size() + 1);
			}
			else if (pack_subs.size() <= 96)
			{
				// Launch a kernel that count the number of element on each chunks
				CUDA_LAUNCH((SparseGridGpuKernels::calc_exist_points_with_boxes<dim,
																			BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
																			96,
																			indexT>),
						 ite,
						 indexBuffer.toKernel(),
						 pack_subs.toKernel(),
						 gridGeometry,
						 dataBuffer.toKernel(),
						 tmp.toKernel(),
						 indexBuffer.size() + 1);
			}
			else if (pack_subs.size() <= 128)
			{
				// Launch a kernel that count the number of element on each chunks
				CUDA_LAUNCH((SparseGridGpuKernels::calc_exist_points_with_boxes<dim,
																			BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
																			128,
																			indexT>),
						 ite,
						 indexBuffer.toKernel(),
						 pack_subs.toKernel(),
						 gridGeometry,
						 dataBuffer.toKernel(),
						 tmp.toKernel(),
						 indexBuffer.size() + 1);
			}
			else
			{
				std::cout << __FILE__ << ":" << __LINE__ << " error no implementation available of packCalculate, create a new case for " << pack_subs.size() << std::endl;
			}
    	}

    	sparsegridgpu_pack_request<AggregateT,prp ...> spq;

    	boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prp)>>(spq);

    	scan_it.resize(pack_subs.size());

    	// scan all
    	for (size_t i = 0 ; i < pack_subs.size() ; i++)
    	{
        	size_t n_pnt = 0;
        	size_t n_cnk = 0;

    		tmp.template get<0>((i+1)*(indexBuffer.size() + 1)-1) = 0;
    		tmp.template get<1>((i+1)*(indexBuffer.size() + 1)-1) = 0;

    		// put a zero at the end
    		tmp.template hostToDevice<0>((i+1)*(indexBuffer.size() + 1)-1,(i+1)*(indexBuffer.size() + 1)-1);
    		tmp.template hostToDevice<1>((i+1)*(indexBuffer.size() + 1)-1,(i+1)*(indexBuffer.size() + 1)-1);

    		openfpm::scan(((indexT *)tmp. template getDeviceBuffer<0>()) + i*(indexBuffer.size() + 1),
    						indexBuffer.size() + 1, (indexT *)tmp. template getDeviceBuffer<0>() + i*(indexBuffer.size() + 1), context);

    		openfpm::scan(((unsigned int *)tmp. template getDeviceBuffer<1>()) + i*(indexBuffer.size() + 1),
    						indexBuffer.size() + 1, (unsigned int *)tmp. template getDeviceBuffer<1>() + i*(indexBuffer.size() + 1), context);

    		tmp.template deviceToHost<0>((i+1)*(indexBuffer.size() + 1)-1,(i+1)*(indexBuffer.size() + 1)-1);
    		tmp.template deviceToHost<1>((i+1)*(indexBuffer.size() + 1)-1,(i+1)*(indexBuffer.size() + 1)-1);

    		scan_it.template get<0>(i) = tmp.template get<0>((i+1)*(indexBuffer.size() + 1)-1);

    		n_pnt = tmp.template get<0>((i+1)*(indexBuffer.size() + 1)-1);
    		n_cnk = tmp.template get<1>((i+1)*(indexBuffer.size() + 1)-1);

    		req += sizeof(indexT) +               // byte required to pack the number of chunk packed
    				2*dim*sizeof(int) +           // starting point + size of the indexing packing
    				  sizeof(indexT)*n_cnk +    				   // byte required to pack the chunk indexes
    				  align_number(sizeof(indexT),(n_cnk+1)*sizeof(unsigned int)) +            // byte required to pack the scan of the chunk point
    				  align_number(sizeof(indexT),n_pnt*(spq.point_size)) +  // byte required to pack data
    				  align_number(sizeof(indexT),n_pnt*sizeof(short int)) + // byte required to pack offsets
    				  align_number(sizeof(indexT),n_pnt*sizeof(unsigned char));  // byte required to pack masks
    	}

    	scan_it.template hostToDevice<0>();

		openfpm::scan((indexT *)scan_it. template getDeviceBuffer<0>(),
								scan_it.size(), (indexT *)scan_it. template getDeviceBuffer<0>(), context);
	}

	/*! \brief Return the mapping vector used to know where the data has been added
	 *
	 * \return the mapping vector
	 *
	 */
	auto getMappingVector() -> decltype(this->blockMap.getMappingVector())
	{
		return this->blockMap.getMappingVector();
	}

	/*! \brief Return the mapping vector used to know where the data has been added
	 *
	 * \return the mapping vector
	 *
	 */
	auto getMergeIndexMapVector() -> decltype(this->blockMap.getMergeIndexMapVector())
	{
		return this->blockMap.getMergeIndexMapVector();
	}

	/*! \brief Pack the object into the memory given an iterator
	 *
	 * \warning the pack does not happen here but it happen in packFinalize, the request is just queued
	 *
	 * \tparam prp properties to pack
	 *
	 * \param mem preallocated memory where to pack the objects
	 * \param sub_it sub grid iterator ( or the elements in the grid to pack )
	 *        \warning in this case this parameter is ignored pack must be called with the same sequence of
	 *                 packRequest
	 *
	 * \param sts pack statistic
	 *
	 */
	template<int ... prp> void pack(ExtPreAlloc<CudaMemory> & mem,
									SparseGridGpu_iterator_sub<dim,self> & sub_it,
									Pack_stat & sts)
	{
		unsigned int i = req_index;

    	sparsegridgpu_pack_request<AggregateT,prp ...> spq;
    	boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prp)>>(spq);

    	auto & indexBuffer = private_get_index_array();
    	auto & dataBuffer = private_get_data_array();

		size_t n_pnt = tmp.template get<0>((i+1)*(indexBuffer.size() + 1)-1);
		size_t n_cnk = tmp.template get<1>((i+1)*(indexBuffer.size() + 1)-1);

		Packer<size_t,CudaMemory>::pack(mem,n_cnk,sts);
		mem.hostToDevice(mem.getOffset(),mem.getOffset()+sizeof(size_t));

		size_t offset1 = mem.getOffsetEnd();

		grid_key_dx<dim> key = sub_it.getStart();

		for (int i = 0 ; i < dim ; i++)
		{Packer<int,CudaMemory>::pack(mem,key.get(i),sts);}

		for (int i = 0 ; i < dim ; i++)
		{Packer<int,CudaMemory>::pack(mem,(int)gridGeometry.getSize()[i],sts);}

		mem.hostToDevice(offset1,offset1+2*dim*sizeof(int));

		// chunk indexes
		mem.allocate(n_cnk*sizeof(indexT));
		index_ptrs.add(mem.getDevicePointer());

		// chunk point scan
		mem.allocate( align_number(sizeof(indexT),(n_cnk+1)*sizeof(unsigned int)) );
		scan_ptrs.add(mem.getDevicePointer());

		// chunk data
		mem.allocate( align_number(sizeof(indexT),n_pnt*(spq.point_size)) );
		data_ptrs.add(mem.getDevicePointer());

		// space for offsets
		mem.allocate( align_number(sizeof(indexT),n_pnt*sizeof(short int) ) );
		offset_ptrs.add(mem.getDevicePointer());

		// space for offsets
		mem.allocate( align_number(sizeof(indexT),n_pnt*sizeof(unsigned char) ) );
		mask_ptrs.add(mem.getDevicePointer());

		req_index++;
	}


	/*! \brief It finalize the queued operations of remove() and copy_to()
	 *
	 * \warning you have to call this function with opt rem_copy_opt::PHASE1 than rem_copy_opt::PHASE2
	 *          and finally rem_copy_opt::PHASE3
	 *
	 * In particular suppose we have 3 grids and we copy sections of the grid 1 into 2 and 3 , sections of grid 2 into
	 * 1 and 3 and sections of grid 3 into 1 and 2. Than we have to first call on all 3 grids removeCopyFinalize with
	 * PHASE1. Than we call removeCopyFinalize with PHASE2 on all 3 grids and finally removeCopyFinalize with PHASE3
	 *
	 * In case we can guarantee that the chunks structure has not changed we can pass the option KEEP_GEOMETRY to make
	 * this function faster
	 *
	 * \param ctx context
	 * \param options
	 *
	 */
	template<unsigned int ... prp>
	void removeCopyToFinalize(mgpu::ofp_context_t & ctx, int opt)
	{
		if ((opt & 0x3) == rem_copy_opt::PHASE1)
		{
			this->template removeCopyToFinalize_phase1<prp ...>(ctx,opt);
		}
		else if ((opt & 0x3) == rem_copy_opt::PHASE2)
		{
			this->template removeCopyToFinalize_phase2<prp ...>(ctx,opt);
		}
		else
		{
			this->template removeCopyToFinalize_phase3<prp ...>(ctx,opt,false);
		}
	}

	/*! \brief Finalize the packing procedure
	 *
	 * \tparam prp properties to pack
	 *
	 * \param mem preallocated memory where to pack the objects
	 * \param sub_it sub grid iterator ( or the elements in the grid to pack )
	 *        \warning in this case this parameter is ignored pack must be called with the same sequence of
	 *                 packRequest
	 *
	 * \param sts pack statistic
	 *
	 */
	template<int ... prp> void packFinalize(ExtPreAlloc<CudaMemory> & mem,
									Pack_stat & sts,
									int opt = 0,
									bool is_pack_remote = false)
	{

    	RestorePackVariableIfKeepGeometry(opt,is_pack_remote);

    	if (pack_subs.size() <= 32)
    	{
    		pack_sg_implement<32,prp...>(mem,sts,opt,is_pack_remote);
    	}
    	else if (pack_subs.size() <= 64)
    	{
    		pack_sg_implement<64, prp...>(mem,sts,opt,is_pack_remote);
    	}
    	else
    	{
    		std::cout << __FILE__ << ":" << __LINE__ << " error no implementation available of packCalculate, create a new case for " << pack_subs.size() << std::endl;
    	}

    	savePackVariableIfNotKeepGeometry(opt,is_pack_remote);
	}

	/*! \brief In this case it does nothing
	 *
	 * \note this function exist to respect the interface to work as distributed
	 *
	 */
	void removeAddUnpackReset()
	{
		rem_sects.clear();

    	auto & vad = BMG::blockMap.private_get_vct_add_data();
    	auto & vai = BMG::blockMap.private_get_vct_add_index();

    	vad.clear();
    	vai.clear();

		// Clear variables
		offset_ptrs_cp.clear();
		scan_ptrs_cp.clear();
		n_cnk_cp.clear();
		n_pnt_cp.clear();
		data_base_ptr_cp.clear();
		box_cp.clear();
		n_shifts_cp.clear();
		convert_blk.clear();
		tmp2.clear();
	}

	/*! \brief Remove the points we queues to remove
	 *
	 * \see
	 *
	 * \param context modern gpu context
	 *
	 */
	void removePoints(mgpu::ofp_context_t& context)
	{
    	auto & indexBuffer = private_get_index_array();
    	auto & dataBuffer = private_get_data_array();

		// first we remove
		if (rem_sects.size() != 0)
		{
			rem_sects.template hostToDevice<0,1>();

			tmp.resize(indexBuffer.size() + 1);

			tmp.template get<1>(tmp.size()-1) = 0;
			tmp.template hostToDevice<1>(tmp.size()-1,tmp.size()-1);

			auto ite = indexBuffer.getGPUIterator();

			if (has_work_gpu(ite) == true)
			{
				// mark all the chunks that must remove points
				CUDA_LAUNCH((SparseGridGpuKernels::calc_remove_points_chunks_boxes<dim,
															 BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask,
															 blockEdgeSize>),ite,indexBuffer.toKernel(),rem_sects.toKernel(),
																  gridGeometry,dataBuffer.toKernel(),
																  tmp.toKernel());

				// scan
				openfpm::scan((unsigned int *)tmp.template getDeviceBuffer<1>(),tmp.size(),(unsigned int *)tmp.template getDeviceBuffer<1>(),context);

				tmp.template deviceToHost<1>(tmp.size()-1,tmp.size()-1);

				// get the number of chunks involved
				size_t nr_cnk = tmp.template get<1>(tmp.size()-1);

				tmp3.resize(nr_cnk);

				// collect the chunks involved in the remove
				ite = indexBuffer.getGPUIterator();

				if (has_work_gpu(ite) == false)	{return;}

				CUDA_LAUNCH((SparseGridGpuKernels::collect_rem_chunks),ite,tmp.toKernel(),tmp3.toKernel());

				// Launch to remove points

				ite = tmp3.getGPUIterator();

				ite.wthr.x = tmp3.size();
				ite.wthr.y = 1;
				ite.wthr.z = 1;
				ite.thr.x = getBlockSize();
				ite.thr.y = 1;
				ite.thr.z = 1;

				if (has_work_gpu(ite) == false)	{return;}

				CUDA_LAUNCH((SparseGridGpuKernels::remove_points<dim,
																BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>),
																ite,indexBuffer.toKernel(),
																gridGeometry,
																dataBuffer.toKernel(),
																tmp3.toKernel(),
																rem_sects.toKernel());

				tmp3.clear();
			}
		}
	}

	/*! \brief This function remove the points we queue to remove and it flush all the added/unpacked data
	 *
	 * \note this function exist to respect the interface to work as distributed
	 *
	 */
	template<unsigned int ... prp>
	void removeAddUnpackFinalize(mgpu::ofp_context_t& context, int opt)
	{
    	removePoints(context);

		removeCopyToFinalize_phase3<prp ...>(context,opt,true);
	}

	/*! \brief Reset the queue to remove and copy section of grids
	 *
	 *
	 */
	void copyRemoveReset()
	{
		rem_sects.clear();
		copySect.clear();
		offset_ptrs_cp.clear();
		scan_ptrs_cp.clear();
		data_base_ptr_cp.clear();
		n_cnk_cp.clear();
		n_pnt_cp.clear();
	    n_shifts_cp.clear();
		convert_blk.clear();
		box_cp.clear();
		data_base_ptr_cp.clear();

		tmp2.clear();
	}

	/*! \brief Remove all the points in this region
	 *
	 * \warning does not remove the chunks only the points
	 *
	 * \param box_src box to kill the points
	 *
	 */
	void remove(const Box<dim,int> & section_to_delete)
	{
		rem_sects.add(section_to_delete);
	}

	/*! \brief This is a multiresolution sparse grid so is a compressed format
	 *
	 * \return true
	 *
	 */
	static constexpr bool isCompressed()
	{
		return true;
	}

	/*! \brief It queue a copy
	 *
	 * \param grid_src source grid
	 * \param box_src gource box
	 * \param box_dst destination box
	 *
	 */
	void copy_to(self & grid_src,
		         const Box<dim,size_t> & box_src,
			     const Box<dim,size_t> & box_dst)
	{
		// first we launch a kernel to count the number of points we have

		sparse_grid_section<self> sgs(*this,box_src,box_dst);

		grid_src.copySect.add(sgs);
	}

	/*! \brief unpack the sub-grid object
	 *
	 * \tparam prp properties to unpack
	 *
	 * \param mem preallocated memory from where to unpack the object
	 * \param sub sub-grid iterator
	 * \param obj object where to unpack
	 *
	 */
	template<unsigned int ... prp, typename S2>
	void unpack(ExtPreAlloc<S2> & mem,
				SparseGridGpu_iterator_sub<dim,self> & sub_it,
				Unpack_stat & ps,
				mgpu::ofp_context_t &context,
				rem_copy_opt opt = rem_copy_opt::NONE_OPT)
	{
		////////////////////////////////////////////////////////////

		if ((opt & rem_copy_opt::KEEP_GEOMETRY) == false)
		{
			this->template addAndConvertPackedChunkToTmp<prp ...>(mem,sub_it,ps,context);

			// readjust mem
		}
		else
		{
			// we have to increment ps by the right amount
	    	sparsegridgpu_pack_request<AggregateT,prp ...> spq;
	    	boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prp)>>(spq);

			// First get the number of chunks

			size_t n_cnk;

			// Unpack the number of chunks
			mem.deviceToHost(ps.getOffset(),ps.getOffset() + sizeof(size_t) + 2*dim*sizeof(int));
			Unpacker<size_t,S2>::unpack(mem,n_cnk,ps);

			// Unpack origin of the chunk indexing
			for (int i = 0 ; i < dim ; i++)
			{
				int tmp;
				Unpacker<int,S2>::unpack(mem,tmp,ps);
			}

			for (int i = 0 ; i < dim ; i++)
			{
				int tmp;
				Unpacker<int,S2>::unpack(mem,tmp,ps);
			}

			size_t actual_offset = n_cnk*sizeof(indexT);
			unsigned int * scan = (unsigned int *)((unsigned char *)mem.getDevicePointer() + ps.getOffset() + n_cnk*sizeof(indexT));

			mem.deviceToHost(ps.getOffset() + actual_offset + n_cnk*sizeof(unsigned int),
							 ps.getOffset() + actual_offset + n_cnk*sizeof(unsigned int) + sizeof(unsigned int));

			// Unpack number of points
			// calculate the number of total points
			size_t n_pnt = *(unsigned int *)((unsigned char *)mem.getPointer() + ps.getOffset() + actual_offset + n_cnk*sizeof(unsigned int));
			actual_offset += align_number(sizeof(indexT),(n_cnk+1)*sizeof(unsigned int));

			void * data_base_ptr = (void *)((unsigned char *)mem.getDevicePointer() + ps.getOffset() + actual_offset );

			actual_offset += align_number(sizeof(indexT),n_pnt*(spq.point_size));
			short int * offsets = (short int *)((unsigned char *)mem.getDevicePointer() + ps.getOffset() + actual_offset);

			actual_offset += align_number(sizeof(indexT),n_pnt*sizeof(short));
			actual_offset += align_number(sizeof(indexT),n_pnt*sizeof(unsigned char));

			scan_ptrs_cp.add(scan);
			offset_ptrs_cp.add(offsets);
			data_base_ptr_cp.add(data_base_ptr);

			ps.addOffset(actual_offset);
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*! \brief Return a SparseGrid iterator
     *
     * \return a SparseGrid iterator
     *
     */
    decltype(self::type_of_iterator()) getIterator() const
    {
    	return decltype(self::type_of_iterator())(*this);
    }

    /*! \brief Return a SparseGrid iterator only on a sub-set of elements
     *
     * \return a SparseGrid iterator on a subset of elements
     *
     */
    decltype(self::type_of_subiterator()) getIterator(const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop, int is_to_init = 1) const
    {
    	return decltype(self::type_of_subiterator())(*this,start,stop,is_to_init);
    }

    /*! \brief Return the index array of the blocks
     *
     * \return the index arrays of the blocks
     *
     */
    auto private_get_add_index_array() -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.private_get_vct_add_index()) &
    {
    	return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.private_get_vct_add_index();
    }

    /*! \brief Return the index array of the blocks
     *
     * \return the index arrays of the blocks
     *
     */
    auto private_get_add_index_array() const -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.private_get_vct_add_index()) &
    {
    	return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.private_get_vct_add_index();
    }

    /*! \brief Return the index array of the blocks
     *
     * \return the index arrays of the blocks
     *
     */
    auto private_get_index_array() const -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer()) &
    {
    	return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
    }

	auto getSegmentToOutMap() -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getSegmentToOutMap())
	{
		return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getSegmentToOutMap();
	}

	auto getSegmentToOutMap() const -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getSegmentToOutMap())
	{
		return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getSegmentToOutMap();
	}

	auto getSegmentToMergeIndexMap() -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getSegmentToMergeIndexMap())
	{
		return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getSegmentToMergeIndexMap();
	}

	auto getSegmentToMergeIndexMap() const -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getSegmentToMergeIndexMap())
	{
		return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getSegmentToMergeIndexMap();
	}

    /*! \brief Return the index array of the blocks
     *
     * \return the index arrays of the blocks
     *
     */
    auto private_get_data_array() -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer()) &
    {
    	return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();
    }

    /*! \brief Return the index array of the blocks
     *
     * \return the index arrays of the blocks
     *
     */
    auto private_get_index_array() -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer()) &
    {
    	return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getIndexBuffer();
    }

    /*! \brief Return the index array of the blocks
     *
     * \return the index arrays of the blocks
     *
     */
    auto private_get_data_array() const -> decltype(BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer()) &
    {
    	return BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::blockMap.getDataBuffer();
    }

    /*! \brief Return the index array of the blocks
     *
     * \return the index arrays of the blocks
     *
     */
    auto private_get_neighborhood_array() -> decltype(nn_blocks) &
    {
    	return nn_blocks;
    }


#if defined(OPENFPM_DATA_ENABLE_IO_MODULE) || defined(PERFORMANCE_TEST) || defined(VTKWRITER_HPP_)

	/*! \brief write the sparse grid into VTK
	 *
	 * \param out VTK output
	 *
	 */
	template<typename Tw = float> bool write(const std::string & output)
	{
		Point<dim,double> spacing;
		Point<dim,double> offset;

		spacing.one();
		offset.zero();

		return write_with_spacing_offset(output,spacing,offset);
	}

	/*! \brief write the sparse grid into VTK
	 *
	 * \param out VTK output
	 *
	 */
	template<typename Tw = float>
	bool write_with_spacing_offset(const std::string & output, Point<dim,double> spacing, Point<dim,double> offset)
	{
		file_type ft = file_type::BINARY;

		auto & bm = this->private_get_blockMap();

		auto & index = bm.getIndexBuffer();
		auto & data = bm.getDataBuffer();

		openfpm::vector<Point<dim,Tw>> tmp_pos;
		openfpm::vector<typename aggregate_add<AggregateT>::type> tmp_prp;

		// copy position and properties

		auto it = index.getIterator();

		while(it.isNext())
		{
			auto key = it.get();

			Point<dim,Tw> p;

			for (size_t i = 0 ; i < gridGeometry.getBlockSize() ; i++)
			{
				if (data.template get<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>(key)[i] != 0)
				{
					// Get the block index
					grid_key_dx<dim,int> keyg = gridGeometry.InvLinId(index.template get<0>(key),i);

					for (size_t k = 0 ; k < dim ; k++)
					{p.get(k) = keyg.get(k)*spacing[k] + offset[k]*spacing[k];}

					tmp_pos.add(p);

					tmp_prp.add();
					copy_prop_to_vector_block<decltype(data.get_o(key)),decltype(tmp_prp.last())>
					cp(data.get_o(key),tmp_prp.last(),key,i);

					boost::mpl::for_each_ref< boost::mpl::range_c<int,0,AggregateT::max_prop> >(cp);

					tmp_prp.last().template get<AggregateT::max_prop>() = data.template get<BlockMapGpu<AggregateInternalT, threadBlockSize, indexT, layout_base>::pMask>(key)[i];
				}
			}

			++it;
		}

		// VTKWriter for a set of points
		VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim,Tw>>, openfpm::vector<typename aggregate_add<AggregateT>::type>>, VECTOR_POINTS> vtk_writer;
		vtk_writer.add(tmp_pos,tmp_prp,tmp_pos.size());

		openfpm::vector<std::string> prp_names;

		// Write the VTK file
		return vtk_writer.write(output,prp_names,"sparse_grid","",ft);
	}

	/*! \brief write the sparse grid into VTK
	 *
	 * \param out VTK output
	 *
	 */
	template<typename Tw = float> bool write_debug(const std::string & output, Point<dim,double> spacing, Point<dim,double> offset)
	{
		//! subdomains_X.vtk domain for the local processor (X) as union of sub-domain
		VTKWriter<openfpm::vector<SpaceBox<dim, double>>, VECTOR_BOX> vtk_box1;

		openfpm::vector<SpaceBox<dim,double>> chunks_box;

		auto & ids = private_get_index_array();

		fill_chunks_boxes(chunks_box,ids,spacing,offset);

		vtk_box1.add(chunks_box);
		vtk_box1.write(std::string("chunks_") + output + std::string(".vtk"));

		//write data

		write_with_spacing_offset(std::string("data_") + output + std::string(".vtk"),spacing,offset);

		return true;
	}

#endif
};

template<unsigned int dim,
		 typename AggregateT,
		 unsigned int blockEdgeSize = default_edge<dim>::type::value,
		 unsigned int threadBlockSize = 128,
		 typename indexT=long int,
		 template<typename> class layout_base=memory_traits_inte,
		 typename linearizer = grid_zmb<dim, blockEdgeSize>>
using SparseGridGpu_z = SparseGridGpu<dim,AggregateT,blockEdgeSize,threadBlockSize,indexT,layout_base,linearizer>;

#endif //OPENFPM_PDATA_SPARSEGRIDGPU_HPP
