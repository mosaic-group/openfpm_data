/*
 * SparseGrid_iterator_block.hpp
 *
 *  Created on: Feb 25, 2020
 *      Author: i-bird
 */

#ifndef SPARSEGRID_ITERATOR_BLOCK_HPP_
#define SPARSEGRID_ITERATOR_BLOCK_HPP_

#include "Grid/iterators/grid_skin_iterator.hpp"
#include "SparseGrid_chunk_copy.hpp"


template<int c,bool is_neg = c < 0>
struct fix_neg_to_one
{
	typedef boost::mpl::int_<c> type;
};

template<int c>
struct fix_neg_to_one<c,true>
{
	typedef boost::mpl::int_<1> type;
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to set a grid info
 *
 * \tparam grid_sm_type
 * \tparam vector_blocks_ext
 *
 */
template<unsigned int dim, typename vector_blocks_ext>
struct calc_loc
{
	grid_key_dx<dim> & k;

	/*! \brief constructor
	 *
	 * \param v set of pointer buffers to set
	 *
	 */
	calc_loc(grid_key_dx<dim> & k)
	:k(k)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& val)
	{
		k.set_d(T::value,boost::mpl::at<typename vector_blocks_ext::type,boost::mpl::int_<T::value>>::type::value*k.get(T::value));
	}
};



/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to set a grid info
 *
 * \tparam grid_sm_type
 * \tparam vector_blocks_ext
 *
 */
template<unsigned int dim, typename header_type, typename vector_blocks_ext>
struct fill_chunk_block
{
	//! sizes
	Box<dim,size_t> cnk_box;

	unsigned int chunk_id;

	header_type & header;

	/*! \brief constructor
	 *
	 * \param v set of pointer buffers to set
	 *
	 */
	inline fill_chunk_block(header_type & header, unsigned int chunk_id)
	:chunk_id(chunk_id),header(header)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& val)
	{
		cnk_box.setLow(T::value,header.get(chunk_id).pos.get(T::value));
		cnk_box.setHigh(T::value,header.get(chunk_id).pos.get(T::value) + boost::mpl::at<typename vector_blocks_ext::type,T>::type::value - 1);
	}
};

template<unsigned int prop, unsigned int stencil_size, unsigned int dim, typename vector_blocks_exts, typename vector_ext>
struct loadBlock_impl
{
	template<unsigned int N1, typename T, typename SparseGridType>
	static void loadBlock(T arr[N1], SparseGridType & sgt, int chunk_id, unsigned char mask[N1])
	{
		get_block_sizes<dim,stencil_size,vector_blocks_exts,vector_ext> gbs;

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,dim> >(gbs);

		grid_sm<dim,void> g_block(gbs.sz_ext);
		grid_sm<dim,void> g_in_block(gbs.sz_block);
		grid_sm<dim,void> g_block_arr(gbs.sz_tot);
;
		grid_key_dx_iterator<dim> it_in_block(g_in_block);

		auto & data = sgt.private_get_data();
		auto & header_mask = sgt.private_get_header_mask();

		auto & h = header_mask.get(chunk_id);

		auto & ref_block = data.template get<prop>(chunk_id);

		while(it_in_block.isNext())
		{
			auto p = it_in_block.get();

			grid_key_dx<dim> arr_p;

			for (int i = 0 ; i < dim ; i++)
			{arr_p.set_d(i,p.get(i)+stencil_size);}

			size_t id = g_block_arr.LinId(arr_p);
			size_t idi = g_in_block.LinId(p);

			arr[id] = ref_block[idi];
			mask[id] = exist_sub(h,idi);

			++it_in_block;
		}
	}

	template<unsigned int N1, typename T, typename SparseGridType>
	static void loadBlock(T arr[N1], SparseGridType & sgt, int chunk_id)
	{
		get_block_sizes<dim,stencil_size,vector_blocks_exts,vector_ext> gbs;

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,dim> >(gbs);

		grid_sm<dim,void> g_block(gbs.sz_ext);
		grid_sm<dim,void> g_in_block(gbs.sz_block);
		grid_sm<dim,void> g_block_arr(gbs.sz_tot);

		grid_key_dx_iterator<dim> it_in_block(g_in_block);

		auto & data = sgt.private_get_data();
		auto & header_mask = sgt.private_get_header_mask();

		auto & ref_block = data.template get<prop>(chunk_id);

		while(it_in_block.isNext())
		{
			auto p = it_in_block.get();

			grid_key_dx<dim> arr_p;

			for (int i = 0 ; i < dim ; i++)
			{arr_p.set_d(i,p.get(i)+stencil_size);}

			arr[g_block_arr.LinId(arr_p)] = ref_block[g_in_block.LinId(p)];

			++it_in_block;
		}
	}

	template<unsigned int N1, typename T, typename SparseGridType>
	static void storeBlock(T arr[N1], SparseGridType & sgt, int chunk_id)
	{
		get_block_sizes<dim,stencil_size,vector_blocks_exts,vector_ext> gbs;

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,dim> >(gbs);

		grid_sm<dim,void> g_block(gbs.sz_ext);
		grid_sm<dim,void> g_in_block(gbs.sz_block);
		grid_sm<dim,void> g_block_arr(gbs.sz_tot);

		grid_key_dx_iterator<dim> it_in_block(g_in_block);

		auto & data = sgt.private_get_data();
		auto & header_mask = sgt.private_get_header_mask();

		auto & ref_block = data.template get<prop>(chunk_id);

		while(it_in_block.isNext())
		{
			auto p = it_in_block.get();

			ref_block[g_in_block.LinId(p)] = arr[g_in_block.LinId(p)];

			++it_in_block;
		}
	}


	/*! \brief load the border
	 *
	 *
	 *
	 */
	template<unsigned int N1, typename T, typename SparseGridType,typename backgroundType>
	static void loadBorder(T arr[N1],
			         SparseGridType & sgt,
			         size_t chunk_id,
			         openfpm::vector<unsigned int> & bord,
			         openfpm::vector<grid_key_dx<SparseGridType::dims>> & block_skin,
			         openfpm::vector<unsigned int> & chunk_ids ,
			         openfpm::vector<short int> & offsets,
			         unsigned char mask[N1],
			         backgroundType & background,
			         openfpm::vector<unsigned int> & maps_blk)
	{
		typedef typename generate_array_vector<size_t,typename vector_blocks_exts::type>::result size;

		auto & data = sgt.private_get_data();
		auto & header_mask = sgt.private_get_header_mask();
		auto & header_inf = sgt.private_get_header_inf();

		auto & hm = header_mask.get(chunk_id);
		auto & hc = header_inf.get(chunk_id);

		maps_blk.resize(block_skin.size());

		for (int i = 0 ; i < maps_blk.size() ; i++)
		{
			grid_key_dx<dim> p;

			for (int j = 0 ; j < dim ; j++)
			{p.set_d(j,block_skin.get(i).get(j) + hc.pos.get(j) / size::data[j] - 1);}

			maps_blk.get(i) = sgt.getChunk(p);
		}

		for (int i = 0 ; i < bord.size(); i++)
		{
			size_t ac = maps_blk.get(chunk_ids.get(i));

			size_t b = bord.get(i);
			size_t off = offsets.template get<0>(i);

			auto & h = header_mask.get(ac);

			arr[b] = (ac == data.size()-1)?background.template get<prop>():data.template get<prop>(ac)[off];
			mask[b] = (ac == data.size()-1)?0:exist_sub(h,off);
		}
	}
};

/*! \brief optimized 3d version
 *
 */
template<unsigned int prop, unsigned int stencil_size, typename vector_blocks_exts, typename vector_ext>
struct loadBlock_impl<prop,stencil_size,3,vector_blocks_exts,vector_ext>
{
	template<unsigned int N1, typename T, typename SparseGridType>
	inline static void loadBlock(T arr[N1], SparseGridType & sgt, int chunk_id, unsigned char mask[N1])
	{
		auto & data = sgt.private_get_data();
		auto & header_mask = sgt.private_get_header_mask();

		auto & h = header_mask.get(chunk_id);

		// Faster version

		auto & chunk = data.template get<prop>(chunk_id);

		copy_xyz<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,false>::template copy<N1>(arr,mask,h,chunk);
	}


	template<unsigned int N1, typename T, typename SparseGridType>
	inline static void storeBlock(T arr[N1], SparseGridType & sgt, int chunk_id)
	{

		auto & data = sgt.private_get_data();
		auto & header_mask = sgt.private_get_header_mask();

		// Faster version

		auto & chunk = data.template get<prop>(chunk_id);

		copy_xyz<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,false>::template store<N1>(arr,chunk);

	}



	/*! \brief load the border
	 *
	 *
	 *
	 */
	template<bool findNN, typename NNType, unsigned int N1, typename T, typename SparseGridType,typename backgroundType>
	inline static void loadBorder(T arr[N1],
			         SparseGridType & sgt,
			         size_t chunk_id,
			         openfpm::vector<unsigned int> & bord,
			         openfpm::vector<grid_key_dx<SparseGridType::dims>> & block_skin,
			         openfpm::vector<unsigned int> & chunk_ids ,
			         openfpm::vector<short int> & offsets,
			         unsigned char mask[N1],
			         backgroundType & background,
			         openfpm::vector<unsigned int> & maps_blk)
	{
		typedef typename generate_array_vector<size_t,typename vector_blocks_exts::type>::result size;

		auto & data = sgt.private_get_data();
		auto & header_mask = sgt.private_get_header_mask();
		auto & NNlist = sgt.private_get_nnlist();

		auto & h = header_mask.get(chunk_id);


		typedef typename generate_array_vector<size_t,typename vector_blocks_exts::type>::result size;

		typedef typename boost::mpl::at<typename vector_blocks_exts::type,boost::mpl::int_<0>>::type sz0;
		typedef typename boost::mpl::at<typename vector_blocks_exts::type,boost::mpl::int_<1>>::type sz1;
		typedef typename boost::mpl::at<typename vector_blocks_exts::type,boost::mpl::int_<2>>::type sz2;

		grid_key_dx<3> p;

		bool exist;
		long int r;
		if (findNN == false)
		{
			p = sgt.getChunkPos(chunk_id) + grid_key_dx<3>({0,0,1});
			r = sgt.getChunk(p,exist);
			NNlist.template get<0>(chunk_id*NNType::nNN) = (exist)?r:-1;
		}
		else
		{
			r = NNlist.template get<0>(chunk_id*NNType::nNN);
			exist = (r != -1);
		}
		if (exist == true)
		{
			auto & h = header_mask.get(r);
			copy_xy_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template copy<0,stencil_size+sz2::value,N1>(arr,mask,h,data.get(r));
		}
		else
		{
			copy_xy_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template mask_null<stencil_size+sz2::value,N1>(mask);
		}
		if (findNN == false)
		{
			p = sgt.getChunkPos(chunk_id) + grid_key_dx<3>({0,0,-1});
			r = sgt.getChunk(p,exist);
			NNlist.template get<0>(chunk_id*NNType::nNN+1) = (exist)?r:-1;
		}
		else
		{
			r = NNlist.template get<0>(chunk_id*NNType::nNN+1);
			exist = (r != -1);
		}
		if (exist == true)
		{
			auto & h = header_mask.get(r);
			copy_xy_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template copy<sz2::value - stencil_size,0,N1>(arr,mask,h,data.get(r));
		}
		else
		{
			copy_xy_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template mask_null<0,N1>(mask);
		}

		if (findNN == false)
		{
			p = sgt.getChunkPos(chunk_id) + grid_key_dx<3>({0,1,0});
			r = sgt.getChunk(p,exist);
			NNlist.template get<0>(chunk_id*NNType::nNN+2) = (exist)?r:-1;
		}
		else
		{
			r = NNlist.template get<0>(chunk_id*NNType::nNN+2);
			exist = (r != -1);
		}
		if (exist == true)
		{
			auto & h = header_mask.get(r);
			copy_xz_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template copy<0,stencil_size+sz1::value,N1>(arr,mask,h,data.get(r));
		}
		else
		{
			copy_xz_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template mask_null<stencil_size+sz1::value,N1>(mask);
		}
		if (findNN == false)
		{
			p = sgt.getChunkPos(chunk_id) + grid_key_dx<3>({0,-1,0});
			r = sgt.getChunk(p,exist);
			NNlist.template get<0>(chunk_id*NNType::nNN+3) = (exist)?r:-1;
		}
		else
		{
			r = NNlist.template get<0>(chunk_id*NNType::nNN+3);
			exist = (r != -1);
		}
		if (exist == true)
		{
			auto & h = header_mask.get(r);
			copy_xz_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template copy<sz1::value-stencil_size,0,N1>(arr,mask,h,data.get(r));
		}
		else
		{
			copy_xz_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template mask_null<0,N1>(mask);
		}

		if (findNN == false)
		{
			p = sgt.getChunkPos(chunk_id) + grid_key_dx<3>({1,0,0});
			r = sgt.getChunk(p,exist);
			NNlist.template get<0>(chunk_id*NNType::nNN+4) = (exist)?r:-1;
		}
		else
		{
			r = NNlist.template get<0>(chunk_id*NNType::nNN+4);
			exist = (r != -1);
		}
		if (exist == true)
		{
			auto & h = header_mask.get(r);
			copy_yz_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template copy<0,sz0::value+stencil_size,N1>(arr,mask,h,data.get(r));
		}
		else
		{
			copy_yz_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template mask_null<sz0::value+stencil_size,N1>(mask);
		}
		if (findNN == false)
		{
			p = sgt.getChunkPos(chunk_id) + grid_key_dx<3>({-1,0,0});
			r = sgt.getChunk(p,exist);
			NNlist.template get<0>(chunk_id*NNType::nNN+5) = (exist)?r:-1;
		}
		else
		{
			r = NNlist.template get<0>(chunk_id*NNType::nNN+5);
			exist = (r != -1);
		}
		if (exist == true)
		{
			auto & h = header_mask.get(r);
			copy_yz_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template copy<sz0::value-stencil_size,0,N1>(arr,mask,h,data.get(r));
		}
		else
		{
			copy_yz_3<is_layout_inte<typename SparseGridType::memory_traits >::type::value,prop,stencil_size,typename vector_blocks_exts::type,NNType::is_cross>::template mask_null<0,N1>(mask);
		}
	}
};


/*! \brief Grid key sparse iterator on a sub-part of the domain
 *
 *
 */
template<unsigned dim,
		 unsigned int stencil_size,
		 typename SparseGridType,
		 typename vector_blocks_exts,
		 typename vector_ext = typename vmpl_create_constant<dim,1>::type>
class grid_key_sparse_dx_iterator_block_sub
{
	//! SparseGrid
	SparseGridType & spg;

	//! point to the actual chunk
	size_t chunk_id;

	//! Starting point
	grid_key_dx<dim> start_;

	//! Stop point
	grid_key_dx<dim> stop_;

	//! Sub-grid box
	Box<dim,size_t> bx;

	//! border
    openfpm::vector<unsigned int> bord;

    //! chunks ids
    openfpm::vector<unsigned int> chunk_shifts;

    //! offsets
    openfpm::vector<short int> offsets;

    //! blocks skin
    openfpm::vector<grid_key_dx<dim>> block_skin;

    // chunk header container
    mheader<SparseGridType::chunking_type::size::value> * hm;
    cheader<dim> * hc;

    // temporary buffer for Load border
    openfpm::vector<unsigned int> maps_blk;

    //! background value
    typename SparseGridType::background_type & background;

    //!iteration block
    Box<dim,size_t> block_it;

	/*! \brief Everytime we move to a new chunk we calculate on which indexes we have to iterate
	 *
	 *
	 */
	void SelectValid()
	{
		auto & header = spg.private_get_header_inf();
		auto & header_mask = spg.private_get_header_mask();

		while (chunk_id < header.size())
		{
			auto & mask = header_mask.get(chunk_id).mask;

			fill_chunk_block<dim,decltype(header),vector_blocks_exts> fcb(header,chunk_id);

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,dim>>(fcb);

			if (bx.Intersect(fcb.cnk_box,block_it) == true)
			{
				block_it -= header.get(chunk_id).pos.toPoint();
				break;
			}
			else
			{chunk_id += 1;}
		}
	}

public:

	// we create first a vector with

	typedef typename vmpl_sum_constant<2*stencil_size,typename vector_blocks_exts::type>::type stop_border_vmpl;
	typedef typename vmpl_create_constant<dim,stencil_size>::type start_border_vmpl;

	typedef typename generate_array_vector<size_t,typename vector_blocks_exts::type>::result size;
	typedef typename generate_array_vector<size_t,start_border_vmpl>::result start_b_;
	typedef typename generate_array_vector<size_t,stop_border_vmpl>::result stop_b_;

	static const int sizeBlock = vector_blocks_exts::size::value;
	static const int sizeBlockBord = vmpl_reduce_prod<stop_border_vmpl>::type::value;

	/*! \brief Default constructor
	 *
	 * \warning extremely unsafe
	 * If you use this constructor before use the iterator you should call reinitialize first
	 *
	 */
	grid_key_sparse_dx_iterator_block_sub()	{};

	grid_key_sparse_dx_iterator_block_sub(SparseGridType & spg,
								const grid_key_dx<dim> & start,
								const grid_key_dx<dim> & stop,
								typename SparseGridType::background_type & background)
	:spg(spg),chunk_id(1),
	 start_(start),stop_(stop),background(background)
	{
		// Create border coeficents
		get_block_sizes<dim,stencil_size,vector_blocks_exts,vector_ext> gbs;

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,dim> >(gbs);

		Box<dim,int> skinb;
		Box<dim,int> skinbb;

		size_t bc[dim];
		for (int i = 0 ; i < dim ; i ++)
		{
			skinb.setLow(i,0);
			skinb.setHigh(i,gbs.sz_tot[i]-1);
			skinbb.setLow(i,0);
			skinbb.setHigh(i,gbs.sz_ext_b[i]-1);
			bc[i] = NON_PERIODIC;
		}

		grid_sm<dim,void> g_smb(gbs.sz_ext_b);

		// Create block skin index

		openfpm::vector<unsigned int> b_map;
		grid_skin_iterator_bc<3> gsi_b(g_smb,skinbb,skinbb,bc);

		b_map.resize(g_smb.size());

		while (gsi_b.isNext())
		{
			auto p = gsi_b.get();

			block_skin.add(p);

			b_map.get(g_smb.LinId(p)) = block_skin.size() - 1;

			++gsi_b;
		}

		grid_sm<dim,void> g_sm(gbs.sz_tot);
		grid_skin_iterator_bc<3> gsi(g_sm,skinb,skinb,bc);

		while (gsi.isNext())
		{
			auto p = gsi.get();

			grid_key_dx<dim> sh;

			if (bord.size() == 1109)
			{
				int debug = 0;
				debug++;
			}

			bord.add(g_sm.LinId(p));

			short offset = 0;
			int stride = 1;
			for (int i = 0 ; i < dim ; i++)
			{

				if (p.get(i) < stencil_size)
				{offset += (gbs.sz_block[i]-1)*stride;}
				else if (p.get(i) >= gbs.sz_tot[i] - stencil_size)
				{offset += 0;}
				else
				{offset += (p.get(i)-stencil_size)*stride;}

				sh.set_d(i,(p.get(i) + (gbs.sz_block[i] - stencil_size)) / gbs.sz_block[i]);
				stride *= gbs.sz_block[i];
			}

			offsets.add(offset);

			size_t bid = g_smb.LinId(sh);
			chunk_shifts.add(b_map.get(bid));

			++gsi;
		}

		for (size_t i = 0 ; i < dim ; i++)
		{
			bx.setLow(i,start.get(i));
			bx.setHigh(i,stop.get(i));
		}

		SelectValid();
	}

	/*! \brief Reinitialize the iterator
	 *
	 * it re-initialize the iterator with the passed grid_key_dx_iterator_sub
	 * the actual position of the grid_key_dx_iterator_sub is ignored
	 *
	 * \param g_s_it grid_key_dx_iterator_sub
	 *
	 */
	inline void reinitialize(const grid_key_sparse_dx_iterator_sub<dim,vector_blocks_exts::size::value> & g_s_it)
	{
		spg = g_s_it.spg;
		chunk_id = g_s_it.chunk_id;
		start_ = g_s_it.start_;
		stop_ = g_s_it.stop_;
		bx = g_s_it.bx;
	}

	inline grid_key_sparse_dx_iterator_block_sub<dim,stencil_size,SparseGridType,vector_blocks_exts> & operator++()
	{
		auto & header = spg.private_get_header_inf();

		chunk_id++;

		if (chunk_id < header.size())
		{
			SelectValid();
		}

		return *this;
	}

	/*! \brief Return true if there is a next grid point
	 *
	 * \return true if there is the next grid point
	 *
	 */
	bool isNext()
	{
		auto & header = spg.private_get_header_inf();

		return chunk_id < header.size();
	}

	/*! \brief Return the starting point for the iteration
	 *
	 * \return the starting point
	 *
	 */
	const grid_key_dx<dim> & getStart() const
	{
		return start_;
	}

	/*! \brief Return the stop point for the iteration
	 *
	 * \return the stop point
	 *
	 */
	const grid_key_dx<dim> & getStop() const
	{
		return stop_;
	}


	template<unsigned int prop, typename T>
	void loadBlock(T arr[sizeBlock])
	{
		auto & header_mask = spg.private_get_header_mask();
		auto & header_inf = spg.private_get_header_inf();

		loadBlock_impl<prop,stencil_size,dim,vector_blocks_exts,vector_ext>::template loadBlock<prop>(arr,spg,chunk_id);

		hm = &header_mask.get(chunk_id);
		hc = &header_inf.get(chunk_id);
	}

	template<unsigned int prop,typename T>
	void loadBlock(T arr[sizeBlock], unsigned char mask[sizeBlock])
	{
		auto & header_mask = spg.private_get_header_mask();
		auto & header_inf = spg.private_get_header_inf();

		loadBlock_impl<prop,stencil_size,dim,vector_blocks_exts,vector_ext>::template loadBlock<prop>(arr,spg,chunk_id,mask);

		hm = &header_mask.get(chunk_id);
		hc = &header_inf.get(chunk_id);
	}

	template<unsigned int prop,typename T>
	void storeBlock(T arr[sizeBlock])
	{
		auto & header_mask = spg.private_get_header_mask();
		auto & header_inf = spg.private_get_header_inf();

		loadBlock_impl<prop,stencil_size,dim,vector_blocks_exts,vector_ext>::template storeBlock<prop>(arr,spg,chunk_id);

		hm = &header_mask.get(chunk_id);
		hc = &header_inf.get(chunk_id);
	}


	template<unsigned int prop, typename NNtype, bool findNN, typename T>
	void loadBlockBorder(T arr[sizeBlockBord],unsigned char mask[sizeBlockBord])
	{
		auto & header_mask = spg.private_get_header_mask();
		auto & header_inf = spg.private_get_header_inf();

		loadBlock_impl<prop,stencil_size,dim,vector_blocks_exts,vector_ext>::template loadBlock<sizeBlockBord>(arr,spg,chunk_id,mask);
		loadBlock_impl<prop,stencil_size,dim,vector_blocks_exts,vector_ext>::template loadBorder<findNN,NNtype,sizeBlockBord>(arr,spg,chunk_id,bord,block_skin,chunk_shifts,offsets,mask,background,maps_blk);

		hm = &header_mask.get(chunk_id);
		hc = &header_inf.get(chunk_id);
	}


	/*! \brief starting point of the computation block
	 *
	 * \param i coordinate
	 *
	 */
	constexpr int start_b(int i) const
	{
		return block_it.getLow(i) + stencil_size;
	}

	/*! \brief stopping point of the computation block
	 *
	 * \param i coordinate
	 *
	 */
	constexpr int stop_b(int i) const
	{
		return block_it.getHigh(i) + 1 + stencil_size;
	}

	/*! \brief starting point of the computation block
	 *
	 * \param i coordinate
	 *
	 */
	constexpr int start(int i) const
	{
		return block_it.getLow(i);
	}

	/*! \brief stopping point of the computation block
	 *
	 * \param i coordinate
	 *
	 */
	constexpr int stop(int i) const
	{
		return block_it.getHigh(i) + 1;
	}

	/*! \brief linearize an arbitrary set of index
	 *
	 * linearize an arbitrary set of index
	 *
	 */
	template<typename a, typename ...lT>
	__device__ __host__ inline size_t Lin(a v,lT...t) const
	{
#ifdef SE_CLASS1
		if (sizeof...(t)+1 > dim)
		{
			std::cerr << "Error incorrect grid cannot linearize more index than its dimensionality" << "\n";
		}
#endif

		return v*vmpl_reduce_prod_stop<typename vector_blocks_exts::type,(int)dim - (int)sizeof...(t) - 2>::type::value + Lin(t...);
	}

	//! Linearize a set of index
	template<typename a> __device__ __host__ inline size_t Lin(a v) const
	{
		return v*vmpl_reduce_prod_stop<typename vector_blocks_exts::type,(int)dim - 2>::type::value;
	}

	/*! \brief linearize an arbitrary set of index
	 *
	 * linearize an arbitrary set of index
	 *
	 */
	template<typename a, typename ...lT>
	__device__ __host__ inline size_t LinB(a v,lT...t) const
	{
#ifdef SE_CLASS1
		if (sizeof...(t)+1 > dim)
		{
			std::cerr << "Error incorrect grid cannot linearize more index than its dimensionality" << "\n";
		}
#endif

		return v*vmpl_reduce_prod_stop<stop_border_vmpl,(int)dim - (int)sizeof...(t) - 2>::type::value + LinB(t...);
	}

	//! Linearize a set of index
	template<typename a> __device__ __host__ inline size_t LinB(a v) const
	{
		return v*vmpl_reduce_prod_stop<stop_border_vmpl,(int)dim - 2>::type::value;
	}

	/*! \brief linearize an arbitrary set of index
	 *
	 * linearize an arbitrary set of index
	 *
	 */
	template<typename a, typename ...lT>
	__device__ __host__ inline size_t LinB_off(a v,lT...t) const
	{
#ifdef SE_CLASS1
		if (sizeof...(t)+1 > dim)
		{
			std::cerr << "Error incorrect grid cannot linearize more index than its dimensionality" << "\n";
		}
#endif

		return (v-stencil_size)*vmpl_reduce_prod_stop<typename vector_blocks_exts::type,(int)dim - (int)sizeof...(t) - 2>::type::value + LinB_off(t...);
	}

	//! Linearize a set of index
	template<typename a> __device__ __host__ inline size_t LinB_off(a v) const
	{
		return (v-stencil_size)*(vmpl_reduce_prod_stop<typename vector_blocks_exts::type,(int)dim - 2>::type::value);
	}

	/*! Check if the point exist
	 *
	 * \param args index to linearize
	 *
	 */
	template<typename ... ArgsType>
	bool exist(ArgsType ... args)
	{
		size_t l = LinB_off(args ...);

		return spg.exist_sub(*hm,l);
	}

	/*! \brief Return the chunk id
	 *
	 * \return the chunk id
	 *
	 */
	int getChunkId()
	{
		return chunk_id;
	}
};


#endif /* SPARSEGRID_ITERATOR_BLOCK_HPP_ */
