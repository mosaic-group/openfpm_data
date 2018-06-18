/*
 * SparseGrid.hpp
 *
 *  Created on: Oct 22, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_HPP_
#define OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_HPP_

#include "memory_ly/memory_array.hpp"
#include "memory_ly/memory_c.hpp"
#include "memory_ly/memory_conf.hpp"
#include "hash_map/hopscotch_map.h"
#include "hash_map/hopscotch_set.h"
#include "Vector/map_vector.hpp"
#include "util/variadic_to_vmpl.hpp"
#include "data_type/aggregate.hpp"
#include "SparseGridUtil.hpp"
#include "SparseGrid_iterator.hpp"
// We do not want parallel writer
#define NO_PARALLEL
#include "VTKWriter/VTKWriter.hpp"


template<typename Tsrc,typename Tdst>
class copy_prop_to_vector
{
	//! source
	Tsrc src;

	//! destination
	Tdst dst;

	size_t pos;

public:

	copy_prop_to_vector(Tsrc src, Tdst dst,size_t pos)
	:src(src),dst(dst),pos(pos)
	{}

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename std::remove_reference<decltype(dst.template get<T::value>())>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(src.template get<T::value>()[pos],dst.template get<T::value>());
	}

};


template<unsigned int dim, typename Tsrc,typename Tdst>
class copy_sparse_to_sparse
{
	//! source
	const Tsrc & src;

	//! destination
	Tdst & dst;

	//! source position
	grid_key_dx<dim> & pos_src;

	//! destination position
	grid_key_dx<dim> & pos_dst;

public:

	copy_sparse_to_sparse(const Tsrc & src, Tdst & dst,grid_key_dx<dim> & pos_src, grid_key_dx<dim> & pos_dst)
	:src(src),dst(dst),pos_src(pos_src),pos_dst(pos_dst)
	{}

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename std::remove_reference<decltype(dst.template insert<T::value>(pos_dst))>::type copy_rtype;

		meta_copy<copy_rtype>::meta_copy_(src.template get<T::value>(pos_src),dst.template insert<T::value>(pos_dst));
	}

};

template< template<typename,typename> class op,unsigned int dim, typename Tsrc,typename Tdst, unsigned int ... prp>
class copy_sparse_to_sparse_op
{
	//! source
	const Tsrc & src;

	//! destination
	Tdst & dst;

	//! source position
	grid_key_dx<dim> & pos_src;

	//! destination position
	grid_key_dx<dim> & pos_dst;

	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

public:

	copy_sparse_to_sparse_op(const Tsrc & src, Tdst & dst,grid_key_dx<dim> & pos_src, grid_key_dx<dim> & pos_dst)
	:src(src),dst(dst),pos_src(pos_src),pos_dst(pos_dst)
	{}

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type idx_type;
		typedef typename std::remove_reference<decltype(dst.template insert<idx_type::value>(pos_dst))>::type copy_rtype;

		if (dst.existPoint(pos_dst) == false)
		{meta_copy_op<replace_,copy_rtype>::meta_copy_op_(src.template get<idx_type::value>(pos_src),dst.template insert<idx_type::value>(pos_dst));}
		else
		{meta_copy_op<op,copy_rtype>::meta_copy_op_(src.template get<idx_type::value>(pos_src),dst.template insert<idx_type::value>(pos_dst));}
	}

};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy a boost::mpl::vector into runtime array
 *
 */

template<unsigned int dim, typename mpl_v>
struct copy_sz
{
	//! sz site_t
	size_t (& sz)[dim];


	/*! \brief constructor
	 *
	 * \param sz runtime sz to fill
	 *
	 */
	inline copy_sz(size_t (& sz)[dim])
	:sz(sz)
	{
	};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		sz[T::value] = boost::mpl::at<mpl_v,boost::mpl::int_<T::value>>::type::value;
	}
};

#define BIT_SHIFT_SIZE_T 6

template<unsigned int dim, typename T, typename S, typename chunking>
class sgrid_cpu<dim,T,S,typename memory_traits_lin<T>::type,chunking>
{
	//! cache pointer
	mutable size_t cache_pnt;

	//! background values
	T background;

	//! cache
	mutable long int cache[SGRID_CACHE];

	//! cached id
	mutable long int cached_id[SGRID_CACHE];

	//! Map to convert from grid coordinates to chunk
	tsl::hopscotch_map<size_t, size_t> map;

	//! indicate which element in the chunk are really filled
	openfpm::vector<cheader<dim,chunking::size::value>> header;

	//Definition of the chunks
	typedef typename v_transform_two<Ft_chunk,boost::mpl::int_<chunking::size::value>,typename T::type>::type chunk_def;

	//! vector of chunks
	openfpm::vector<aggregate_bfv<chunk_def>> chunks;

	//! grid size information
	grid_sm<dim,void> g_sm;

	//! grid size information with shift
	grid_sm<dim,void> g_sm_shift;

	//! conversion position in the chunks
	grid_key_dx<dim> pos_chunk[chunking::size::value];

	//! size of the chunk
	size_t sz_cnk[dim];

	openfpm::vector<size_t> empty_v;

	/*! \brief Remove
	 *
	 *
	 */
	template<unsigned int n_ele>
	inline void remove_from_chunk(size_t sub_id,
			 	 	 	 	 	  size_t & nele,
								  size_t (& mask)[n_ele / (sizeof(size_t)*8) + (n_ele % (sizeof(size_t)*8) != 0) + 1])
	{
		// set the mask to null
		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		size_t swt = mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check;

		nele = (swt)?nele-1:nele;

		mask[sub_id >> BIT_SHIFT_SIZE_T] &= ~mask_check;
	}

	/*! \brief reconstruct the map
	 *
	 *
	 *
	 */
	inline void reconstruct_map()
	{
		// reconstruct map

		map.clear();
		for (size_t i = 0 ; i < header.size() ; i++)
		{
			grid_key_dx<dim> kh = header.get(i).pos;
			grid_key_dx<dim> kl;

			// shift the key
			key_shift<dim,chunking>::shift(kh,kl);

			long int lin_id = g_sm_shift.LinId(kh);

			map[lin_id] = i;
		}
	}

	/*! \brief Eliminate empty chunks
	 *
	 * \warning Because this operation is time consuming it perform the operation once
	 *          we reach a critical size in the list of the empty chunks
	 *
	 */
	inline void remove_empty()
	{
		if (empty_v.size() >= FLUSH_REMOVE)
		{
			// eliminate double entry

			empty_v.sort();
			empty_v.unique();

			// Because chunks can be refilled the empty list can contain chunks that are
			// filled so before remove we have to check that they are really empty

			for (int i = empty_v.size() - 1 ; i >= 0  ; i--)
			{
				if (header.get(empty_v.get(i)).nele != 0)
				{empty_v.remove(i);}
			}

			header.remove(empty_v);
			chunks.remove(empty_v);

			// reconstruct map

			reconstruct_map();

			empty_v.clear();

			// cache must be cleared

			clear_cache();
		}
	}

	/*! \brief add on cache
	 *
	 * \param lin_id linearized id
	 * \param active_cnk active chunk
	 *
	 */
	inline void add_on_cache(size_t lin_id, size_t active_cnk) const
	{
		// Add on cache the chunk
		cache[cache_pnt] = lin_id;
		cached_id[cache_pnt] = active_cnk;
		cache_pnt++;
		cache_pnt = (cache_pnt >= SGRID_CACHE)?0:cache_pnt;
	}

	/*! \brief reset the cache
	 *
	 *
	 */
	inline void clear_cache()
	{
		cache_pnt = 0;
		for (size_t i = 0 ; i < SGRID_CACHE ; i++)
		{cache[i] = -1;}
	}

	/*! \brief set the grid shift from size
	 *
	 * \param sz size of the grid
	 * \param g_sm_shift chunks grid size
	 *
	 */
	void set_g_shift_from_size(const size_t (& sz)[dim], grid_sm<dim,void> & g_sm_shift)
	{
		grid_key_dx<dim> cs;
		grid_key_dx<dim> unused;

		for (size_t i = 0 ; i < dim ; i++)
		{cs.set_d(i,sz[i]);}

		key_shift<dim,chunking>::shift(cs,unused);

		size_t sz_i[dim];

		for (size_t i = 0 ; i < dim ; i++)
		{sz_i[i] = cs.get(i) + 1;}

		g_sm_shift.setDimensions(sz_i);
	}

	/*! \brief initialize
	 *
	 *
	 */
	void init()
	{
		for (size_t i = 0 ; i < SGRID_CACHE ; i++)
		{cache[i] = -1;}

		// fill pos_g

		copy_sz<dim,typename chunking::type> cpsz(sz_cnk);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,dim> >(cpsz);

		grid_sm<dim,void> gs(sz_cnk);

		grid_key_dx_iterator<dim> it(gs);
		size_t cnt = 0;

		while (it.isNext())
		{
			auto key = it.get();

			for (size_t i = 0 ; i < dim ; i++)
			{
				pos_chunk[cnt].set_d(i,key.get(i));
			}

			++cnt;
			++it;
		}
	}

	/*! \brief Before insert data you have to do this
	 *
	 * \param v1 grid key where you want to insert data
	 * \param active_cnk output which chunk is competent for this data
	 * \param sub_id element inside the chunk
	 *
	 */
	inline bool pre_insert(const grid_key_dx<dim> & v1, size_t & active_cnk, size_t & sub_id)
	{
		bool exist = true;
		active_cnk = 0;

		grid_key_dx<dim> kh = v1;
		grid_key_dx<dim> kl;

		// shift the key
		key_shift<dim,chunking>::shift(kh,kl);

		long int lin_id = g_sm_shift.LinId(kh);

		size_t id = 0;
		for (size_t k = 0 ; k < SGRID_CACHE; k++)
		{id += (cache[k] == lin_id)?k+1:0;}

		if (id == 0)
		{
			// we do not have it in cache we check if we have it in the map

			auto fnd = map.find(lin_id);
			if (fnd == map.end())
			{
				// we do not have it in the map create a chunk

				map[lin_id] = chunks.size();
				chunks.add();
				header.add();
				header.last().pos = kh;
				header.last().nele = 0;

				// set the mask to null
				auto & h = header.last().mask;

				for (size_t i = 0 ; i < chunking::size::value / (sizeof(size_t)*8) + (chunking::size::value % (sizeof(size_t)*8) != 0) + 1 ; i++)
				{h[i] = 0;}

				key_shift<dim,chunking>::cpos(header.last().pos);

				active_cnk = chunks.size() - 1;
			}
			else
			{
				// we have it in the map

				active_cnk = fnd->second;
			}

			// Add on cache the chunk
			cache[cache_pnt] = lin_id;
			cached_id[cache_pnt] = active_cnk;
			cache_pnt++;
			cache_pnt = (cache_pnt >= SGRID_CACHE)?0:cache_pnt;
		}
		else
		{
			active_cnk = cached_id[id-1];
			cache_pnt = id;
			cache_pnt = (cache_pnt == SGRID_CACHE)?0:cache_pnt;
		}

		sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

		// the chunk is in cache, solve

		// we notify that we added one element
		auto & h = header.get(active_cnk);

		// we set the mask

		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		exist = h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check;
		h.nele = (exist)?h.nele:h.nele + 1;
		h.mask[sub_id >> BIT_SHIFT_SIZE_T] |= mask_check;

		return exist;
	}

	inline void remove_point(const grid_key_dx<dim> & v1)
	{
		size_t active_cnk = 0;

		grid_key_dx<dim> kh = v1;
		grid_key_dx<dim> kl;

		// shift the key
		key_shift<dim,chunking>::shift(kh,kl);

		long int lin_id = g_sm_shift.LinId(kh);

		size_t id = 0;
		for (size_t k = 0 ; k < SGRID_CACHE; k++)
		{id += (cache[k] == lin_id)?k+1:0;}

		if (id == 0)
		{
			// we do not have it in cache we check if we have it in the map

			auto fnd = map.find(lin_id);
			if (fnd == map.end())
			{return;}
			else
			{active_cnk = fnd->second;}

			// Add on cache the chunk
			cache[cache_pnt] = lin_id;
			cached_id[cache_pnt] = active_cnk;
			cache_pnt++;
			cache_pnt = (cache_pnt >= SGRID_CACHE)?0:cache_pnt;
		}
		else
		{
			active_cnk = cached_id[id-1];
			cache_pnt = id;
			cache_pnt = (cache_pnt == SGRID_CACHE)?0:cache_pnt;
		}

		size_t sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

		// eliminate the element

		// set the mask to null
		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		auto & h = header.get(active_cnk);
		size_t swt = h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check;

		h.nele = (swt)?h.nele-1:h.nele;

		h.mask[sub_id >> BIT_SHIFT_SIZE_T] &= ~mask_check;

		if (h.nele == 0 && swt != 0)
		{
			// Add the chunks in the empty list
			empty_v.add(active_cnk);
		}
	}

public:

	//! it define that this data-structure is a grid
	typedef int yes_i_am_grid;

	//! base_key for the grid
	typedef grid_key_dx<dim> base_key;

	//! expose the dimansionality as a static const
	static constexpr unsigned int dims = dim;

	//! value that the grid store
	//! The object type the grid is storing
	typedef T value_type;

	//! sub-grid iterator type
	typedef grid_key_sparse_dx_iterator_sub<dim, chunking::size::value> sub_grid_iterator_type;

	/*! \brief Trivial constructor
	 *
	 *
	 */
	inline sgrid_cpu()
	:cache_pnt(0)
	{
		init();
	}

	/*! \brief create a sparse grid from another grid
	 *
	 * \param g the grid to copy
	 *
	 */
	inline sgrid_cpu(const sgrid_cpu & g) THROW
	{
		this->operator=(g);
	}

	/*! \brief create a sparse grid from another grid
	 *
	 * \param g the grid to copy
	 *
	 */
	inline sgrid_cpu(const sgrid_cpu && g) THROW
	{
		this->operator=(g);
	}

	/*! \brief Constructor for sparse grid
	 *
	 * \param sz size in each direction of the sparse grid
	 *
	 */
	sgrid_cpu(const size_t (& sz)[dim])
	:cache_pnt(0),g_sm(sz)
	{
		// calculate the chunks grid

		set_g_shift_from_size(sz,g_sm_shift);

		// fill pos_g

		init();
	}

	/*! \brief Get the background value
	 *
	 * \return background value
	 *
	 */
	T & getBackgroundValue()
	{
		return background;
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

	/*! \brief Insert a full element (with all properties)
	 *
	 * \param v1 where you want to insert the element
	 *
	 */
	inline auto insert_o(const grid_key_dx<dim> & v1, size_t & ele_id) -> decltype(chunks.template get_o(0))
	{
		size_t active_cnk;

		pre_insert(v1,active_cnk,ele_id);

		return chunks.template get_o(active_cnk);
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(chunks.template get<p>(0)[0])>
	inline r_type insert(const grid_key_dx<dim> & v1)
	{
		size_t active_cnk = 0;

		grid_key_dx<dim> kh = v1;
		grid_key_dx<dim> kl;

		// shift the key
		key_shift<dim,chunking>::shift(kh,kl);

		long int lin_id = g_sm_shift.LinId(kh);

		size_t id = 0;
		for (size_t k = 0 ; k < SGRID_CACHE; k++)
		{id += (cache[k] == lin_id)?k+1:0;}

		if (id == 0)
		{
			// we do not have it in cache we check if we have it in the map

			auto fnd = map.find(lin_id);
			if (fnd == map.end())
			{
				// we do not have it in the map create a chunk

				map[lin_id] = chunks.size();
				chunks.add();
				header.add();
				header.last().pos = kh;
				header.last().nele = 0;

				// set the mask to null
				auto & h = header.last().mask;

				for (size_t i = 0 ; i < chunking::size::value / (sizeof(size_t)*8) + (chunking::size::value % (sizeof(size_t)*8) != 0) + 1 ; i++)
				{h[i] = 0;}

				key_shift<dim,chunking>::cpos(header.last().pos);

				active_cnk = chunks.size() - 1;
			}
			else
			{
				// we have it in the map

				active_cnk = fnd->second;
			}

			// Add on cache the chunk
			cache[cache_pnt] = lin_id;
			cached_id[cache_pnt] = active_cnk;
			cache_pnt++;
			cache_pnt = (cache_pnt >= SGRID_CACHE)?0:cache_pnt;
		}
		else
		{
			active_cnk = cached_id[id-1];
			cache_pnt = id;
			cache_pnt = (cache_pnt == SGRID_CACHE)?0:cache_pnt;
		}

		size_t sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

		// the chunk is in cache, solve

		// we notify that we added one element
		auto & h = header.get(active_cnk);

		// we set the mask

		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		h.nele = (h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check)?h.nele:h.nele + 1;
		h.mask[sub_id >> BIT_SHIFT_SIZE_T] |= mask_check;

		return chunks.template get<p>(active_cnk)[sub_id];
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(chunks.template get<p>(0)[0])>
	inline r_type insert(const grid_key_sparse_lin_dx & v1)
	{
		size_t active_cnk = v1.getChunk();
		size_t sub_id = v1.getPos();

		// the chunk is in cache, solve

		// we notify that we added one element
		auto & h = header.get(active_cnk);

		// we set the mask

		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		h.nele = (h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check)?h.nele:h.nele + 1;
		h.mask[sub_id >> BIT_SHIFT_SIZE_T] |= mask_check;

		return chunks.template get<p>(active_cnk)[sub_id];
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p>
	inline auto get(const grid_key_dx<dim> & v1) const -> decltype(openfpm::as_const(chunks.template get<p>(0))[0])
	{
		size_t active_cnk;
		grid_key_dx<dim> kh = v1;
		grid_key_dx<dim> kl;

		// shift the key
		key_shift<dim,chunking>::shift(kh,kl);

		long int lin_id = g_sm_shift.LinId(kh);

		size_t id = 0;
		for (size_t k = 0 ; k < SGRID_CACHE; k++)
		{id += (cache[k] == lin_id)?k+1:0;}

		if (id == 0)
		{
			auto it = map.find(lin_id);

			if (it == map.end())
			{return background.template get<p>();}

			add_on_cache(lin_id,it->second);

			active_cnk = it->second;
		}
		else
		{
			active_cnk = cached_id[id-1];
		}

		size_t sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

		// we check the mask
		auto & h = header.get(active_cnk);

		// We check the mask
		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		if ((h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check) == 0)
		{return background.template get<p>();}

		return chunks.template get<p>(active_cnk)[sub_id];
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p>
	inline auto get(const grid_key_dx<dim> & v1) -> decltype(openfpm::as_const(chunks.template get<p>(0))[0])
	{
		size_t active_cnk;
		grid_key_dx<dim> kh = v1;
		grid_key_dx<dim> kl;

		// shift the key
		key_shift<dim,chunking>::shift(kh,kl);

		long int lin_id = g_sm_shift.LinId(kh);

		size_t id = 0;
		for (size_t k = 0 ; k < SGRID_CACHE; k++)
		{id += (cache[k] == lin_id)?k+1:0;}

		if (id == 0)
		{
			auto it = map.find(lin_id);

			if (it == map.end())
			{return background.template get<p>();}

			add_on_cache(lin_id,it->second);

			active_cnk = it->second;
		}
		else
		{
			active_cnk = cached_id[id-1];
		}

		size_t sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

		// we check the mask
		auto & h = header.get(active_cnk);

		// We check the mask
		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		if ((h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check) == 0)
		{return background.template get<p>();}

		return chunks.template get<p>(active_cnk)[sub_id];
	}

	/*! \brief Check if the point exist
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the true if the point exist
	 *
	 */
	inline bool existPoint(const grid_key_dx<dim> & v1) const
	{
		size_t active_cnk;
		grid_key_dx<dim> kh = v1;
		grid_key_dx<dim> kl;

		// shift the key
		key_shift<dim,chunking>::shift(kh,kl);

		long int lin_id = g_sm_shift.LinId(kh);

		size_t id = 0;
		for (size_t k = 0 ; k < SGRID_CACHE; k++)
		{id += (cache[k] == lin_id)?k+1:0;}

		if (id == 0)
		{
			auto it = map.find(lin_id);

			if (it == map.end())
			{return false;}

			add_on_cache(lin_id,it->second);

			active_cnk = it->second;
		}
		else
		{
			active_cnk = cached_id[id-1];
		}

		size_t sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

		// we check the mask
		auto & h = header.get(active_cnk);

		// We check the mask
		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		if ((h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check) == 0)
		{return false;}

		return true;
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p>
	inline auto get(const grid_key_sparse_lin_dx & v1) -> decltype(chunks.template get<p>(0)[0])
	{
		return chunks.template get<p>(v1.getChunk())[v1.getPos()];
	}

	/*! \brief Return a Domain iterator
	 *
	 * \return return the domain iterator
	 *
	 */
	grid_key_sparse_dx_iterator<dim,chunking::size::value>
	getIterator(size_t opt = 0) const
	{
		return grid_key_sparse_dx_iterator<dim,chunking::size::value>(&header,&pos_chunk);
	}

	/*! \brief Return an iterator over a sub-grid
	 *
	 * \return return an iterator over a sub-grid
	 *
	 */
	grid_key_sparse_dx_iterator_sub<dim,chunking::size::value>
	getIterator(const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop, size_t opt = 0) const
	{
		return grid_key_sparse_dx_iterator_sub<dim,chunking::size::value>(header,pos_chunk,start,stop,sz_cnk);
	}

	/*! \brief Return the internal grid information
	 *
	 * Return the internal grid information
	 *
	 * \return the internal grid
	 *
	 */
	const grid_sm<dim,void> & getGrid() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return g_sm;
	}

	/*! \brief Remove the point
	 *
	 * \param v1 element to remove
	 *
	 */
	void remove(const grid_key_dx<dim> & v1)
	{
		remove_point(v1);
		remove_empty();
	}

	/*! \brief Remove the point but do not flush the remove
	 *
	 * \param v1 element to remove
	 *
	 */
	void remove_no_flush(const grid_key_dx<dim> & v1)
	{
		remove_point(v1);
	}

	/*! \brief Remove the point
	 *
	 * \param v1 element to remove
	 *
	 */
	void flush_remove()
	{
		remove_empty();
	}

	/*! \brief Resize the grid
	 *
	 * The old information is retained on the new grid if the new grid is bigger.
	 * When smaller the data are cropped
	 *
	 * \param sz reference to an array of dimension dim
	 *
	 */
	void resize(const size_t (& sz)[dim])
	{
		bool is_bigger = true;

		// we check if we are resizing bigger, because if is the case we do not have to do
		// much

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (sz[i] < g_sm.size(i))
			{is_bigger = false;}
		}

		g_sm.setDimensions(sz);

		// set g_sm_shift

		set_g_shift_from_size(sz,g_sm_shift);

		clear_cache();

		if (is_bigger == true)
		{

			// if we resize bigger we do not have to do anything in the headers
			// and in the chunks we just have to update g_sm and reset the cache
			// and reconstruct the map. So we reconstruct the map and we just
			// finish

			reconstruct_map();

			return;
		}

		// create a box that is as big as the grid

		Box<dim,size_t> gs_box;

		for (size_t i = 0 ; i < dim ; i++)
		{
			gs_box.setLow(i,0);
			gs_box.setHigh(i,g_sm.size(i));
		}

		// we take a list of all chunks to remove
		openfpm::vector<size_t> rmh;

		// in this case we have to crop data, we go through all the headers

		for (size_t i = 0 ; i < header.size() ; i++)
		{
			Box<dim,size_t> cnk;

			for (size_t j = 0 ; j < dim ; j++)
			{
				cnk.setLow(j,header.get(i).pos.get(j));
				cnk.setHigh(j,sz_cnk[j] + header.get(i).pos.get(j));
			}

			// if the chunk is not fully contained in the new smaller sparse grid
			// we have to crop it
			if (!cnk.isContained(gs_box))
			{
				// We check if the chunks is fully out or only partially in
				// cheking the intersection between the new grid and the box
				// enclosing the chunk as it was before
				Box<dim,size_t> inte;

				if (gs_box.Intersect(cnk,inte))
				{
					// part of the chunk is in, part is out

					// shift P1 to the origin
					// this is the valid box everything out must me reset
					inte -= inte.getP1();

					size_t mask_nele;
					short unsigned int mask_it[chunking::size::value];

					auto & mask = header.get(i).mask;
					auto & n_ele = header.get(i).nele;

					// ok so the box is not fully contained so we must crop data

					fill_mask(mask_it,mask,mask_nele);

					// now we have the mask of all the filled elements

					for (size_t j = 0 ; j < mask_nele ; j++)
					{
						if (!inte.isInside(pos_chunk[mask_it[j]].toPoint()))
						{
							// if is not inside, the point must be deleted

							remove_from_chunk<chunking::size::value>(mask_it[j],n_ele,mask);
						}
					}
				}
				else
				{
					// the chunk is completely out and must be removed completely
					// we add it to the list of the chunks to remove

					rmh.add(i);
				}
			}
		}

		header.remove(rmh,0);
		chunks.remove(rmh,0);

		reconstruct_map();
	}

	/*! \brief Calculate the memory size required to pack n elements
	 *
	 * Calculate the total size required to store n-elements in a vector
	 *
	 * \param n number of elements
	 * \param e unused
	 *
	 * \return the size of the allocation number e
	 *
	 */
	template<int ... prp> static inline size_t packMem(size_t n, size_t e)
	{
		if (sizeof...(prp) == 0)
		{return n * sizeof(typename T::type);}

		typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

		return n * sizeof(prp_object);
	}

	/*! \brief Insert an allocation request
	 *
	 * \tparam prp set of properties to pack
	 *

	 * \param sub sub-grid iterator
	 * \param vector of requests
	 *
	 */
	template<int ... prp> inline
	void packRequest(size_t & req) const
	{
		grid_sm<dim,void> gs_cnk(sz_cnk);

		// For sure we have to pack the number of chunk we want to pack

		req += sizeof(size_t);
		req += dim*sizeof(size_t);

		// Here we have to calculate the number of points to pack

		for (size_t i = 0 ; i < header.size() ; i++)
		{
			auto & h = header.get(i);

			size_t mask_nele;
			short unsigned int mask_it[chunking::size::value];

			fill_mask(mask_it,h.mask,mask_nele);

			for (size_t j = 0 ; j < mask_nele ; j++)
			{
				// If all of the aggregate properties do not have a "pack()" member
				if (has_pack_agg<T,prp...>::result::value == false)
				{
					// here we count how many chunks must be sent

					size_t alloc_ele = this->packMem<prp...>(1,0);
					req += alloc_ele;
				}
				else
				{
					//Call a pack request
					call_aggregatePackRequestChunking<decltype(chunks.template get_o(i)),
																  S,prp ... >
																  ::call_packRequest(chunks.template get_o(i),mask_it[j],req);
				}
			}

			// There are point to send. So we have to save the mask chunk
			req += sizeof(header.get(i).mask);
			// the chunk position
			req += sizeof(header.get(i).pos);
			// and the number of element
			req += sizeof(header.get(i).nele);
		}
	}

	/*! \brief Insert an allocation request
	 *
	 * \tparam prp set of properties to pack
	 *

	 * \param sub sub-grid to pack
	 * \param vector of requests
	 *
	 */
	template<int ... prp> inline
	void packRequest(grid_key_sparse_dx_iterator_sub<dim,chunking::size::value> & sub_it,
					 size_t & req) const
	{
		grid_sm<dim,void> gs_cnk(sz_cnk);

		// For sure we have to pack the number of chunk we want to pack

		req += sizeof(size_t);
		req += dim*sizeof(size_t);

		// Here we have to calculate the number of points to pack

		Box<dim,size_t> section_to_pack;

		for (size_t i = 0; i < dim ; i++)
		{
			section_to_pack.setLow(i,sub_it.getStart().get(i));
			section_to_pack.setHigh(i,sub_it.getStop().get(i));
		}

		for (size_t i = 0 ; i < header.size() ; i++)
		{
			auto & h = header.get(i);

			Box<dim,size_t> bc;

			for (size_t j = 0 ; j < dim ; j++)
			{
				bc.setLow(j,header.get(i).pos.get(j));
				bc.setHigh(j,header.get(i).pos.get(j) + sz_cnk[j] - 1);
			}

			// now we intersect the chunk box with the box

			Box<dim,size_t> inte;
			bool stp = bc.Intersect(section_to_pack,inte);

			if (stp == true)
			{
				// If it is intersect ok we have to check if there are points to pack
				// we shift inte to be relative to the chunk origin

				inte -= header.get(i).pos.toPoint();

				// we iterate all the points

				size_t old_req = req;
				grid_key_dx_iterator_sub<dim> sit(gs_cnk,inte.getKP1(),inte.getKP2());

				while (sit.isNext())
				{
					auto key = sit.get();

					size_t sub_id = gs_cnk.LinId(key);

					size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

					if (h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check)
					{
						// If all of the aggregate properties do not have a "pack()" member
						if (has_pack_agg<T,prp...>::result::value == false)
						{
							// here we count how many chunks must be sent

							size_t alloc_ele = this->packMem<prp...>(1,0);
							req += alloc_ele;
						}
						//If at least one property has "pack()"
						else
						{
							//Call a pack request
							call_aggregatePackRequestChunking<decltype(chunks.template get_o(i)),
																  S,prp ... >
																  ::call_packRequest(chunks.template get_o(i),sub_id,req);
						}
					}

					++sit;
				}

				if (old_req != req)
				{
					// There are point to send. So we have to save the mask chunk
					req += sizeof(header.get(i).mask);
					// the chunk position
					req += sizeof(header.get(i).pos);
					// and the number of element
					req += sizeof(header.get(i).nele);
				}
			}
		}
	}

	/*! \brief Pack the object into the memory given an iterator
	 *
	 * \tparam prp properties to pack
	 *
	 * \param mem preallocated memory where to pack the objects
	 * \param sub_it sub grid iterator ( or the elements in the grid to pack )
	 * \param sts pack statistic
	 *
	 */
	template<int ... prp> void pack(ExtPreAlloc<S> & mem,
									grid_key_sparse_dx_iterator_sub<dims,chunking::size::value> & sub_it,
									Pack_stat & sts)
	{
		grid_sm<dim,void> gs_cnk(sz_cnk);

		// Here we allocate a size_t that indicate the number of chunk we are packing,
		// because we do not know a priory, we will fill it later

		mem.allocate(sizeof(size_t));
		size_t * number_of_chunks = (size_t *)mem.getPointer();

		// Pack the size of the grid

		for (size_t i = 0 ; i < dim ; i++)
		{Packer<size_t,S>::pack(mem,getGrid().size(i),sts);}

		// Here we have to calculate the number of points to pack

		Box<dim,size_t> section_to_pack;

		for (size_t i = 0; i < dim ; i++)
		{
			section_to_pack.setLow(i,sub_it.getStart().get(i));
			section_to_pack.setHigh(i,sub_it.getStop().get(i));
		}

		size_t n_packed_chunk = 0;

		for (size_t i = 0 ; i < header.size() ; i++)
		{
			auto & h = header.get(i);

			Box<dim,size_t> bc;

			for (size_t j = 0 ; j < dim ; j++)
			{
				bc.setLow(j,header.get(i).pos.get(j));
				bc.setHigh(j,header.get(i).pos.get(j) + sz_cnk[j] - 1);
			}

			// now we intersect the chunk box with the box

			Box<dim,size_t> inte;
			bool stp = bc.Intersect(section_to_pack,inte);

			if (stp == true)
			{
				// This flag indicate if something has been packed from this chunk
				bool has_packed = false;

				size_t mask_to_pack[chunking::size::value / (sizeof(size_t)*8) + (chunking::size::value % (sizeof(size_t)*8) != 0) + 1];
				memset(mask_to_pack,0,sizeof(mask_to_pack));
				mem.allocate_nocheck(sizeof(header.get(i).mask) + sizeof(header.get(i).pos) + sizeof(header.get(i).nele));

				// here we get the pointer of the memory in case we have to pack the header
				// and we also shift the memory pointer by an offset equal to the header
				// to pack
				unsigned char * ptr_start = (unsigned char *)mem.getPointer();

				// If it is intersect ok we have to check if there are points to pack
				// we shift inte intp the chunk origin

				inte -= header.get(i).pos.toPoint();

				// we iterate all the points

				grid_key_dx_iterator_sub<dim> sit(gs_cnk,inte.getKP1(),inte.getKP2());

				while (sit.isNext())
				{
					auto key = sit.get();

					size_t sub_id = gs_cnk.LinId(key);

					size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

					if (h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check)
					{
						Packer<decltype(chunks.template get_o(i)),
									S,
									PACKER_ENCAP_OBJECTS_CHUNKING>::template pack<T,prp...>(mem,chunks.template get_o(i),sub_id,sts);

						mask_to_pack[sub_id >> BIT_SHIFT_SIZE_T] |= mask_check;
						has_packed = true;

					}

					++sit;
				}

				if (has_packed == true)
				{
					unsigned char * ptr_final = (unsigned char *)mem.getPointer();
					unsigned char * ptr_final_for = (unsigned char *)mem.getPointerEnd();

					// Ok we packed something so we have to pack the header
					 size_t shift = ptr_final - ptr_start;

					 mem.shift_backward(shift);

					 // The position of the chunks

					 grid_key_dx<dim> pos = header.get(i).pos - sub_it.getStart();

					 Packer<decltype(header.get(i).mask),S>::pack(mem,mask_to_pack,sts);
					 Packer<decltype(header.get(i).pos),S>::pack(mem,pos,sts);
					 Packer<decltype(header.get(i).nele),S>::pack(mem,header.get(i).nele,sts);

					 size_t shift_for = ptr_final_for - (unsigned char *)mem.getPointer();

					 mem.shift_forward(shift_for);

					 n_packed_chunk++;
				}
				else
				{
					// This just reset the last allocation
					mem.shift_backward(0);
				}
			}
		}

		// Now we fill the number of packed chunks
		*number_of_chunks = n_packed_chunk;
	}

	/*! \brief Pack the object into the memory given an iterator
	 *
	 * \tparam prp properties to pack
	 *
	 * \param mem preallocated memory where to pack the objects
	 * \param sub_it sub grid iterator ( or the elements in the grid to pack )
	 * \param sts pack statistic
	 *
	 */
	template<int ... prp> void pack(ExtPreAlloc<S> & mem,
									Pack_stat & sts)
	{
		grid_sm<dim,void> gs_cnk(sz_cnk);

		// Here we allocate a size_t that indicate the number of chunk we are packing,
		// because we do not know a priory, we will fill it later

		Packer<size_t,S>::pack(mem,header.size(),sts);

		for (size_t i = 0 ; i < dim ; i++)
		{Packer<size_t,S>::pack(mem,getGrid().size(i),sts);}

		// Here we have to calculate the number of points to pack

		for (size_t i = 0 ; i < header.size() ; i++)
		{
			auto & h = header.get(i);

			Packer<decltype(header.get(i).mask),S>::pack(mem,h.mask,sts);
			Packer<decltype(header.get(i).pos),S>::pack(mem,h.pos,sts);
			Packer<decltype(header.get(i).nele),S>::pack(mem,h.nele,sts);

			// we iterate all the points

			size_t mask_nele;
			short unsigned int mask_it[chunking::size::value];

			fill_mask(mask_it,h.mask,mask_nele);

			for (size_t j = 0 ; j < mask_nele ; j++)
			{
				Packer<decltype(chunks.template get_o(i)),
								S,
								PACKER_ENCAP_OBJECTS_CHUNKING>::template pack<T,prp...>(mem,chunks.template get_o(i),mask_it[j],sts);
			};
		}
	}

	/*! \brief It does materially nothing
	 *
	 */
	void setMemory()
	{}


	/*! \number of element inserted
	 *
	 * \warning this function is not as fast as the size in other structures
	 *
	 * \return the total number of elements inserted
	 *
	 */
	size_t size() const
	{
		size_t tot = 0;

		for (size_t i = 0 ; i < header.size() ; i++)
		{
			tot += header.get(i).nele;
		}

		return tot;
	}

	/*! \brief Remove all the points in this region
	 *
	 * \param box_src box to kill the points
	 *
	 */
	void remove(Box<dim,size_t> & section_to_delete)
	{
		grid_sm<dim,void> gs_cnk(sz_cnk);

		for (size_t i = 0 ; i < header.size() ; i++)
		{
			auto & h = header.get(i);

			if (i == 9672)
			{
				int debug = 0;
				debug++;
			}

			Box<dim,size_t> bc;

			for (size_t j = 0 ; j < dim ; j++)
			{
				bc.setLow(j,header.get(i).pos.get(j));
				bc.setHigh(j,header.get(i).pos.get(j) + sz_cnk[j] - 1);
			}

			// now we intersect the chunk box with the box

			Box<dim,size_t> inte;
			bool stp = bc.Intersect(section_to_delete,inte);

			if (stp == true)
			{
				// If it is intersect ok we have to check if there are points to pack
				// we shift inte intp the chunk origin

				inte -= header.get(i).pos.toPoint();

				// we iterate all the points

				grid_key_dx_iterator_sub<dim> sit(gs_cnk,inte.getKP1(),inte.getKP2());

				while (sit.isNext())
				{
					auto key = sit.get();

					size_t sub_id = gs_cnk.LinId(key);
					size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

					size_t swt = header.get(i).mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check;

					h.nele = (swt)?h.nele-1:h.nele;
					h.mask[sub_id >> BIT_SHIFT_SIZE_T] &= ~mask_check;

					if (h.nele == 0 && swt != 0)
					{
						// Add the chunks in the empty list
						empty_v.add(i);
					}

					++sit;
				}
			}
		}

		remove_empty();
	}

	void copy_to(const sgrid_cpu<dim,T,S,typename memory_traits_lin<T>::type,chunking> & grid_src,
		         const Box<dim,size_t> & box_src,
			     const Box<dim,size_t> & box_dst)
	{
		auto it = grid_src.getIterator(box_src.getKP1(),box_src.getKP2());

		while (it.isNext())
		{
			auto key_src = it.get();
			grid_key_dx<dim> key_dst = key_src + box_dst.getKP1();
			key_dst -= box_src.getKP1();

			typedef typename std::remove_const<typename std::remove_reference<decltype(grid_src)>::type>::type gcopy;

			copy_sparse_to_sparse<dim,gcopy,gcopy> caps(grid_src,*this,key_src,key_dst);
			boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(caps);

			++it;
		}
	}

	template<template <typename,typename> class op, unsigned int ... prp >
	void copy_to_op(const sgrid_cpu<dim,T,S,typename memory_traits_lin<T>::type,chunking> & grid_src,
		         const Box<dim,size_t> & box_src,
			     const Box<dim,size_t> & box_dst)
	{
		auto it = grid_src.getIterator(box_src.getKP1(),box_src.getKP2());

		while (it.isNext())
		{
			auto key_src = it.get();
			grid_key_dx<dim> key_dst = key_src + box_dst.getKP1();
			key_dst -= box_src.getKP1();

			typedef typename std::remove_const<typename std::remove_reference<decltype(grid_src)>::type>::type gcopy;

			copy_sparse_to_sparse_op<op,dim,gcopy,gcopy,prp ...> caps(grid_src,*this,key_src,key_dst);
			boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(caps);

			++it;
		}
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
				grid_key_sparse_dx_iterator_sub<dims,chunking::size::value> & sub_it,
				Unpack_stat & ps)
	{
		short unsigned int mask_it[chunking::size::value];

		// first we unpack the number of chunks

		size_t n_chunks;

		Unpacker<size_t,S2>::unpack(mem,n_chunks,ps);

		size_t sz[dim];
		for (size_t i = 0 ; i < dim ; i++)
		{Unpacker<size_t,S2>::unpack(mem,sz[i],ps);}

		openfpm::vector<cheader<dim,chunking::size::value>> header_tmp;
		openfpm::vector<aggregate_bfv<chunk_def>> chunks_tmp;

		header_tmp.resize(n_chunks);
		chunks_tmp.resize(n_chunks);

		for (size_t i = 0 ; i < n_chunks ; i++)
		{
			auto & h = header_tmp.get(i);

			Unpacker<decltype(header.get(i).mask),S2>::unpack(mem,h.mask,ps);
			Unpacker<decltype(header.get(i).pos),S2>::unpack(mem,h.pos,ps);
			Unpacker<decltype(header.get(i).nele),S2>::unpack(mem,h.nele,ps);

			// fill the mask_it

			fill_mask(mask_it,h.mask,h.nele);

			// now we unpack the information
			size_t active_cnk;
			size_t ele_id;

			for (size_t k = 0 ; k < h.nele ; k++)
			{
				// construct v1
				grid_key_dx<dim> v1;
				for (size_t i = 0 ; i < dim ; i++)
				{v1.set_d(i,h.pos.get(i) + pos_chunk[mask_it[k]].get(i) + sub_it.getStart().get(i));}

				pre_insert(v1,active_cnk,ele_id);

				Unpacker<decltype(chunks.template get_o(mask_it[k])),
							S2,
							PACKER_ENCAP_OBJECTS_CHUNKING>::template unpack<T,prp...>(mem,chunks.template get_o(active_cnk),ele_id,ps);

			}
		}
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
				Unpack_stat & ps)
	{
		this->clear();

		grid_key_dx<dim> start;
		grid_key_dx<dim> stop;

		// We preunpack some data
		Unpack_stat ps_tmp = ps;

		size_t unused;
		Unpacker<size_t,S2>::unpack(mem,unused,ps_tmp);

		size_t sz[dim];
		for (size_t i = 0 ; i < dim ; i++)
		{Unpacker<size_t,S2>::unpack(mem,sz[i],ps_tmp);}

		g_sm.setDimensions(sz);
		for (size_t i = 0 ; i < dim ; i++)
		{
			start.set_d(i,0);
			stop.set_d(i,getGrid().size(i)-1);
		}

		set_g_shift_from_size(sz,g_sm_shift);

		auto sub_it = this->getIterator(start,stop);

		unpack<prp...>(mem,sub_it,ps);
	}

	/*! \brief unpack the sub-grid object applying an operation
	 *
	 * \tparam op operation
	 * \tparam prp properties to unpack
	 *
	 * \param mem preallocated memory from where to unpack the object
	 * \param sub sub-grid iterator
	 * \param obj object where to unpack
	 *
	 */
	template<template<typename,typename> class op, typename S2, unsigned int ... prp>
	void unpack_with_op(ExtPreAlloc<S2> & mem,
						grid_key_sparse_dx_iterator_sub<dim,chunking::size::value> & sub2,
						Unpack_stat & ps)
	{
		short unsigned int mask_it[chunking::size::value];

		// first we unpack the number of chunks

		size_t n_chunks;

		Unpacker<size_t,S2>::unpack(mem,n_chunks,ps);

		size_t sz[dim];
		for (size_t i = 0 ; i < dim ; i++)
		{Unpacker<size_t,S2>::unpack(mem,sz[i],ps);}

		openfpm::vector<cheader<dim,chunking::size::value>> header_tmp;
		openfpm::vector<aggregate_bfv<chunk_def>> chunks_tmp;

		header_tmp.resize(n_chunks);
		chunks_tmp.resize(n_chunks);

		for (size_t i = 0 ; i < n_chunks ; i++)
		{
			auto & h = header_tmp.get(i);

			Unpacker<decltype(header.get(i).mask),S2>::unpack(mem,h.mask,ps);
			Unpacker<decltype(header.get(i).pos),S2>::unpack(mem,h.pos,ps);
			Unpacker<decltype(header.get(i).nele),S2>::unpack(mem,h.nele,ps);

			// fill the mask_it

			fill_mask(mask_it,h.mask,h.nele);

			// now we unpack the information
			size_t active_cnk;
			size_t ele_id;

			for (size_t k = 0 ; k < h.nele ; k++)
			{
				// construct v1
				grid_key_dx<dim> v1;
				for (size_t i = 0 ; i < dim ; i++)
				{v1.set_d(i,h.pos.get(i) + pos_chunk[mask_it[k]].get(i) + sub2.getStart().get(i));}

				bool exist = pre_insert(v1,active_cnk,ele_id);

				if (exist == false)
				{
					Unpacker<decltype(chunks.template get_o(mask_it[k])),
								S2,
								PACKER_ENCAP_OBJECTS_CHUNKING>::template unpack_op<replace_,prp...>(mem,chunks.template get_o(active_cnk),ele_id,ps);
				}
				else
				{
					Unpacker<decltype(chunks.template get_o(mask_it[k])),
								S2,
								PACKER_ENCAP_OBJECTS_CHUNKING>::template unpack_op<op,prp...>(mem,chunks.template get_o(active_cnk),ele_id,ps);
				}

			}
		}
	}

	/*! \brief This is a meta-function return which type of sub iterator a grid produce
	 *
	 * \return the type of the sub-grid iterator
	 *
	 */
	template <typename stencil = no_stencil>
	static grid_key_sparse_dx_iterator_sub<dim, chunking::size::value>
	type_of_subiterator()
	{
		return  grid_key_sparse_dx_iterator_sub<dim, chunking::size::value>();
	}

	/*! \brief This is a meta-function return which type of sub iterator a grid produce
	 *
	 * \return the type of the sub-grid iterator
	 *
	 */
	static grid_key_sparse_dx_iterator<dim, chunking::size::value>
	type_of_iterator()
	{
		return  grid_key_sparse_dx_iterator<dim, chunking::size::value>();
	}

	/*! \brief Here we convert the linearized sparse key into the grid_key_dx
	 *
	 * \param key_out output key
	 * \param key_in input key
	 *
	 */
	void convert_key(grid_key_dx<dim> & key_out, const grid_key_sparse_lin_dx & key_in) const
	{
		auto & ph = header.get(key_in.getChunk()).pos;
		auto & pos_h = pos_chunk[key_in.getPos()];

		for (size_t i = 0 ; i < dim ; i++)
		{
			key_out.set_d(i,ph.get(i) + pos_h.get(i));
		}
	}

	/*! \brief copy an sparse grid
	 *
	 * \param tmp sparse grdi to copy
	 *
	 */
	sgrid_cpu & operator=(const sgrid_cpu & sg)
	{
		cache_pnt = sg.cache_pnt;
		meta_copy<T>::meta_copy_(sg.background,background);

		for (size_t i = 0 ; i < SGRID_CACHE ; i++)
		{
			cache[i] = sg.cache[i];
			cached_id[i] = sg.cached_id[i];
		}

		//! Map to convert from grid coordinates to chunk
		map = sg.map;
		header = sg.header;
		chunks = sg.chunks;
		g_sm = sg.g_sm;
		g_sm_shift = sg.g_sm_shift;

		for (size_t i = 0 ; i < chunking::size::value ; i++)
		{
			pos_chunk[i] = sg.pos_chunk[i];
		}


		for (size_t i = 0 ; i < dim ; i++)
		{sz_cnk[i] = sg.sz_cnk[i];}

		empty_v = sg.empty_v;

		return *this;
	}

	/*! \brief copy an sparse grid
	 *
	 * \param tmp sparse grdi to copy
	 *
	 */
	sgrid_cpu & operator=(sgrid_cpu && sg)
	{
		cache_pnt = sg.cache_pnt;
		meta_copy<T>::meta_copy_(sg.background,background);

		for (size_t i = 0 ; i < SGRID_CACHE ; i++)
		{
			cache[i] = sg.cache[i];
			cached_id[i] = sg.cached_id[i];
		}

		//! Map to convert from grid coordinates to chunk
		map.swap(sg.map);
		header.swap(sg.header);
		chunks.swap(sg.chunks);
		g_sm = sg.g_sm;
		g_sm_shift = sg.g_sm_shift;

		for (size_t i = 0 ; i < chunking::size::value ; i++)
		{
			pos_chunk[i] = sg.pos_chunk[i];
		}


		for (size_t i = 0 ; i < dim ; i++)
		{sz_cnk[i] = sg.sz_cnk[i];}

		empty_v = sg.empty_v;

		return *this;
	}

	/*! \brief Get the number of inserted points
	 *
	 * \return the number of inserted points
	 *
	 */
	size_t size_inserted()
	{
		return size();
	}

	/*! \brief delete all the points
	 *
	 *
	 */
	void clear()
	{
		header.clear();
		chunks.clear();

		clear_cache();
		reconstruct_map();
	}

	/*! \brief This is an internal function to clear the cache
	 *
	 *
	 *
	 */
	void internal_clear_cache()
	{
		clear_cache();
	}

	/*! \brief write the sparse grid into VTK
	 *
	 * \param out VTK output
	 *
	 */
	template<typename Tw = float> bool write(const std::string & output)
	{
		file_type ft = file_type::BINARY;

		openfpm::vector<Point<dim,Tw>> tmp_pos;
		openfpm::vector<T> tmp_prp;

		// copy position and properties

		auto it = getIterator();

		while(it.isNext())
		{
			auto key = it.getKey();
			auto keyg = it.getKeyF();

			Point<dim,Tw> p;

			for (size_t i = 0 ; i < dim ; i++)
			{p.get(i) = keyg.get(i);}

			tmp_pos.add(p);

			tmp_prp.add();
			copy_prop_to_vector<decltype(chunks.template get_o(key.getChunk())),decltype(tmp_prp.last())>
			cp(chunks.template get_o(key.getChunk()),tmp_prp.last(),key.getPos());

			boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);

			++it;
		}

		// VTKWriter for a set of points
		VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim,Tw>>, openfpm::vector<T>>, VECTOR_POINTS> vtk_writer;
		vtk_writer.add(tmp_pos,tmp_prp,tmp_pos.size());

		openfpm::vector<std::string> prp_names;

		// Write the VTK file
		return vtk_writer.write(output,prp_names,"sparse_grid",ft);
	}

	//Functions to check if the packing object is complex
	static bool pack()
	{
		return false;
	}

	static bool packRequest()
	{
		return false;
	}

	static bool packMem()
	{
		return false;
	}
};


#endif /* OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_HPP_ */
