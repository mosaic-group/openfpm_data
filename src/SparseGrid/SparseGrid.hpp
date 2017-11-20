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
	inline void add_on_cache(size_t lin_id, size_t active_cnk)
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

		size_t sz_cnk[dim];

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
	inline void pre_insert(const grid_key_dx<dim> & v1, size_t & active_cnk, size_t & sub_id)
	{
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

		h.nele = (h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check)?h.nele:h.nele + 1;
		h.mask[sub_id >> BIT_SHIFT_SIZE_T] |= mask_check;
	}

public:

	//! it define that this data-structure is a grid
	typedef int yes_i_am_grid;

	//! expose the dimansionality as a static const
	static constexpr unsigned int dims = dim;

	/*! \brief Trivial constructor
	 *
	 *
	 */
	inline sgrid_cpu()
	:cache_pnt(0)
	{}

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
	template <unsigned int p>
	inline auto get(const grid_key_dx<dim> & v1) const -> decltype(openfpm::as_const(chunks.template get<p>(0))[0])
	{
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

			size_t sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

			add_on_cache(sub_id,it->second);

			// Add on cache the chunk
/*			cache[cache_pnt] = lin_id;
			cached_id[cache_pnt] = it->second;
			cache_pnt++;
			cache_pnt = (cache_pnt >= SGRID_CACHE)?0:cache_pnt;*/

			return chunks.template get<p>(it->second)[sub_id];
		}

		size_t sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

		return chunks.template get<p>(cached_id[id-1])[sub_id];
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

			size_t sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

			add_on_cache(sub_id,it->second);

			// Add on cache the chunk
/*			cache[cache_pnt] = lin_id;
			cached_id[cache_pnt] = it->second;
			cache_pnt++;
			cache_pnt = (cache_pnt >= SGRID_CACHE)?0:cache_pnt;*/

			return chunks.template get<p>(it->second)[sub_id];
		}

		size_t sub_id = sublin<dim,typename chunking::shift_c>::lin(kl);

		return chunks.template get<p>(cached_id[id-1])[sub_id];
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
	grid_key_sparse_dx_iterator<dim,chunking::size::value> getDomainIterator()
	{
		return grid_key_sparse_dx_iterator<dim,chunking::size::value>(header,pos_chunk);
	}

	/*! \brief Return the internal grid information
	 *
	 * Return the internal grid information
	 *
	 * \return the internal grid
	 *
	 */
	const grid_sm<dim,T> & getGrid() const
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

			remove_empty();
		}
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

		size_t sz_cnk[dim];

		copy_sz<dim,typename chunking::type> cpsz(sz_cnk);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,dim> >(cpsz);

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

	/*! \number of element inserted
	 *
	 * \warning this function is not as fast as the size in other structures
	 *
	 * \return the total number of elements inserted
	 *
	 */
	size_t size()
	{
		size_t tot = 0;

		for (size_t i = 0 ; i < header.size() ; i++)
		{
			tot += header.get(i).nele;
		}

		return tot;
	}

	void copy_to(const sgrid_cpu<dim,T,S,chunking> & grid_src,
		         const Box<dim,size_t> & box_src,
			     const Box<dim,size_t> & box_dst)
	{
		size_t sz_cnk[dim];

		copy_sz<dim,typename chunking::type> cpsz(sz_cnk);
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,dim> >(cpsz);

		// we go across all the source chunks and we identify all the chunks in src that must be copied

		for (size_t i = 0 ; i < grid_src.header.size() ; i++)
		{
			Box<dim,size_t> cnk;
			Box<dim,size_t> inte;

			for (size_t j = 0 ; j < dim ; j++)
			{
				cnk.setLow(j,grid_src.header.get(i).pos.get(j));
				cnk.setHigh(j,sz_cnk[j] + grid_src.header.get(i).pos.get(j));
			}

			// now we check if a chunk is inside box_src

			if (cnk.Intersect(box_src,inte))
			{
				size_t mask_nele;
				short int mask_it[chunking::data::size];
				auto & h = grid_src.header.get(i).mask;

				// Fill the mask_it for the source grid

				fill_mask_it(mask_it,h.mask,mask_nele);

				// calculate the index of the element to insert
/*				for (size_t i)

				// now we insert all the elements

				insert_o()

				// It insersect we must move the information*/

/*				if (inte == cnk)
				{
					// we can just copy the chunk

					header.add(grid_src.header.get(i));
					chunks.add(grid_src.chunks.get(i));
				}
				else
				{
					// we only have to copy part of the chunk


				}*/
			}
		}
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

		auto it = getDomainIterator();

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
};


#endif /* OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_HPP_ */
