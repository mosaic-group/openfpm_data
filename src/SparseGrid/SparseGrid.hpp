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

			empty_v.clear();

			// cache must be cleared

			cache_pnt = 0;
			for (size_t i = 0 ; i < SGRID_CACHE ; i++)
			{cache[i] = -1;}
		}
	}

public:

	/*! \brief Constructor for sparse grid
	 *
	 * \param sz size in each direction of the sparse grid
	 *
	 */
	sgrid_cpu(const size_t (& sz)[dim])
	:cache_pnt(0),g_sm(sz)
	{
		for (size_t i = 0 ; i < SGRID_CACHE ; i++)
		{cache[i] = -1;}

		// calculate the chunks grid

		grid_key_dx<dim> cs;
		grid_key_dx<dim> unused;

		for (size_t i = 0 ; i < dim ; i++)
		{cs.set_d(i,sz[i]);}

		key_shift<dim,chunking>::shift(cs,unused);

		size_t sz_i[dim];

		for (size_t i = 0 ; i < dim ; i++)
		{sz_i[i] = cs.get(i) + 1;}

		g_sm_shift.setDimensions(sz_i);

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

	/*! \brief Get the background value
	 *
	 * \return background value
	 *
	 */
	T & getBackgroundValue()
	{
		return background;
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

			// Add on cache the chunk
			cache[cache_pnt] = lin_id;
			cached_id[cache_pnt] = it->second;
			cache_pnt++;
			cache_pnt = (cache_pnt >= SGRID_CACHE)?0:cache_pnt;

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
	template <unsigned int p, typename r_type=decltype(chunks.template get<p>(0)[0])>
	inline r_type get(const grid_key_sparse_lin_dx & v1)
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
