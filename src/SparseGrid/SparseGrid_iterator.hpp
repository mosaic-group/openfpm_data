/*
 * SparseGrid_iterator.hpp
 *
 *  Created on: Oct 24, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_ITERATOR_HPP_
#define OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_ITERATOR_HPP_

/*! \brief It store the position in space of the sparse grid
 *
 * linearized index
 *
 */
class grid_key_sparse_lin_dx
{
	//! chunk id
	size_t chunk;

	//! linearized id
	size_t lin_id;

public:

	inline grid_key_sparse_lin_dx(size_t chunk, size_t lin_id)
	:chunk(chunk),lin_id(lin_id)
	{

	}

	/*! \brief Return the chunk id
	 *
	 * \return the chunk id
	 *
	 */
	inline size_t getChunk() const
	{
		return chunk;
	}

	/*! \brief Return the linearized position in the chunk
	 *
	 * \return the linearized position in the chunk
	 *
	 */
	inline size_t getPos() const
	{
		return lin_id;
	}
};

/*! \brief This function fill the set of all non zero elements
 *
 *
 */
template<unsigned int n_ele>
inline void fill_mask(short unsigned int (& mask_it)[n_ele],
		       const size_t (& mask)[n_ele / (sizeof(size_t)*8) + (n_ele % (sizeof(size_t)*8) != 0) + 1],
		       size_t & mask_nele)
{
	mask_nele = 0;
	size_t index = 0;

	for (size_t i = 0 ; i < n_ele / (sizeof(size_t) * 8) ; i++)
	{
		size_t mi = mask[i];
		size_t tot_idx = 0;

		do
		{
#if defined(__GNUC__) || defined(__clang__)
			index = __builtin_ffsl(mi);
#elif defined(__INTEL_COMPILER)
			_BitScanForward64(&index,mi)
			index += 1;
#else
			index = 0;
			if (mi != 0)
			{
				while (mi >> index & 0x1 == 0)	{index++;}
				index += 1;
			}
#endif
			if (index != 0)
			{
				mask_it[mask_nele] = (index - 1 + tot_idx) + i*sizeof(size_t)*8;
				mask_nele++;

				mi = (index == 64)?0:mi >> index;
				tot_idx += index;
			}
		}
		while (index != 0);
	}
}

/*! \brief This function fill the set of all non zero elements
 *
 *
 */
template<unsigned int dim, unsigned int n_ele>
inline void fill_mask_box(short unsigned int (& mask_it)[n_ele],
		       const size_t (& mask)[n_ele / (sizeof(size_t)*8) + (n_ele % (sizeof(size_t)*8) != 0) + 1],
		       size_t & mask_nele,
			   Box<dim,size_t> & bx,
			   const grid_key_dx<dim> (& loc_grid)[n_ele])
{
	mask_nele = 0;
	size_t index = 0;

	for (size_t i = 0 ; i < n_ele / (sizeof(size_t) * 8) ; i++)
	{
		size_t mi = mask[i];
		size_t tot_idx = 0;

		do
		{
#if defined(__GNUC__) || defined(__clang__)
			index = __builtin_ffsl(mi);
#elif defined(__INTEL_COMPILER)
			_BitScanForward64(&index,mi)
			index += 1;
#else
			index = 0;
			if (mi != 0)
			{
				while (mi >> index & 0x1 == 0)	{index++;}
				index += 1;
			}
#endif
			if (index != 0)
			{
				size_t id = (index - 1 + tot_idx) + i*sizeof(size_t)*8;

				bool is_inside = true;
				// we check the point is inside inte
				for (size_t i = 0 ; i < dim ; i++)
				{
					if (loc_grid[id].get(i) < (long int)bx.getLow(i) ||
						loc_grid[id].get(i) > (long int)bx.getHigh(i))
					{
						is_inside = false;

						break;
					}
				}

				if (is_inside == true)
				{
					mask_it[mask_nele] = id;
					mask_nele++;
				}

				mi = (index == 64)?0:mi >> index;
				tot_idx += index;
			}
		}
		while (index != 0);
	}
}

/*! \brief This structure contain the information of a chunk
 *
 * \tparam dim dimensionality of the chunk
 * \tparam n_ele number of elements in the chunk
 *
 */
template<unsigned int dim, unsigned int n_ele>
struct cheader
{
	//! where is located the chunk
	grid_key_dx<dim> pos;

	//! how many elements in the chunk are set
	size_t nele;

	//! which elements in the chunks are set
	size_t mask[n_ele / (sizeof(size_t)*8) + (n_ele % (sizeof(size_t)*8) != 0) + 1];
};


/*! \brief Grid key sparse iterator on a sub-part of the domain
 *
 *
 */
template<unsigned dim, unsigned int n_ele>
class grid_key_sparse_dx_iterator_sub
{
	//! It store the information of each chunk
	const openfpm::vector<cheader<dim,n_ele>> * header;

	//! linearized id to position in coordinates conversion
	const grid_key_dx<dim> (* lin_id_pos)[n_ele];

	//! point to the actual chunk
	size_t chunk_id;

	//! Number of points in mask_it
	size_t mask_nele;

	//! Actual point in mask it
	size_t mask_it_pnt;

	//! Starting point
	grid_key_dx<dim> start;

	//! Stop point
	grid_key_dx<dim> stop;

	//! Sub-grid box
	Box<dim,size_t> bx;

	//! Chunk size
	size_t sz_cnk[dim];

	//! set of index in the chunk on which we have to iterate
	short unsigned int mask_it[n_ele];

	/*! \brief Everytime we move to a new chunk we calculate on which indexes we have to iterate
	 *
	 *
	 */
	void SelectValidAndFill_mask_it()
	{
		mask_it_pnt = 0;
		mask_nele = 0;

		while (mask_nele == 0 && chunk_id < header->size())
		{
			auto & mask = header->get(chunk_id).mask;

			Box<dim,size_t> cnk_box;

			for (size_t i = 0 ; i < dim ; i++)
			{
				cnk_box.setLow(i,header->get(chunk_id).pos.get(i));
				cnk_box.setHigh(i,header->get(chunk_id).pos.get(i) + sz_cnk[i] - 1);
			}

			Box<dim,size_t> inte;

			if (bx.Intersect(cnk_box,inte) == true)
			{
				inte -= header->get(chunk_id).pos.toPoint();
				fill_mask_box<dim,n_ele>(mask_it,mask,mask_nele,inte,*lin_id_pos);
			}
			else
			{mask_nele = 0;}

			chunk_id = (mask_nele == 0)?chunk_id + 1:chunk_id;

		}
	}

public:

	/*! \brief Default constructor
	 *
	 * \warning extremely unsafe
	 * If you use this constructor before use the iterator you should call reinitialize first
	 *
	 */
	grid_key_sparse_dx_iterator_sub()	{};

	grid_key_sparse_dx_iterator_sub(const openfpm::vector<cheader<dim,n_ele>> & header,
								const grid_key_dx<dim> (& lin_id_pos)[n_ele],
								const grid_key_dx<dim> & start,
								const grid_key_dx<dim> & stop,
								const size_t (& sz_cnk)[dim])
	:header(&header),lin_id_pos(&lin_id_pos),chunk_id(0),
	 mask_nele(0),mask_it_pnt(0),start(start),stop(stop)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			this->sz_cnk[i] = sz_cnk[i];

			bx.setLow(i,start.get(i));
			bx.setHigh(i,stop.get(i));
		}

		SelectValidAndFill_mask_it();
	}

	/*! \brief Reinitialize the iterator
	 *
	 * it re-initialize the iterator with the passed grid_key_dx_iterator_sub
	 * the actual position of the grid_key_dx_iterator_sub is ignored
	 *
	 * \param g_s_it grid_key_dx_iterator_sub
	 *
	 */
	inline void reinitialize(const grid_key_sparse_dx_iterator_sub<dim,n_ele> & g_s_it)
	{
		header = g_s_it.header;
		lin_id_pos = g_s_it.lin_id_pos;
		chunk_id = g_s_it.chunk_id;
		mask_nele = g_s_it.mask_nele;
		mask_it_pnt = g_s_it.mask_it_pnt;
		start = g_s_it.start;
		stop = g_s_it.stop;
		bx = g_s_it.bx;
		for (size_t i = 0 ; i < dim ; i++)
		{sz_cnk[i] = g_s_it.sz_cnk[i];}

		memcpy(mask_it,g_s_it.mask_it,sizeof(short unsigned int)*n_ele);
	}

	inline grid_key_sparse_dx_iterator_sub<dim,n_ele> & operator++()
	{
		mask_it_pnt++;

		if (mask_it_pnt < mask_nele)
		{
			return *this;
		}

		chunk_id++;
		mask_it_pnt = 0;

		if (chunk_id < header->size())
		{
			SelectValidAndFill_mask_it();
		}

		return *this;
	}


	/*! \brief Return the position of the actual point in coordinates
	 *
	 * \return the position
	 *
	 */
	inline grid_key_dx<dim> get() const
	{
		grid_key_dx<dim> k_pos;

		size_t lin_id = mask_it[mask_it_pnt];

		for (size_t i = 0 ; i < dim ; i++)
		{
			k_pos.set_d(i,(*lin_id_pos)[lin_id].get(i) + header->get(chunk_id).pos.get(i));
		}

		return k_pos;
	}

	template<typename> struct Debug;

	/*! \brief Return the actual point linearized
	 *
	 * The disadvantage is that in general you cannot use move in this type of key
	 *
	 * \return the actual point
	 *
	 */
	inline grid_key_sparse_lin_dx getKeyF()
	{
		return grid_key_sparse_lin_dx(chunk_id,mask_it[mask_it_pnt]);
	}

	/*! \brief Return true if there is a next grid point
	 *
	 * \return true if there is the next grid point
	 *
	 */
	bool isNext()
	{
		return chunk_id < header->size();
	}

	/*! \brief Return the starting point for the iteration
	 *
	 * \return the starting point
	 *
	 */
	const grid_key_dx<dim> & getStart() const
	{
		return start;
	}

	/*! \brief Return the stop point for the iteration
	 *
	 * \return the stop point
	 *
	 */
	const grid_key_dx<dim> & getStop() const
	{
		return stop;
	}

	/*! \brief Return the private member header
	 *
	 * \return header
	 *
	 */
	const openfpm::vector<cheader<dim,n_ele>> * private_get_header() const
	{return header;}

	/*! \brief Return the private member lin_id_pos
	 *
	 * \return lin_id_pos
	 *
	 */
	const grid_key_dx<dim> (* private_get_lin_id_pos() const)[n_ele]
	{return lin_id_pos;}

	/*! \brief Return the private member chunk_id
	 *
	 * \return chunk_id
	 *
	 */
	size_t private_get_chunk_id() const
	{return chunk_id;}

	/*! \brief Return the private member mask_nele
	 *
	 * \return mask_nele
	 *
	 */
	size_t private_get_mask_nele() const
	{return mask_nele;}


	/*! \brief Return the private member mask_it_pnt
	 *
	 * \return mask_it_pnt
	 *
	 */
	size_t private_get_mask_it_pnt() const
	{return mask_it_pnt;}

	/*! \brief Return the private member mask_it
	 *
	 * \return mask_it
	 *
	 */
	const short unsigned int (& private_get_mask_it() const) [n_ele]
	{
		return mask_it;
	}
};

/*! \brief Grid key sparse iterator
 *
 *
 */
template<unsigned dim, unsigned int n_ele>
class grid_key_sparse_dx_iterator
{
	//! It store the information of each chunk
	const openfpm::vector<cheader<dim,n_ele>> * header;

	//! linearized id to position in coordinates conversion
	const grid_key_dx<dim> (* lin_id_pos)[n_ele];

	//! point to the actual chunk
	size_t chunk_id;

	//! Number of points in mask_it
	size_t mask_nele;

	//! Actual point in mask it
	size_t mask_it_pnt;

	//! set of index in the chunk on which we have to iterate
	short unsigned int mask_it[n_ele];

	/*! \brief Everytime we move to a new chunk we calculate on which indexes we have to iterate
	 *
	 *
	 */
	void SelectValidAndFill_mask_it()
	{
		mask_nele = 0;
		mask_it_pnt = 0;

		while (mask_nele == 0 && chunk_id < header->size())
		{
			auto & mask = header->get(chunk_id).mask;

			fill_mask<n_ele>(mask_it,mask,mask_nele);

			chunk_id = (mask_nele == 0)?chunk_id + 1:chunk_id;

		}
	}

public:

	/*! \brief Default constructor
	 *
	 * \warning extremely unsafe
	 * If you use this constructor before use the iterator you should call reinitialize first
	 *
	 */
	grid_key_sparse_dx_iterator()	{};

	grid_key_sparse_dx_iterator(const openfpm::vector<cheader<dim,n_ele>> * header,
								const grid_key_dx<dim> (* lin_id_pos)[n_ele])
	:header(header),lin_id_pos(lin_id_pos),chunk_id(0),mask_nele(0),mask_it_pnt(0)
	{
		SelectValidAndFill_mask_it();
	}

	inline grid_key_sparse_dx_iterator<dim,n_ele> & operator++()
	{
		mask_it_pnt++;

		if (mask_it_pnt < mask_nele)
		{
			return *this;
		}

		chunk_id++;
		mask_it_pnt = 0;

		if (chunk_id < header->size())
		{
			SelectValidAndFill_mask_it();
		}

		return *this;
	}

	/*! \brief Reinitialize the iterator
	 *
	 * it re-initialize the iterator with the passed grid_key_dx_iterator_sub
	 * the actual position of the grid_key_dx_iterator_sub is ignored
	 *
	 * \param g_s_it grid_key_dx_iterator
	 *
	 */
	inline void reinitialize(const grid_key_sparse_dx_iterator<dim,n_ele> & g_s_it)
	{
		header = g_s_it.header;
		lin_id_pos = g_s_it.lin_id_pos;
		chunk_id = g_s_it.chunk_id;
		mask_nele = g_s_it.mask_nele;
		mask_it_pnt = g_s_it.mask_it_pnt;

		memcpy(mask_it,g_s_it.mask_it,sizeof(short unsigned int)*n_ele);
	}

	/*! \brief Reinitialize the iterator
	 *
	 * it re-initialize the iterator with the passed grid_key_dx_iterator_sub
	 * the actual position of the grid_key_dx_iterator_sub is ignored
	 *
	 * \param g_s_it grid_key_dx_iterator
	 *
	 */
	inline void reinitialize(const grid_key_sparse_dx_iterator_sub<dim,n_ele> & g_s_it)
	{
		header = g_s_it.private_get_header();
		lin_id_pos = g_s_it.private_get_lin_id_pos();
		chunk_id = 0;

		SelectValidAndFill_mask_it();
	}

	/*! \brief Return the actual point
	 *
	 * \warning the key that this function return cannot be used with function like move
	 *
	 * \return the actual point
	 *
	 */
	inline grid_key_sparse_lin_dx getKeyF()
	{
		return grid_key_sparse_lin_dx(chunk_id,mask_it[mask_it_pnt]);
	}

	/*! \brief Return the position of the actual point in coordinates
	 *
	 * \return the position
	 *
	 */
	inline grid_key_dx<dim> get()
	{
		grid_key_dx<dim> k_pos;

		size_t lin_id = mask_it[mask_it_pnt];

		for (size_t i = 0 ; i < dim ; i++)
		{
			k_pos.set_d(i,(*lin_id_pos)[lin_id].get(i) + header->get(chunk_id).pos.get(i));
		}

		return k_pos;
	}

	/*! \brief Return true if there is a next grid point
	 *
	 * \return true if there is the next grid point
	 *
	 */
	bool isNext()
	{
		return chunk_id < header->size();
	}
};

#endif /* OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_ITERATOR_HPP_ */
