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
		       size_t (& mask)[n_ele / (sizeof(size_t)*8) + (n_ele % (sizeof(size_t)*8) != 0) + 1],
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

/*! \brief Grid key sparse iterator
 *
 *
 */
template<unsigned dim, unsigned int n_ele>
class grid_key_sparse_dx_iterator
{
	//! It store the information of each chunk
	openfpm::vector<cheader<dim,n_ele>> & header;

	//! linearized id to position in coordinates conversion
	grid_key_dx<dim> (& lin_id_pos)[n_ele];

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

		while (mask_nele == 0 && chunk_id < header.size())
		{
			auto & mask = header.get(chunk_id).mask;

			fill_mask<n_ele>(mask_it,mask,mask_nele);

			chunk_id = (mask_nele == 0)?chunk_id + 1:chunk_id;

		}
	}

public:

	grid_key_sparse_dx_iterator(openfpm::vector<cheader<dim,n_ele>> & header,
								grid_key_dx<dim> (& lin_id_pos)[n_ele])
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

		if (chunk_id < header.size())
		{
			SelectValidAndFill_mask_it();
		}

		return *this;
	}

	/*! \brief Return the actual point
	 *
	 * \return the actual point
	 *
	 */
	inline grid_key_sparse_lin_dx getKey()
	{
		return grid_key_sparse_lin_dx(chunk_id,mask_it[mask_it_pnt]);
	}

	/*! \brief Return the position of the actual point in coordinates
	 *
	 * \return the position
	 *
	 */
	inline grid_key_dx<dim> getKeyF()
	{
		grid_key_dx<dim> k_pos;

		size_t lin_id = mask_it[mask_it_pnt];

		for (size_t i = 0 ; i < dim ; i++)
		{
			k_pos.set_d(i,lin_id_pos[lin_id].get(i) + header.get(chunk_id).pos.get(i));
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
		return chunk_id < header.size();
	}
};

#endif /* OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_ITERATOR_HPP_ */
