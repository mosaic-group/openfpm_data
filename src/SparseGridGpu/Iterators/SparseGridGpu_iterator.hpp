/*
 * SparseGridGpu_iterator.hpp
 *
 *  Created on: Sep 4, 2019
 *      Author: i-bird
 */

#ifndef SPARSEGRIDGPU_ITERATOR_HPP_
#define SPARSEGRIDGPU_ITERATOR_HPP_

/*! \brief Element index contain a data chunk index and a point index
 *
 * \tparam SparseGridGpu type
 *
 */
template<typename SparseGridGpu_type>
class sparse_grid_gpu_index
{
	//! chunk position id
	unsigned int cnk_pos_id;

	//! data id
	unsigned int data_id;

	//! SparseGridGpu used to add functionalities
	const SparseGridGpu_type & sparseGrid;

public:

	/*! \brief Constructor from SparseGridGpu
	 *
	 *
	 */
	inline sparse_grid_gpu_index(const SparseGridGpu_type & sparseGrid, int cnk_pos_id, int data_id)
	:sparseGrid(sparseGrid),cnk_pos_id(cnk_pos_id),data_id(data_id)
	{}

	/*! \brief Convert to a point this index
	 *
	 * \see toPointS
	 *
	 * \return a point unsigned long int
	 *
	 */
	inline Point<SparseGridGpu_type::dims,size_t> toPoint() const
	{
		auto indexCnk = sparseGrid.getIndex(cnk_pos_id);

		auto coord = sparseGrid.getCoord(indexCnk*sparseGrid.getBlockSize() + data_id);

		Point<SparseGridGpu_type::dims,size_t> p;

		for (size_t i = 0; i < SparseGridGpu_type::dims ; i++)
		{
			p.get(i) = coord.get(i);
		}

		return p;
	}

	/*! \brief Get chunk position id
	 *
	 * Return the position of the chunk in the chunks array \see SparseGridGpu \see private_get_data_array
	 *
	 * \return Get chunk position id
	 *
	 */
	int get_cnk_pos_id() const
	{
		return cnk_pos_id;
	}

	/*! \brief Get chunk local index (the returned index < getblockSize())
	 *
	 * \return Get chunk position id
	 *
	 */
	int get_data_id() const
	{
		return data_id;
	}

	/*! \brief Get chunk position id
	 *
	 * Return the position of the chunk in the chunks array \see SparseGridGpu \see private_get_data_array
	 *
	 * \return Get chunk position id
	 *
	 */
	unsigned int & get_cnk_pos_id_ref()
	{
		return cnk_pos_id;
	}

	/*! \brief Get chunk local index (the returned index < getblockSize())
	 *
	 * \return Get chunk position id
	 *
	 */
	unsigned int & get_data_id_ref()
	{
		return data_id;
	}

	/*! \brief Set chunk position id
	 *
	 *
	 * \param cnk_pos_id chunk position id
	 *
	 */
	void set_cnk_pos_id(int cnk_pos_id)
	{
		this->cnk_pos_id = cnk_pos_id;
	}

	/*! \brief Set chunk local index (the returned index < getblockSize())
	 *
	 * \param data_id local id
	 *
	 */
	void set_data_id(int data_id)
	{
		this->data_id = data_id;
	}

	/*! \brief return toPoint() + p
	 *
	 * \param p the point p
	 *
	 * \return toPoint() + p
	 *
	 */
	inline grid_key_dx<SparseGridGpu_type::dims> operator+(const Point<SparseGridGpu_type::dims,size_t> & p)
	{
		grid_key_dx<SparseGridGpu_type::dims> ret;
		Point<SparseGridGpu_type::dims,size_t> key = toPoint();

		for (int i = 0 ; i < SparseGridGpu_type::dims ; i++)
		{
			ret.set_d(i,key.get(i) + p.get(i));
		}

		return ret;
	}

	/*! \brief set the position of this key in direction i to n
	 *
	 * \param i direction
	 * \param n position
	 *
	 */
	inline void set_d(size_t i,int n)
	{
		auto indexCnk = sparseGrid.getIndex(cnk_pos_id);

		auto coord = sparseGrid.getCoord(indexCnk*sparseGrid.getBlockSize() + data_id);

		coord.set_d(i,n);

		auto key = sparseGrid.get_sparse(coord);

		cnk_pos_id = key.cnk_pos_id;
		data_id = key.data_id;
	}

	/* \brief Return the coordinate in direction i
	 *
	 * \return the coordinates in direction i
	 *
	 */
	inline unsigned int get(unsigned int i) const
	{
		auto indexCnk = sparseGrid.getIndex(cnk_pos_id);

		return sparseGrid.getCoord(indexCnk*sparseGrid.getBlockSize() + data_id).get(i);
	}
};

template<unsigned int dim, typename SparseGridType>
class SparseGridGpu_iterator
{
	//! actual chunk
	unsigned int chunk;

	//! actual point inside the chunk
	unsigned int pnt;

	//! original SparseGrid
	const SparseGridType * sparseGrid;

	//! Select the first valid point chunk
	void SelectValid()
	{
		while (pnt < sparseGrid->getBlockSize() && sparseGrid->getMask(chunk,pnt) == 0)
		{
			pnt++;
		}

		while (pnt == sparseGrid->getBlockSize() && chunk < sparseGrid->numBlocks())
		{
			chunk++;
			pnt = 0;
			while (pnt < sparseGrid->getBlockSize() && sparseGrid->getMask(chunk,pnt) == 0)
			{
				pnt++;
			}
		}
	}

public:


	inline SparseGridGpu_iterator()
	:chunk(0),
	 pnt(0),
	 sparseGrid(NULL)
	{
	}


	/*! \brief Constructor
	 *
	 * \param ids vector of ids
	 * \para data vector of chunk data
	 *
	 */
	inline SparseGridGpu_iterator(const SparseGridType & sparseGrid)
	:chunk(0),
	 pnt(0),
	 sparseGrid(&sparseGrid)
	{
		SelectValid();
	}

	/*! \brief Check if there is the next element
	 *
	 * Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */
	bool isNext() const
	{
		return chunk < sparseGrid->numBlocks();
	}

	/*! \brief Get the next element
	 *
	 * Get the next element
	 *
	 * \return the next grid_key
	 *
	 */
	inline SparseGridGpu_iterator<dim,SparseGridType> & operator++()
	{
		++pnt;

		if (pnt >= sparseGrid->getBlockSize())
		{
			++chunk;
			pnt = 0;
		}

		SelectValid();

		return *this;
	}

	/*! \brief return the actual point
	 *
	 * \return the index of the actual point
	 *
	 */
	inline sparse_grid_gpu_index<SparseGridType> get() const
	{
		sparse_grid_gpu_index<SparseGridType> spgi(*sparseGrid,chunk,pnt);

		return spgi;
	}

	/*! \brief Reinitialize the iterator
	 *
	 * \param it_sub subiterator
	 *
	 */
	inline void reinitialize(const SparseGridGpu_iterator<dim,SparseGridType> & it)
	{
		this->operator=(it);
	}

	SparseGridGpu_iterator<dim,SparseGridType> & operator=(const SparseGridGpu_iterator<dim,SparseGridType> & it)
	{
		chunk = it.chunk;
		sparseGrid = it.sparseGrid;
		pnt = it.pnt;

		return *this;
	}
};


#endif /* SPARSEGRIDGPU_ITERATOR_HPP_ */
