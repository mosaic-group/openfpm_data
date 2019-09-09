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
	int cnk_pos_id;

	//! data id
	int data_id;

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
		auto indexCnk = sparseGrid.private_get_index_array().template get<0>(cnk_pos_id);

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
};

template<unsigned int dim, typename SparseGridType>
class SparseGridGpu_iterator
{
	//! actual chunk
	unsigned int chunk;

	//! actual point inside the chunk
	unsigned int pnt;

	grid_key_dx<dim> a_cnk;

	//! original SparseGrid
	const SparseGridType & sparseGrid;

	//! array type for the indexes
	typedef typename std::remove_reference<decltype(sparseGrid.private_get_index_array())>::type index_array_type;

	//! array type for the data
	typedef typename std::remove_reference<decltype(sparseGrid.private_get_data_array())>::type data_array_type;

	//! vector of the chunk indexes
	const decltype(sparseGrid.private_get_index_array()) & ids;

	//! vector containing each chunks datas
	const decltype(sparseGrid.private_get_data_array()) & data;

	//Get the chunk type
	typedef typename boost::mpl::at<typename data_array_type::value_type::type,boost::mpl::int_<0>>::type chunk_type;

	//Get the chunk type
	typedef boost::mpl::int_<boost::mpl::size<typename data_array_type::value_type::type>::type::value-1> pMask;

	//! Select the first valid point chunk
	void SelectValid()
	{
		while (pnt < chunk_type::size && data.template get<pMask::value>(chunk)[pnt] == 0)
		{
			pnt++;
		}

		while (pnt == chunk_type::size && chunk < ids.size())
		{
			chunk++;
			pnt = 0;
			while (pnt < chunk_type::size && data.template get<pMask::value>(chunk)[pnt] == 0)
			{
				pnt++;
			}
		}
	}

public:

	/*! \brief Constructor
	 *
	 * \param ids vector of ids
	 * \para data vector of chunk data
	 *
	 */
	inline SparseGridGpu_iterator(const SparseGridType & sparseGrid)
	:chunk(0),
	 pnt(0),
	 sparseGrid(sparseGrid),
	 ids(sparseGrid.private_get_index_array()),
	 data(sparseGrid.private_get_data_array())
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
		return chunk < ids.size();
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

		if (pnt >= chunk_type::size)
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
	inline sparse_grid_gpu_index<SparseGridType> get()
	{
		sparse_grid_gpu_index<SparseGridType> spgi(sparseGrid,chunk,pnt);

		return spgi;
	}
};


#endif /* SPARSEGRIDGPU_ITERATOR_HPP_ */
