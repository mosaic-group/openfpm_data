/*
 * SparseGridGpu_iterator_sub.hpp
 *
 *  Created on: Jun 24, 2019
 *      Author: i-bird
 */

#ifndef SPARSEGRIDGPU_ITERATOR_SUB_HPP_
#define SPARSEGRIDGPU_ITERATOR_SUB_HPP_

#include "SparseGridGpu_iterator.hpp"

template<unsigned int dim, typename SparseGridType>
class SparseGridGpu_iterator_sub
{
	//! actual chunk
	unsigned int chunk;

	//! original SparseGrid
	const SparseGridType * sparseGrid;

	//! chunk coordinates
	grid_key_dx<dim,int> chunk_coord;

	//! array type for the indexes
	typedef typename std::remove_reference<decltype(sparseGrid->private_get_index_array())>::type index_array_type;

	//! array type for the data
	typedef typename std::remove_reference<decltype(sparseGrid->private_get_data_array())>::type data_array_type;

	//! vector of the chunk indexes
	const typename std::remove_reference<decltype(sparseGrid->private_get_index_array())>::type * ids;

	//! vector containing each chunks datas
	const typename std::remove_reference<decltype(sparseGrid->private_get_data_array())>::type * data;

	//Get the chunk type
	typedef typename boost::mpl::at<typename data_array_type::value_type::type,boost::mpl::int_<0>>::type chunk_type;

	//Get the chunk type
	typedef boost::mpl::int_<boost::mpl::size<typename data_array_type::value_type::type>::type::value-1> pMask;

	// subset in grid coordinates
	Box<dim,int> sub_set;

	//! start point
	mutable grid_key_dx<dim> start;

	//! stop point
	mutable grid_key_dx<dim> stop;

	// subset in local chunk coordinates
	Box<dim,int> res;

	// in chunk iterator
	grid_key_dx_iterator_sub<dim> in_chunk_it;

	// chunk size
	grid_sm<dim,void> chunk_sz;

	/*! \brief initialize the chunk interator
	 *
	 *
	 */
	void initialize_chunk_it()
	{
		if (chunk >= ids->size())
		{return;}

		// compute if the chunk intersect the start - stop box
		chunk_coord = sparseGrid->getCoord(ids->template get<0>(chunk)*sparseGrid->getBlockSize());

		Box<dim,int> box;

		for (int i = 0 ; i < dim ; i++)
		{
			box.setLow(i,chunk_coord.get(i));
			box.setHigh(i,chunk_coord.get(i) + sparseGrid->getBlockEdgeSize() - 1);
		}


		if (sub_set.Intersect(box,res) == true)
		{
			// remove the offset
			for (int i = 0 ; i < dim ; i++)
			{
				res.setLow(i,res.getLow(i) - chunk_coord.get(i));
				res.setHigh(i,res.getHigh(i) - chunk_coord.get(i));
			}

			in_chunk_it.reinitialize(grid_key_dx_iterator_sub<dim>(chunk_sz,res.getKP1(),res.getKP2()));

			while (in_chunk_it.isNext() == true && (data->template get<pMask::value>(chunk)[chunk_sz.LinId(in_chunk_it.get())]&1ul) == 0)
			{
				++in_chunk_it;
			}
		}
	}

	//! Select the first valid point chunk
	void SelectValid()
	{
		while (in_chunk_it.isNext() == true && (data->template get<pMask::value>(chunk)[chunk_sz.LinId(in_chunk_it.get())]&1ul) == 0)
		{
			++in_chunk_it;
		}

		while (in_chunk_it.isNext() == false && chunk < ids->size())
		{
			chunk++;

			initialize_chunk_it();
		}
	}

	/*! \brief Initialize chunk_sz member
	 *
	 *
	 */
	void initialize_chunk_sz()
	{
		size_t sz[dim];

		for (int i = 0 ; i < dim ; i++)
		{sz[i] = sparseGrid->getBlockEdgeSize();}

		chunk_sz.setDimensions(sz);
	}

public:

	/*! \brief Default constructor
	 *
	 */
	inline SparseGridGpu_iterator_sub()
	:chunk(0),
	 sparseGrid(NULL),
	 ids(NULL),
	 data(NULL)
	{
		initialize_chunk_sz();
	}

	/*! \brief Constructor
	 *
	 * \param sparseGrid original sparse grid
	 * \param start starting point
	 * \param stop stop point
	 *
	 */
	inline SparseGridGpu_iterator_sub(const SparseGridType & sparseGrid,const grid_key_dx<dim> & start,const  grid_key_dx<dim> & stop)
	:chunk(0),
	 sparseGrid(&sparseGrid),
	 ids(&sparseGrid.private_get_index_array()),
	 data(&sparseGrid.private_get_data_array())
	{
		for (int i = 0; i < dim ; i++)
		{
			sub_set.setLow(i,start.get(i));
			sub_set.setHigh(i,stop.get(i));
		}

		initialize_chunk_sz();
		in_chunk_it.invalidate();
		initialize_chunk_it();

		SelectValid();
	}

	/*! \brief Return the starting point
	 *
	 * \return the start point
	 *
	 */
	grid_key_dx<dim> & getStart() const
	{
		start = sub_set.getKP1();

		return start;
	}

	/*! \brief Return the stop point
	 *
	 * \return the stop point
	 *
	 */
	grid_key_dx<dim> & getStop() const
	{
		stop = sub_set.getKP2();

		return stop;
	}

	/*! \brief Reinitialize the iterator
	 *
	 * \param it_sub subiterator
	 *
	 */
	inline void reinitialize(const SparseGridGpu_iterator_sub<dim,SparseGridType> & it_sub)
	{
		this->operator=(it_sub);
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
		return chunk < ids->size();
	}

	/*! \brief Get the next element
	 *
	 * Get the next element
	 *
	 * \return the next grid_key
	 *
	 */
	inline SparseGridGpu_iterator_sub<dim,SparseGridType> & operator++()
	{
		++in_chunk_it;

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
		sparse_grid_gpu_index<SparseGridType> spgi(*sparseGrid,chunk,chunk_sz.LinId(in_chunk_it.get()));

		return spgi;
	}

	SparseGridGpu_iterator_sub<dim,SparseGridType> & operator=(const SparseGridGpu_iterator_sub<dim,SparseGridType> & it_sub)
	{
		chunk = it_sub.chunk;
		sparseGrid = it_sub.sparseGrid;
		chunk_coord = it_sub.chunk_coord;
		ids = it_sub.ids;
		data = it_sub.data;
		sub_set = it_sub.sub_set;
		res = it_sub.res;
		in_chunk_it = it_sub.in_chunk_it;
		chunk_sz = it_sub.chunk_sz;

		return *this;
	}
};


#endif /* SPARSEGRIDGPU_ITERATOR_SUB_HPP_ */
