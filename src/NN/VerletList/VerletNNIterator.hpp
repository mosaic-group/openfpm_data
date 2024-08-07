#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETNNITERATOR_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETNNITERATOR_HPP_


/*! \brief Iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 *
 */
    template<unsigned int dim, typename Ver> class VerletNNIterator
{
	//! start index for the neighborhood
	const typename Ver::Mem_type_type::local_index_type  * start;

	//! stop index for the neighborhood
	const typename Ver::Mem_type_type::local_index_type  * stop;

	//! actual neighborhood
	const typename Ver::Mem_type_type::local_index_type  * ele_id;

	//! verlet list
	Ver & ver;

public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param part_id Particle id
	 * \param ver Verlet-list
	 *
	 */
	inline VerletNNIterator(size_t part_id, Ver & ver)
	:start(&ver.getStart(part_id)),stop(&ver.getStop(part_id)),ver(ver)
	{ele_id = start;}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 * \return true if there is the next element
	 *
	 */
	inline bool isNext()
	{
		if (ele_id < stop)
			return true;
		return false;
	}

	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	inline VerletNNIterator & operator++()
	{
		ele_id++;

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline typename Ver::Mem_type_type::local_index_type get()
	{
		return ver.get_lin(ele_id);
	}


	/*! \brief Resets the iterator to the starting position
	 *
	 *
	 */
	inline void reset()
	{
		ele_id = start;
	}


};


#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETNNITERATOR_HPP_ */
