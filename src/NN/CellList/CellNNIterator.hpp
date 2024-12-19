/*
 * CellNNIterator.hpp
 *
 *  Created on: Mar 26, 2015
 *      Author: i-bird
 */

#ifndef CELLNNITERATOR_FULL_HPP_
#define CELLNNITERATOR_FULL_HPP_

#include "util/mathutil.hpp"
#include "NN/CellList/NNc_array.hpp"

#define FULL openfpm::math::pow(3,dim)
#define SYM  openfpm::math::pow(3,dim)/2 + 1
#define CRS openfpm::math::pow(2,dim)


/*! \brief Iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells
 *
 * \note to calculate quantities that involve a total reduction (like energies) use the CellIteratorSymRed
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 *
 */
template<unsigned int dim, typename Cell,int NNc_size>
class CellNNIterator
{
protected:

	//! actual element id
	const typename Cell::Mem_type_type::local_index_type * start_id;

	//! stop id to read the end of the cell
	const typename Cell::Mem_type_type::local_index_type * stop_id;

	//! Actual NNc_id;
	size_t NNc_id;

	//! Center cell, or cell for witch we are searching the NN-cell
	const long int cell;

	//! actual cell id = NNc[NNc_id]+cell stored for performance reason
	size_t cell_id;

	//! Cell list
	Cell & cl;

	//! NN cell id
	const NNc_array<dim,NNc_size> & NNc;

	/*! \brief Select non-empty cell
	 *
	 */
	__attribute__((always_inline)) inline void selectValid()
	{
		while (start_id == stop_id)
		{
			NNc_id++;

			// No more Cell
			if (NNc_id >= NNc_size) return;

			cell_id = NNc[NNc_id] + cell;

			start_id = &cl.getStartId(cell_id);
			stop_id = &cl.getStopId(cell_id);
		}
	}

private:


public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 *
	 */
	__attribute__((always_inline)) inline CellNNIterator(size_t cell, const NNc_array<dim,NNc_size> &NNc, Cell & cl)
	:NNc_id(0),cell(cell),cell_id(NNc[NNc_id] + cell),cl(cl),NNc(NNc)
	{
		start_id = &cl.getStartId(cell_id);
		stop_id = &cl.getStopId(cell_id);
		selectValid();
	}

	/*! \brief Check if there is the next element
	 *
	 * \return true if there is the next element
	 *
	 */
	__attribute__((always_inline)) inline bool isNext()
	{
		if (NNc_id >= NNc_size)
			return false;
		return true;
	}

	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	__attribute__((always_inline)) inline CellNNIterator & operator++()
	{
		start_id++;

		selectValid();

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	__attribute__((always_inline)) inline const typename Cell::Mem_type_type::local_index_type & get()
	{
		return cl.get_lin(start_id);
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	__attribute__((always_inline)) inline const typename Cell::Mem_type_type::local_index_type & get() const
	{
		return cl.get_lin(start_id);
	}
};


/*! \brief Symmetric local iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells.
 *
 * \note if we query the neighborhood of p and q is the neighborhood of p
 *          when we will query the neighborhood of q p is not present. This is
 *          useful to implement formula like \f$  \sum_{q = neighborhood(p) and p <= q} \f$
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 *
 */
template<unsigned int dim, typename Cell, typename vector_pos_type,int NNc_size>
class CellNNIteratorSymLocal
{
	//! index of the particle p
	size_t p;

	//! Position of the particle p
	const vector_pos_type & v;

	//! actual element id
	const typename Cell::Mem_type_type::local_index_type * start_id;

	//! stop id to read the end of the cell
	const typename Cell::Mem_type_type::local_index_type * stop_id;

	//! Actual NNc_id;
	size_t NNc_id;

	//! Center cell, or cell for witch we are searching the NN-cell
	const long int cell;

	//! actual cell id = NNc[NNc_id]+cell stored for performance reason
	size_t cell_id;

	//! Cell list
	Cell & cl;

	//! NN cell id
	const NNc_array<dim,NNc_size> & NNc;

	/*! Select the next valid element
	 *
	 */
	__attribute__((always_inline)) inline void IterateOwnCell()
	{
		while (start_id < stop_id)
		{

			size_t q = cl.get_lin(start_id);

			for (long int i = dim-1 ; i >= 0 ; i--)
			{
				if (v.template get<0>(p)[i] < v.template get<0>(q)[i])
					return;
				else if (v.template get<0>(p)[i] > v.template get<0>(q)[i])
					goto next;
			}

			if (q >= p)
				return;
next:
			start_id++;
		}
	}


	/*! Select the next valid element
	 *
	 */
	__attribute__((always_inline)) inline void selectValid()
	{
		if (NNc[NNc_id] == 0)
			IterateOwnCell();

		while (start_id == stop_id)
		{
			// No more Cell
			if (++NNc_id >= NNc_size) return;

			cell_id = NNc[NNc_id] + cell;

			if (NNc[NNc_id] == 0)
			{
				start_id = &(cl.getStartId(cell_id));
				stop_id = &(cl.getStopId(cell_id));
				IterateOwnCell();

			} else if (NNc[NNc_id] < 0) {
				start_id = &(cl.getGhostId(cell_id));
				stop_id = &(cl.getStopId(cell_id));
			} else {
				start_id = &(cl.getStartId(cell_id));
				stop_id = &(cl.getStopId(cell_id));
			}
		}
	}

public:
	__attribute__((always_inline)) inline const typename Cell::Mem_type_type::local_index_type & get()
	{
		return cl.get_lin(start_id);
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	__attribute__((always_inline)) inline const typename Cell::Mem_type_type::local_index_type & get() const
	{
		return cl.get_lin(start_id);
	}


	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param p index of the particle from which we are searching the neighborhood particles
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 *
	 */
	__attribute__((always_inline)) inline CellNNIteratorSymLocal(size_t cell, size_t p,
		const NNc_array<dim,NNc_size> &NNc, Cell & cl,
		const vector_pos_type & v)
	: NNc_id(0),cell(cell),
	cell_id(NNc[NNc_id] + cell),
	cl(cl),NNc(NNc),
	p(p),v(v)
	{
		if (NNc[NNc_id] < 0) {
			start_id = &(cl.getGhostId(cell_id));
			stop_id = &(cl.getStopId(cell_id));
		} else {
			start_id = &(cl.getStartId(cell_id));
			stop_id = &(cl.getStopId(cell_id));
		}

		selectValid();
	}

	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	__attribute__((always_inline)) inline CellNNIteratorSymLocal<dim,Cell,vector_pos_type,NNc_size> & operator++()
	{
		start_id++;

		selectValid();

		return *this;
	}

	/*! \brief Check if there is the next element
	 *
	 * \return true if there is the next element
	 *
	 */
	__attribute__((always_inline)) inline bool isNext()
	{
		if (NNc_id >= NNc_size)
			return false;
		return true;
	}
};

/*! \brief Symmetric iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells.
 *
 * \note if we query the neighborhood of p and q is the neighborhood of p
 *          when we will query the neighborhood of q p is not present. This is
 *          useful to implement formula like \f$  \sum_{q = neighborhood(p) and p <= q} \f$
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 *
 */
template<unsigned int dim, typename Cell, typename vector_pos_type,int NNc_size>
class CellNNIteratorSym : public CellNNIterator<dim,Cell,NNc_size>
{
	//! index of the particle p
	size_t p;

	//! Position of the particle p
	const vector_pos_type & v;

	/*! Select the next valid element
	 *
	 */
	__attribute__((always_inline)) inline void selectValid()
	{
		if (this->NNc[this->NNc_id] == 0)
		{
			while (this->start_id < this->stop_id)
			{
				size_t q = this->cl.get_lin(this->start_id);
				for (long int i = dim-1 ; i >= 0 ; i--)
				{
					if (v.template get<0>(p)[i] < v.template get<0>(q)[i])
						return;
					else if (v.template get<0>(p)[i] > v.template get<0>(q)[i])
						goto next;
				}
				if (q >= p)	return;
next:
				this->start_id++;
			}

			CellNNIterator<dim,Cell,NNc_size>::selectValid();
		}
		else
		{
			CellNNIterator<dim,Cell,NNc_size>::selectValid();
		}
	}

public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param p index of the particle from which we are searching the neighborhood particles
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 *
	 */
	__attribute__((always_inline)) inline CellNNIteratorSym(size_t cell, size_t p, const NNc_array<dim,NNc_size> &NNc, Cell & cl, const vector_pos_type & v)
	:CellNNIterator<dim,Cell,NNc_size>(cell,NNc,cl),p(p),v(v)
	{
		selectValid();
	}


	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	__attribute__((always_inline)) inline CellNNIteratorSym<dim,Cell,vector_pos_type,NNc_size> & operator++()
	{
		this->start_id++;

		selectValid();

		return *this;
	}
};

/*! \brief Symmetric iterator for the neighborhood of the cell structures
 *
 * In general you never create it directly but you get it from the CellList structures
 *
 * It iterate across all the element of the selected cell and the near cells.
 *
 * \note if we query the neighborhood of p and q is the neighborhood of p
 *          when we will query the neighborhood of q p is not present. This is
 *          useful to implement formula like \f$  \sum_{q = neighborhood(p) and p <= q} \f$
 *
 * \tparam dim dimensionality of the space where the cell live
 * \tparam Cell cell type on which the iterator is working
 * \tparam NNc_size neighborhood size
 *
 */
template<unsigned int dim, typename Cell, typename vector_pos_type , int NNc_size>
class CellNNIteratorSymMP : public CellNNIterator<dim,Cell,NNc_size>
{
	//! index of the particle p
	size_t p;

	//! Phase vector for particle p
	const vector_pos_type & v_p1;

	//! Phase vector for particle q
	const vector_pos_type & v_p2;

	/*! Select the next valid element
	 *
	 */
	__attribute__((always_inline)) inline void selectValid()
	{
		if (this->NNc[this->NNc_id] == 0)
		{
			while (this->start_id < this->stop_id)
			{
				size_t q = this->cl.get_lin(this->start_id);
				for (long int i = dim-1 ; i >= 0 ; i--)
				{
					if (v_p1.template get<0>(p)[i] < v_p2.template get<0>(q)[i])
						return;
					else if (v_p1.template get<0>(p)[i] > v_p2.template get<0>(q)[i])
						goto next;
				}
				if (q >= p)	return;
next:
				this->start_id++;
			}

			CellNNIterator<dim,Cell,NNc_size>::selectValid();
		}
		else
		{
			CellNNIterator<dim,Cell,NNc_size>::selectValid();
		}
	}

public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param p index of the particle from which we are searching the neighborhood particles
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 *
	 */
	__attribute__((always_inline)) inline CellNNIteratorSymMP(size_t cell,
			                   size_t p,
							   const NNc_array<dim,NNc_size> &NNc,
							   Cell & cl,
							   const vector_pos_type & v_p1,
							   const vector_pos_type & v_p2)
	:CellNNIterator<dim,Cell,NNc_size>(cell,NNc,cl),p(p),v_p1(v_p1),v_p2(v_p2)
	{
		selectValid();
	}


	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	__attribute__((always_inline)) inline CellNNIteratorSymMP<dim,Cell,vector_pos_type,NNc_size> & operator++()
	{
		this->start_id++;

		selectValid();

		return *this;
	}
};

/*! \brief it iterate through the elements of a cell
 *
 * In general you do not create this object you get it from the CellList structures
 *
 * \tparam Cell cell type
 *
 */
template<typename Cell> class CellIterator
{
	//! Cell list
	Cell & cl;

	//! actual element id inside the cell
	size_t ele_id;

	//! selected cell
	const long int cell;

public:

	/*! \brief Cell iterator
	 *
	 * \param cell Cell id
	 * \param cl Cell on which iterate
	 *
	 */
	__attribute__((always_inline)) inline CellIterator(const size_t cell, Cell & cl)
	:cl(cl),ele_id(0),cell(cell)
	{
	}

	/*! \brief Check if there is the next element
	 *
	 * \return true if there are still neighborhood particles
	 *
	 */
	__attribute__((always_inline)) inline bool isNext()
	{
		return cl.getNelements(cell) > ele_id;
	}

	/*! \brief take the next neoghborhood particle
	 *
	 * \return itself
	 *
	 */
	__attribute__((always_inline)) inline CellIterator & operator++()
	{
		ele_id++;

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	__attribute__((always_inline)) inline typename Cell::value_type & get()
	{
		return cl.get(cell,ele_id);
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	__attribute__((always_inline)) inline const typename Cell::value_type & get() const
	{
		return cl.get(cell,ele_id);
	}
};

#include "NN/CellList/CellNNIteratorRuntime.hpp"

#endif /* CELLNNITERATOR_FULL_HPP_ */
