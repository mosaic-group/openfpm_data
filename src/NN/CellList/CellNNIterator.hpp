/*
 * CellNNIterator.hpp
 *
 *  Created on: Mar 26, 2015
 *      Author: i-bird
 */

#ifndef CELLNNITERATOR_FULL_HPP_
#define CELLNNITERATOR_FULL_HPP_

#include "util/mathutil.hpp"

#define FULL openfpm::math::pow(3,dim)
#define SYM  openfpm::math::pow(3,dim)/2 + 1
#define CRS openfpm::math::pow(2,dim)

#define NO_CHECK 1
#define SAFE 2

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
 * \tparam impl implementation specific options NO_CHECK do not do check on access, SAFE do check on access
 *
 */
template<unsigned int dim, typename Cell,unsigned int NNc_size, unsigned int impl>
class CellNNIterator
{
protected:

	//! actual element id
	size_t ele_id;

	//! Actual NNc_id;
	size_t NNc_id;

	//! Center cell, or cell for witch we are searching the NN-cell
	const long int cell;

	//! actual cell id = NNc[NNc_id]+cell stored for performance reason
	size_t cell_id;

	//! Cell list
	Cell & cl;

	/*! \brief Select non-empty cell
	 *
	 */
	inline void selectValid()
	{
		while (ele_id >= cl.getNelements(cell_id))
		{
			NNc_id++;

			// No more Cell
			if (NNc_id >= NNc_size) return;

			cell_id = NNc[NNc_id] + cell;

			ele_id = 0;
		}
	}

private:


	//! NN cell id
	const long int (& NNc)[NNc_size];

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
	inline CellNNIterator(size_t cell, const long int (&NNc)[NNc_size], Cell & cl)
	:ele_id(0),NNc_id(0),cell(cell),cell_id(NNc[NNc_id] + cell),cl(cl),NNc(NNc)
	{
		selectValid();
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	inline bool isNext()
	{
		if (NNc_id >= NNc_size)
			return false;
		return true;
	}

	/*! \brief take the next element
	 *
	 */
	inline CellNNIterator & operator++()
	{
		ele_id++;

		selectValid();

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline typename Cell::value_type & get()
	{
		return cl.get(cell_id,ele_id);
	}
};


/*! \brief Symmetric iterator for the neighborhood of the cell structures for reductions (like calculation of total energy)
 *
 * \note to calculate quantities per particles based (like forces) use the CellIteratorSym
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
 * \tparam impl implementation specific options NO_CHECK do not do check on access, SAFE do check on access
 *
 */
template<unsigned int dim, typename Cell,unsigned int NNc_size, unsigned int impl>
class CellNNIteratorSymRed : public CellNNIterator<dim,Cell,NNc_size,impl>
{
	size_t p;
	size_t g_m;

	const openfpm::vector<Point<dim,typename Cell::stype>> & v_pos;

	// Simulation domain (careful simulation not processor domain)
	const Box<dim,typename Cell::stype> & dom;

	inline void selectValid()
	{
		if (this->NNc_id == 0)
		{
			while (this->ele_id < this->cl.getNelements(this->cell_id))
			{
				size_t ele = this->cl.get(this->cell_id,this->ele_id);

				if (ele == 5954)
				{
					int debug = 0;
					debug++;
				}

				if (ele < g_m)
				{
					if (ele >= p)	return;
					this->ele_id++;
				}
				else
				{
					if (p > g_m)
					{
						this->ele_id++;
						continue;
					}

					bool ret = true;
	//				for (size_t i = 0 ; i < dim ; i++)
	//				{
						if (v_pos.template get<0>(ele)[2] < v_pos.template get<0>(p)[2])
						{
							this->ele_id++;
							ret = false;
						}
					/*	if (v_pos.template get<0>(ele)[1] < v_pos.template get<0>(p)[1])
						{
							this->ele_id++;
							ret = false;
							break;
						}
						if (v_pos.template get<0>(ele)[0] < v_pos.template get<0>(p)[0])
						{
							this->ele_id++;
							ret = false;
							break;
						}*/
	//				}
					if (ret == true)
						return;
				}
			}

			CellNNIterator<dim,Cell,NNc_size,impl>::selectValid();

			while (this->isNext() && p > g_m && this->cl.get(this->cell_id,this->ele_id) > g_m)
			{
				this->ele_id++;
				CellNNIterator<dim,Cell,NNc_size,impl>::selectValid();
			}
		}
		else
		{
			CellNNIterator<dim,Cell,NNc_size,impl>::selectValid();

			while (this->isNext() && p > g_m && this->cl.get(this->cell_id,this->ele_id) > g_m)
			{
				this->ele_id++;
				CellNNIterator<dim,Cell,NNc_size,impl>::selectValid();
			}


		}

		if (this->isNext() == true && p == 2146 && this->get() == 5954)
		{
			int debug = 0;
			debug++;
		}
	}

public:

	/*! \brief
	 *
	 * Cell NN iterator
	 *
	 * \param cell Cell id
	 * \param NNc Cell neighborhood indexes (relative)
	 * \param cl Cell structure
	 * \param v_pos vector of particle positions
	 * \param dom Simulation domain (carefull simulation domain not processor domain)
	 *
	 */
	inline CellNNIteratorSymRed(size_t cell, size_t p,
								const long int (&NNc)[NNc_size],
								Cell & cl,
								size_t g_m,
								const openfpm::vector<Point<dim,typename Cell::stype>> & v_pos,
								const Box<dim,typename Cell::stype> & dom)
	:CellNNIterator<dim,Cell,NNc_size,impl>(cell,NNc,cl),p(p),g_m(g_m),v_pos(v_pos),dom(dom)
	{
		selectValid();
	}


	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	inline CellNNIteratorSymRed<dim,Cell,NNc_size,impl> & operator++()
	{
		this->ele_id++;

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
 * \tparam impl implementation specific options NO_CHECK do not do check on access, SAFE do check on access
 *
 */
template<unsigned int dim, typename Cell,unsigned int NNc_size, unsigned int impl> class CellNNIteratorSym : public CellNNIterator<dim,Cell,NNc_size,impl>
{
	size_t p;

	inline void selectValid()
	{
		if (this->NNc_id == 0)
		{
			while (this->ele_id < this->cl.getNelements(this->cell_id))
			{
				if (this->cl.get(this->cell_id,this->ele_id) >= p)	return;
				this->ele_id++;
			}

			CellNNIterator<dim,Cell,NNc_size,impl>::selectValid();
		}
		else
		{
			CellNNIterator<dim,Cell,NNc_size,impl>::selectValid();
		}
	}

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
	inline CellNNIteratorSym(size_t cell, size_t p, const long int (&NNc)[NNc_size], Cell & cl)
	:CellNNIterator<dim,Cell,NNc_size,impl>(cell,NNc,cl),p(p)
	{
		selectValid();
	}


	/*! \brief take the next element
	 *
	 * \return itself
	 *
	 */
	inline CellNNIteratorSym<dim,Cell,NNc_size,impl> & operator++()
	{
		this->ele_id++;

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
	// Cell list
	Cell & cl;

	// actual element id inside the cell
	size_t ele_id;

	// selected cell
	const long int cell;

public:

	/*! \brief Cell iterator
	 *
	 * \param cell Cell id
	 * \param cl Cell on which iterate
	 *
	 */
	inline CellIterator(const size_t cell, Cell & cl)
	:cl(cl),ele_id(0),cell(cell)
	{
	}

	/*! \brief
	 *
	 * Check if there is the next element
	 *
	 */
	inline bool isNext()
	{
		return cl.getNelements(cell) > ele_id;
	}

	/*! \brief take the next element
	 *
	 */
	inline CellIterator & operator++()
	{
		ele_id++;

		return *this;
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline typename Cell::value_type & get()
	{
		return cl.get(cell,ele_id);
	}

	/*! \brief Get the value of the cell
	 *
	 * \return  the next element object
	 *
	 */
	inline const typename Cell::value_type & get() const
	{
		return cl.get(cell,ele_id);
	}
};

#endif /* CELLNNITERATOR_FULL_HPP_ */
