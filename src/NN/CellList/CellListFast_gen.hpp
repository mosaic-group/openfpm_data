/*
 * CellListFast_hilb.hpp
 *
 *  Created on: May 17, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTFAST_GEN_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTFAST_GEN_HPP_

#define HILBERT 1

#include "CellList.hpp"
#include "ProcKeys.hpp"

/* \brief Cell list implementation with particle iterator over cells
 *
 * \see CellList<dim,T,FAST,transform,base>
 *
 * \tparam dim Dimensionality of the space
 * \tparam T type of the space float, double, complex
 * \tparam Prock indicate the cell iterator over the Cell
 * \tparam Mem_type indicate how the Cell-list are implemented in memory
 * \tparam base Base structure that store the information
 *
 */
template<unsigned int dim,
         typename T,
		 template <unsigned int, typename> class Prock,
		 typename Mem_type = Mem_fast<dim,T>,
		 typename transform = no_transform<dim,T>,
		 typename base=openfpm::vector<size_t>>
class CellList_gen : public CellList<dim,T,Mem_type,transform,base>
{
private:

	// Ghost marker
	size_t g_m = 0;

	//! It is an object that indicate which space filling curve to use for the
	//! iteration across cells
	Prock<dim,CellList_gen<dim,T,Prock,Mem_type,transform,base>> SFC;

	//! Init SFC
	bool init_sfc;

	void initialize_sfc(size_t pad)
	{
		size_t sz[dim];

		//Get grid_sm without padding (gs_small)
		for (size_t i = 0; i < dim ; i++)
			sz[i] = this->getGrid().size(i) - 2*pad;

		grid_sm<dim,void> gs_small(sz);

		size_t a = gs_small.size(0);

		for (size_t i = 1 ; i < dim ; i++)
		{
			if (a < gs_small.size(i))
				a = gs_small.size(i);
		}

		size_t m;

		//Calculate an hilberts curve order
		for (m = 0; ; m++)
		{
			if ((1ul << m) >= a)
				break;
		}

		grid_key_dx_iterator<dim> it(gs_small);

		while (it.isNext())
		{
			auto gk = it.get();

			// Get a key of each cell and add to 'keys' vector
			SFC.get_hkey(*this,gk,m);

			++it;
		}

		// Sort and linearize keys
		SFC.linearize_hkeys(*this,m);
	}

public:


	CellList_gen()
	:CellList<dim,T,Mem_type,transform,base>(),init_sfc(false)
	{};


	/*! \brief Set the ghost marker
	 *
	 *
	 */
	const Prock<dim,CellList_gen<dim,T,Prock,Mem_type,transform,base>> & getCellSFC() const
	{
		return SFC;
	}

	/*! \brief return the celllist iterator (across cells)
	 *
	 * \return an iterator
	 *
	 */
	inline typename Prock<dim,CellList_gen<dim,T,Prock,Mem_type,transform,base>>::Pit getIterator()
	{
		// Initialize SFC
		if (init_sfc == false)
		{
			initialize_sfc(this->getPadding(0));
			init_sfc = true;
		}

		return typename Prock<dim,CellList_gen<dim,T,Prock,Mem_type,transform,base>>::Pit(*this);
	}

	/*! Initialize the cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param g_m_new A ghost marker
	 * \param pad padding cell
	 * \param slot maximum number of slot
	 *
	 */
	void Initialize(const Box<dim,T> & box, const size_t (&div)[dim], size_t g_m_new, const size_t pad = 1, size_t slot=STARTING_NSLOT)
	{
		CellList<dim,T,Mem_type,transform,base>::Initialize(box,div,pad,slot);

		g_m = g_m_new;
	}


	/*! \brief Return cellkeys vector
	 *
	 * \return vector of cell keys
	 *
	 */
	inline const openfpm::vector<size_t> & getKeys()
	{
		return SFC.getKeys();
	}

	/*! \brief return the ghost marker
	 *
	 * \return ghost marker
	 *
	 */
	inline size_t get_gm()
	{
		return g_m;
	}

	/*! \brief Set the ghost marker
	 *
	 *
	 */
	inline void set_gm(size_t g_m)
	{
		this->g_m = g_m;
	}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTFAST_GEN_HPP_ */
