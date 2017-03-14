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

extern "C"
{
#include "hilbertKey.h"
}

/* !Brief Class for a linear (1D-like) order processing of cell keys for CellList_gen implementation
 *
 * \tparam dim Dimansionality of the space
 */
template<unsigned int dim>
class Process_keys_lin
{
public:

	/*! \brief Get a linear (1D-like) key from the coordinates and add to the getKeys vector
	 *
	 * \tparam S Cell list type
	 *
	 * \param obj Cell list object
	 * \param gk grid key
	 * \param m order of a curve
	 */
	template<typename S> static void get_hkey(S & obj, grid_key_dx<dim> gk, size_t m)
	{
		size_t point[dim];

		for (size_t i = 0; i < dim; i++)
		{
			point[i] = gk.get(i) + obj.getPadding(i);
		}

		obj.getKeys().add(obj.getGrid().LinIdPtr(static_cast<size_t *>(point)));
	}

	template<typename S> static void linearize_hkeys(S & obj, size_t m)
	{
		return;
	}
};

/* !Brief Class for an hilbert order processing of cell keys for CellList_gen implementation
 *
 * \tparam dim Dimansionality of the space
 */
template<unsigned int dim>
class Process_keys_hilb
{
public:

	/*! \brief Get an hilbert key from the coordinates and add to the getKeys vector
	 *
	 * \tparam S Cell list type
	 *
	 * \param obj Cell list object
	 * \param gk grid key
	 * \param m order of a curve
	 */
	template<typename S> static void get_hkey(S & obj, grid_key_dx<dim> gk, size_t m)
	{
		//An integer to handle errors
		int err;

		uint64_t point[dim];

		for (size_t i = 0; i < dim; i++)
		{
			point[i] = gk.get(i);
		}

		size_t hkey = getHKeyFromIntCoord(m, dim, point, &err);

		obj.getKeys().add(hkey);
	}

	/*! \brief Get get the coordinates from hilbert key, linearize and add to the getKeys vector
	 *
	 * \tparam S Cell list type
	 *
	 * \param obj Cell list object
	 * \param m order of a curve
	 */
	template<typename S> static void linearize_hkeys(S & obj, size_t m)
	{
		//An integer to handle errors
		int err;

		//Array to handle output
		uint64_t coord[dim];
		size_t coord2[dim];

		obj.getKeys().sort();

		openfpm::vector<size_t> keys_new;

		for(size_t i = 0; i < obj.getKeys().size(); i++)
		{
			getIntCoordFromHKey(coord, m, dim, obj.getKeys().get(i), &err);

			for (size_t j = 0 ; j < dim ; j++)	{coord2[j] = coord[j] + obj.getPadding(j);}

			keys_new.add(obj.getGrid().LinIdPtr(static_cast<size_t *>(coord2)));
		}

		obj.getKeys().swap(keys_new);
	}
};

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
template<unsigned int dim, typename T, typename Prock, typename Mem_type = Mem_fast<dim,T>, typename transform = no_transform<dim,T>, typename base=openfpm::vector<size_t>>
class CellList_gen : public CellList<dim,T,Mem_type,transform,base>, Prock
{
private:
	// Ghost marker
//	size_t g_m = 0;
/*
	// vector for storing the cell keys
	openfpm::vector<size_t> keys;

	// vector for storing the particle keys
	openfpm::vector<size_t> p_keys;

	//Order of an hilbert curve
	size_t m;*/

	/*! \brief Get an hilbert key from the coordinates and add to the getKeys vector
	 *
	 * \param gk grid key
	 *
	 */
/*	void get_hkey(grid_key_dx<dim> gk)
	{
		Prock::template get_hkey<CellList_gen<dim,T,Prock,Mem_type,transform,base>>(*this,gk,m);
	}*/

	/*! \brief Get get the coordinates from hilbert key, linearize and add to the getKeys vector
	 *
	 */
/*	void linearize_hkeys()
	{
		Prock::template linearize_hkeys<CellList_gen<dim,T,Prock,Mem_type,transform,base>>(*this,m);
	}*/

public:

	/*! Initialize the cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param g_m_new A ghost marker
	 * \param pad padding cell
	 * \param slot maximum number of slot
	 *
	 */
/*	void Initialize(const Box<dim,T> & box, const size_t (&div)[dim], size_t g_m_new, const size_t pad = 1, size_t slot=STARTING_NSLOT)
	{
		CellList<dim,T,Mem_type,transform,base>::Initialize(box,div,pad,slot);

		g_m = g_m_new;

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
			this->get_hkey(gk);

			++it;
		}

		// Sort and linearize keys
		this->linearize_hkeys();
	}*/

	CellList_gen()
	:CellList<dim,T,Mem_type,transform,base>()
	{};


	/*! \brief Return cellkeys vector
	 *
	 * \return vector of cell keys
	 *
	 */
/*	inline openfpm::vector<size_t> & getKeys()
	{
		return keys;
	}*/

	/*! \brief return the ghost marker
	 *
	 * \return ghost marker
	 *
	 */
/*	inline size_t get_gm()
	{
		return g_m;
	}*/

	/*! \brief Set the ghost marker
	 *
	 *
	 */
/*	inline void set_gm(size_t g_m)
	{
		this->g_m = g_m;
	}*/

	/*! \brief return the celllist iterator (across cells)
	 *
	 * \return an iterator
	 *
	 */
/*	inline Cell_list_iterator<CellList_gen<dim,T,Prock,Mem_type,transform,base>> getIterator()
	{
		return Cell_list_iterator<CellList_gen<dim,T,Prock,Mem_type,transform,base>>(*this);
	}*/
};


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTFAST_GEN_HPP_ */
