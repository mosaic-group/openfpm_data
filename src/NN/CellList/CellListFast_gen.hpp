/*
 * CellListFast_hilb.hpp
 *
 *  Created on: May 17, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTFAST_GEN_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTFAST_GEN_HPP_

#define HILBERT 1

#include "CellListFast.hpp"
#include "CellListIterator.hpp"

extern "C"
{
#include "hilbertKey.h"
}

template<unsigned int dim>
class Process_keys
{
public:
	static void get_hkey(grid_key_dx<dim> gk, size_t m, openfpm::vector<size_t> & keys)
	{

	}

	template<typename S> static void linearize_hkeys(size_t m, S & obj)
	{

	}
};

template<unsigned int dim>
class Process_keys_hilb
{
public:
	static void get_hkey(grid_key_dx<dim> gk, size_t m, openfpm::vector<size_t> & keys)
	{
		//An integer to handle errors
		int err;

		uint64_t point[dim];

		for (size_t i = 0; i < dim; i++)
		{
			point[i] = gk.get(i);
		}

		size_t hkey = getHKeyFromIntCoord(m, dim, point, &err);

		keys.add(hkey);
	}


	template<typename S> static void linearize_hkeys(size_t m, S & obj)
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

/* !Brief Class for FAST hilbert cell list implementation
 *
 * \see CellList<dim,T,FAST,transform,base>
 *
 * \tparam dim Dimansionality of the space
 * \tparam T type of the space float, double, complex
 * \tparam base Base structure that store the information
 *
 */

template<unsigned int dim, typename T, typename Prock, unsigned int impl=FAST, typename transform = no_transform<dim,T>, typename base=openfpm::vector<size_t>>
class CellList_gen : public CellList<dim,T,impl,transform,base>
{
private:
	// Ghost marker
	size_t g_m = 0;

	// vector for storing the cell keys
	openfpm::vector<size_t> keys;

	// vector for storing the particle keys
	openfpm::vector<size_t> p_keys;

	//Order of an hilbert curve
	size_t m;

	/*! \brief Get an hilbert key from the coordinates and add to the getKeys vector
	 *
	 * \param gk grid key
	 *
	 */
	void get_hkey(grid_key_dx<dim> gk)
	{
		Prock::get_hkey(gk,m,this->getKeys());
/*
		//An integer to handle errors
		int err;

		uint64_t point[dim];

		for (size_t i = 0; i < dim; i++)
		{
			point[i] = gk.get(i);
		}

		size_t hkey = getHKeyFromIntCoord(m, dim, point, &err);

		this->getKeys().add(hkey);
*/
	}

	/*! \brief Get get the coordinates from hilbert key, linearize and add to the getKeys vector
	 *
	 *
	 *
	 */
	void linearize_hkeys()
	{
		Prock::template linearize_hkeys<CellList_gen<dim,T,Prock,impl,transform,base>>(m,*this);
/*
		//An integer to handle errors
		int err;

		//Array to handle output
		uint64_t coord[dim];
		size_t coord2[dim];

		this->getKeys().sort();

		openfpm::vector<size_t> keys_new;

		for(size_t i = 0; i < this->getKeys().size(); i++)
		{
			getIntCoordFromHKey(coord, m, dim, this->getKeys().get(i), &err);

			for (size_t j = 0 ; j < dim ; j++)	{coord2[j] = coord[j] + this->getPadding(j);}

			keys_new.add(this->getGrid().LinIdPtr(static_cast<size_t *>(coord2)));
		}

		this->getKeys().swap(keys_new);
*/
	}

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
	void Initialize(const Box<dim,T> & box, const size_t (&div)[dim], size_t g_m_new, const size_t pad = 1, size_t slot=STARTING_NSLOT)
	{
		CellList<dim,T,impl,transform,base>::Initialize(box,div,pad,slot);

		g_m = g_m_new;

		size_t sz[dim];

		for (size_t i = 0; i < dim ; i++)
			sz[i] = this->getGrid().size(i) - 2*pad;

		grid_sm<dim,void> gs_small(sz);

		size_t a = gs_small.size(0);

		for (size_t i = 1 ; i < dim ; i++)
		{
			if (a < gs_small.size(i))
				a = gs_small.size(i);
		}

		for (m = 0; ; m++)
		{
			if ((1ul << m) >= a)
				break;
		}

		grid_key_dx_iterator<dim> it(gs_small);

		while (it.isNext())
		{
			auto gk = it.get();

			// Get an hilbert key of each cell and add to 'keys' vector
			this->get_hkey(gk);

			++it;
		}
		// Sort and linearize hilbert keys
		this->linearize_hkeys();
	}

	CellList_gen()
	:CellList<dim,T,impl,transform,base>(),m(0)
	{};


	/*! \brief Return cellkeys vector
	 *
	 * \return vector of cell keys
	 *
	 */
	inline openfpm::vector<size_t> & getKeys()
	{
		return keys;
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

	/*! \brief return the celllist iterator (across cells)
	 *
	 * \return an iterator
	 *
	 */
	inline Cell_list_iterator<CellList_gen<dim,T,Prock,impl,transform,base>> getIterator()
	{
		return Cell_list_iterator<CellList_gen<dim,T,Prock,impl,transform,base>>(*this);
	}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTFAST_GEN_HPP_ */
