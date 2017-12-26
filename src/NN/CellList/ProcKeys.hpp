/*
 * ProcKeys.hpp
 *
 *  Created on: Mar 14, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_PROCKEYS_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_PROCKEYS_HPP_

extern "C"
{
#include "hilbertKey.h"
}

/* !Brief Class for a linear (1D-like) order processing of cell keys for CellList_gen implementation
 *
 * \tparam dim Dimansionality of the space
 */
template<unsigned int dim, typename CellList>
class Process_keys_lin
{
	//! stub object
	openfpm::vector<size_t> keys;

public:

	//! Particle Iterator produced by this key generator
	typedef ParticleIt_CellP<CellList> Pit;

	/*! \brief Return cellkeys vector
	 *
	 * \return vector of cell keys
	 *
	 */
	inline const openfpm::vector<size_t> & getKeys() const
	{
		return keys;
	}

	/*! \brief Get a linear (1D-like) key from the coordinates and add to the getKeys vector
	 *
	 * \tparam S Cell list type
	 *
	 * \param obj Cell list object
	 * \param gk grid key
	 * \param m order of a curve
	 */
	template<typename S> void get_hkey(S & obj, grid_key_dx<dim> gk, size_t m)
	{
		size_t point[dim];

		for (size_t i = 0; i < dim; i++)
		{
			point[i] = gk.get(i) + obj.getPadding(i);
		}

		keys.add(obj.getGrid().LinIdPtr(static_cast<size_t *>(point)));
	}

	template<typename S> void linearize_hkeys(S & obj, size_t m)
	{
		return;
	}
};

/*! \brief Class for an hilbert order processing of cell keys for CellList_gen implementation
 *
 * \tparam dim Dimansionality of the space
 */
template<unsigned int dim, typename CellList>
class Process_keys_hilb
{
	//! vector for storing the cell keys
	openfpm::vector<size_t> keys;

	//! vector for storing the particle keys
	openfpm::vector<size_t> p_keys;

	//! Order of an hilbert curve
	size_t m = 0;

public:

	//! Particle Iterator produced by this key generator
	typedef ParticleIt_CellP<CellList> Pit;

	/*! \brief Return cellkeys vector
	 *
	 * \return vector of cell keys
	 *
	 */
	inline const openfpm::vector<size_t> & getKeys() const
	{
		return keys;
	}

	/*! \brief Get an hilbert key from the coordinates and add to the getKeys vector
	 *
	 * \tparam S Cell list type
	 *
	 * \param obj Cell list object
	 * \param gk grid key
	 * \param m order of a curve
	 */
	template<typename S> inline void get_hkey(S & obj, grid_key_dx<dim> gk, size_t m)
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

	/*! \brief Get get the coordinates from hilbert key, linearize and add to the getKeys vector
	 *
	 * \tparam S Cell list type
	 *
	 * \param obj Cell list object
	 * \param m order of a curve
	 */
	template<typename S> inline void linearize_hkeys(S & obj, size_t m)
	{
		//An integer to handle errors
		int err;

		//Array to handle output
		uint64_t coord[dim];
		size_t coord2[dim];

		keys.sort();

		openfpm::vector<size_t> keys_new;

		for(size_t i = 0; i < obj.getKeys().size(); i++)
		{
			getIntCoordFromHKey(coord, m, dim, obj.getKeys().get(i), &err);

			for (size_t j = 0 ; j < dim ; j++)	{coord2[j] = coord[j] + obj.getPadding(j);}

			keys_new.add(obj.getGrid().LinIdPtr(static_cast<size_t *>(coord2)));
		}

		keys.swap(keys_new);
	}
};

/*! \brief Class for an hilbert order processing of cell keys for CellList_gen implementation
 *
 * \tparam dim Dimansionality of the space
 */
template<unsigned int dim, typename CellList>
class Process_keys_grid
{
	//! vector for storing the cell keys
	openfpm::vector<size_t> keys;

	//! vector for storing the particle keys
	openfpm::vector<size_t> p_keys;

	//! Order of an hilbert curve
	size_t m = 0;

public:

	//! Particle Iterator produced by this key generator
	typedef ParticleIt_CellP<CellList> Pit;

	/*! \brief Return cellkeys vector
	 *
	 * \return vector of cell keys
	 *
	 */
	inline const openfpm::vector<size_t> & getKeys() const
	{
		return keys;
	}

	/*! \brief Get an hilbert key from the coordinates and add to the getKeys vector
	 *
	 * \tparam S Cell list type
	 *
	 * \param obj Cell list object
	 * \param gk grid key
	 * \param m order of a curve
	 */
	template<typename S> inline void get_hkey(S & obj, grid_key_dx<dim> gk, size_t m)
	{
		//An integer to handle errors
		int err;

		uint64_t point[dim];

		for (size_t i = 0; i < dim; i++)
		{
			point[i] = gk.get(i);
		}

		size_t hkey = obj.getGrid().LinId(gk);

		keys.add(hkey);
	}

	/*! \brief Get get the coordinates from hilbert key, linearize and add to the getKeys vector
	 *
	 * \tparam S Cell list type
	 *
	 * \param obj Cell list object
	 * \param m order of a curve
	 */
	template<typename S> inline void linearize_hkeys(S & obj, size_t m)
	{
	}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_PROCKEYS_HPP_ */
