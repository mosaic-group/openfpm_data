/*
 *  SFCKeys.hpp
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

/* !Brief Class for a linear (1D-like) order processing of cell keys for CellList implementation
 *
 * \tparam dim Dimansionality of the space
 */
template<unsigned int dim>
class SFCKeysLinear
{
	//! stub object
	openfpm::vector<size_t> keys;

public:
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
	 * \tparam CellList_type Cell list type
	 *
	 * \param cellList Cell list object
	 * \param gridKey grid key
	 * \param m order of a curve
	 */
	template<typename CellList_type> void get_hkey(CellList_type & cellList, grid_key_dx<dim> gridKey, size_t m)
	{
		size_t point[dim];

		for (size_t i = 0; i < dim; i++)
		{
			point[i] = gridKey.get(i) + cellList.getPadding(i);
		}

		keys.add(cellList.getGrid().LinIdPtr(static_cast<size_t *>(point)));
	}

	template<typename CellList_type> void linearize_hkeys(CellList_type & cellList, size_t m)
	{
		return;
	}
};

/*! \brief Class for an hilbert order processing of cell keys for CellList implementation
 *
 * \tparam dim Dimansionality of the space
 */
template<unsigned int dim>
class SFCKeysHilbert
{
	//! vector for storing the cell keys
	openfpm::vector<size_t> keys;

	//! Order of an hilbert curve
	size_t m = 0;

public:
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
	 * \tparam CellList_type Cell list type
	 *
	 * \param cellList Cell list object
	 * \param gridKey grid key
	 * \param m order of a curve
	 */
	template<typename CellList_type> inline void get_hkey(CellList_type & cellList, grid_key_dx<dim> gridKey, size_t m)
	{
		//An integer to handle errors
		int err;

		uint64_t point[dim];

		for (size_t i = 0; i < dim; i++)
		{
			point[i] = gridKey.get(i);
		}

		size_t hkey = getHKeyFromIntCoord(m, dim, point, &err);

		keys.add(hkey);
	}

	/*! \brief Get get the coordinates from hilbert key, linearize and add to the getKeys vector
	 *
	 * \tparam CellList_type Cell list type
	 *
	 * \param cellList Cell list object
	 * \param m order of a curve
	 */
	template<typename CellList_type> inline void linearize_hkeys(CellList_type & cellList, size_t m)
	{
		//An integer to handle errors
		int err;

		//Array to handle output
		uint64_t coord[dim];
		size_t coord2[dim];

		keys.sort();

		openfpm::vector<size_t> keys_new;

		for(size_t i = 0; i < keys.size(); i++)
		{
			getIntCoordFromHKey(coord, m, dim, keys.get(i), &err);

			for (size_t j = 0 ; j < dim ; j++)	{coord2[j] = coord[j] + cellList.getPadding(j);}

			keys_new.add(cellList.getGrid().LinIdPtr(static_cast<size_t *>(coord2)));
		}

		keys.swap(keys_new);
	}
};

/*! \brief Class for an hilbert order processing of cell keys for CellList implementation
 *
 * \tparam dim Dimansionality of the space
 */
template<unsigned int dim>
class SFCKeysGrid
{
	//! vector for storing the cell keys
	openfpm::vector<size_t> keys;

	//! Order of an hilbert curve
	size_t m = 0;

public:
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
	 * \tparam CellList_type Cell list type
	 *
	 * \param cellList Cell list object
	 * \param gridKey grid key
	 * \param m order of a curve
	 */
	template<typename CellList_type> inline void get_hkey(CellList_type & cellList, grid_key_dx<dim> gridKey, size_t m)
	{
		uint64_t point[dim];

		for (size_t i = 0; i < dim; i++)
		{
			point[i] = gridKey.get(i);
		}

		size_t hkey = cellList.getGrid().LinId(gridKey);

		keys.add(hkey);
	}

	/*! \brief Get get the coordinates from hilbert key, linearize and add to the getKeys vector
	 *
	 * \tparam CellList_type Cell list type
	 *
	 * \param cellList Cell list object
	 * \param m order of a curve
	 */
	template<typename CellList_type> inline void linearize_hkeys(CellList_type & cellList, size_t m)
	{
	}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_PROCKEYS_HPP_ */
