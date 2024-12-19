/*
 * CellListIterator_CRS.hpp
 *
 *  Created on: Nov 14, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_PARTICLEITCRS_CELLS_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_PARTICLEITCRS_CELLS_HPP_

#include "NN/CellList/CellNNIterator.hpp"
#include "NN/CellList/CellList_util.hpp"
#include "NN/CellList/NNc_array.hpp"

/*! \brief sub-sub-domain
 *
 * \tparam dim dimensionality
 *
 */
template<unsigned int dim>
struct subsub
{
	//! sub-sub-domain
	grid_key_dx<dim> subsub;

	//! Neighborhood of each sub-sub-domains
	openfpm::vector<grid_key_dx<dim>> NN_subsub;
};

/*! \brief Linearized version of subsub
 *
 * \tparam dim dimensionality
 *
 */
template<unsigned int dim>
struct subsub_lin
{
	//! sub-sub-domain
	size_t subsub;

	//! Neighborhood of each sub-sub-domains (indicate the relative position compared to subsub)
	openfpm::vector<long int> NN_subsub;
};

/*! \brief This iterator iterate across the particles of a Cell-list following the Cell structure
 *
 *
 * \tparam dim Dimensionality
 * \tparam CellListType type of the cell-list
 *
 */
template<unsigned int dim,typename CellListType, typename vector_pos_type> class ParticleItCRS_Cells
{
private:

	//! starting position
	const typename CellListType::Mem_type_type::loc_index * start;

	//! stop position
	const typename CellListType::Mem_type_type::loc_index * stop;

	//! Actual cell
	size_t cid;

	//! Neighborhood
	const long int * NNc;

	//! Neighborhood size
	long int NNc_size;

	//! If 0 we are iterating over the domain, if 1 we are iterating over the
	//! anomalous neighborhood cells, if 2 we terminate
	size_t dom_or_anom;

	//! List of all the domain cells
	const openfpm::vector<size_t> &dom_cell;

	//! List of all anomalous domain cells with neighborhood
	const openfpm::vector<subsub_lin<dim>> &anom_dom_cell;

	//! The array contain the neighborhood of the cell-id in case of symmetric interaction
	//
	//   * * *
	//     x *
	//
	const NNc_array<dim,openfpm::math::pow(3,dim)/2+1> & NNc_sym;

	//! Celllist type
	CellListType & cli;

	/*! \brief Adjust the counters to reach a valid particle element
	 *
	 *
	 */
	void selectValid()
	{
		while (start == stop)
		{
			cid++;

			size_t s_cell;
			if (dom_or_anom == 0)
			{
				if (cid >= dom_cell.size())
				{
					dom_or_anom = 1;
					cid = 0;

					// Terminate if we do not have anom cell
					if (anom_dom_cell.size() == 0)
					{
						dom_or_anom = 2;
						return;
					}

					s_cell = anom_dom_cell.get(cid).subsub;
				}
				else
					s_cell = dom_cell.get(cid);
			}
			else
			{
				if (cid >= anom_dom_cell.size())
				{
					dom_or_anom = 2;
					return;
				}

				s_cell = anom_dom_cell.get(cid).subsub;
			}

			// Get the starting particle
			start = &cli.getStartId(s_cell);

			// Get the stop particle
			stop = &cli.getStopId(s_cell);
		}
	}

public:

	/*! \brief Initialize the iterator
	 *
	 * \param cli Cell-list
	 * \param dom_cell domain cell
	 * \param anom_dom_cell anomalous domain cell
	 * \param NNc_sym symmetric neighborhood
	 *
	 *
	 */
	ParticleItCRS_Cells(CellListType & cli,
					 const openfpm::vector<size_t> & dom_cell,
					 const openfpm::vector<subsub_lin<dim>> & anom_dom_cell,
			         const NNc_array<dim,openfpm::math::pow(3,dim)/2+1> & NNc_sym)
	:cid(0),NNc(NULL),NNc_size(0),dom_or_anom(0),dom_cell(dom_cell),anom_dom_cell(anom_dom_cell),NNc_sym(NNc_sym),cli(cli)
	{
		size_t s_cell;
		if (dom_cell.size() != 0)
			s_cell = dom_cell.get(0);
		else if (anom_dom_cell.size() != 0)
		{
			s_cell = anom_dom_cell.get(0).subsub;
			dom_or_anom = 1;
		}
		else
		{
			dom_or_anom = 2;

			start = NULL;
			stop = NULL;

			return;
		}

		// Get the starting particle
		start = &cli.getStartId(s_cell);

		// Get the stop particle
		stop = &cli.getStopId(s_cell);

		selectValid();
	}

	/*! \brief Increment to the next particle
	 *
	 * \return The actual particle iterator
	 *
	 */
	ParticleItCRS_Cells & operator++()
	{
		++start;

		selectValid();

		return *this;
	}

	/*! \brief Return true if there is the next particle
	 *
	 * \return true if there is a new point
	 *
	 */
	bool isNext()
	{
		return dom_or_anom != 2;
	}

	/*! \brief Get the actual particle id
	 *
	 * \return the actual particle id
	 *
	 */
	size_t get()
	{
		return *start;
	}


	/*! \brief Get the neighborhood iterator according to the CRS scheme
	 *
	 * The CRS scheme use different neighborhood based on where the cell
	 * is positioned in the processor domain
	 *
	 *
	 * * * *
	 *   x *  for a cell x in the center of the domain
	 *
	 *  *
	 *    x   for a cell in the outside right
	 *
	 * \return Return an iterator over the neighborhood particles
	 *
	 */
	typename CellListType::SymNNIterator getNNIteratorCSR(const vector_pos_type & v) const
	{
		if (dom_or_anom == 0)
		{return typename CellListType::SymNNIterator(dom_cell.get(cid),*start,NNc_sym.getPointer(),openfpm::math::pow(3,dim)/2+1,cli,v);}
		else
		{return typename CellListType::SymNNIterator(anom_dom_cell.get(cid).subsub,
					                                    *start,
														&anom_dom_cell.get(cid).NN_subsub.get(0),
														anom_dom_cell.get(cid).NN_subsub.size(),
														cli,
														v);}
	}

	/*! \brief Get the neighborhood iterator according to the CRS scheme Multi-phase case
	 *
	 * The CRS scheme use different neighborhood based on where the cell
	 * is positioned in the processor domain
	 *
	 *
	 * * * *
	 *   x *  for a cell x in the center of the domain
	 *
	 *  *
	 *    x   for a cell in the outside right
	 *
	 * \return Return an iterator over the neighborhood particles
	 *
	 */
	typename CellListType::SymNNIterator getNNIteratorCSRM(const vector_pos_type & pos ,
														   const openfpm::vector<pos_v<vector_pos_type>> & v) const
	{
		if (dom_or_anom == 0)
			return typename CellListType::SymNNIterator(dom_cell.get(cid),CellListType::getV(*start),CellListType::getP(*start),NNc_sym.getPointer(),openfpm::math::pow(3,dim)/2+1,cli,pos,v);
		else
			return typename CellListType::SymNNIterator(anom_dom_cell.get(cid).subsub,CellListType::getV(*start),CellListType::getP(*start),&anom_dom_cell.get(cid).NN_subsub.get(0),anom_dom_cell.get(cid).NN_subsub.size(),cli,pos,v);
	}
};


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_PARTICLEITCRS_CELLS_HPP_ */
