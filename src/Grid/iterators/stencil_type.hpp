/*
 * stencil_type.hpp
 *
 *  Created on: Jun 25, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_ITERATORS_STENCIL_TYPE_HPP_
#define OPENFPM_DATA_SRC_GRID_ITERATORS_STENCIL_TYPE_HPP_

/*! \brief Structure for stencil iterator
 *
 *
 * \tparam Np number of stencil points
 *
 */
template<unsigned int dim, unsigned int Np>
struct stencil_offset_compute
{
	//! Number of stencil points
	static const unsigned int nsp = Np;

	//! set of offsets for the stencil
	long int stencil_offset[Np];

	//! Stencil points
	grid_key_dx<dim> stencil_pnt[Np];


	/*! \brief Set the stencil points
	 *
	 * \param stencil points
	 *
	 */
	template<unsigned int dim2> void set_stencil(const grid_key_dx<dim2> (& stencil_pnt)[Np])
	{
		for (size_t i = 0 ; i < Np ; i++)
		{this->stencil_pnt[i] = stencil_pnt[i];}
	}

	/*! \brief Get the calculated stencil offset
	 *
	 * \tparam id stencil point
	 *
	 * \return the point
	 *
	 */
	template<unsigned int id> size_t getStencil() const
	{
		return stencil_offset[id];
	}

	/*! \brief Calculate the offsets of the stencil points
	 *
	 * \param g object storing the grid information
	 * \param start_p starting point where we calculate the stencil offsets
	 *
	 *
	 */
	template<unsigned int dim2, typename ginfo>
	inline void calc_offsets(const ginfo & g,
							 const grid_key_dx<dim2> & start_p)
	{
		for (size_t i = 0 ; i < Np ; i++)
		{
			grid_key_dx<dim2> zero;
			zero = start_p + stencil_pnt[i];

			stencil_offset[i] = g.LinId(zero);
		}
	}

	/*! \brief Calculate the offsets of the stencil points
	 *
	 * \tparam dim dimensionality
	 * \tparam ginfo grid information object
	 *
	 * \param g information of the grid
	 * \param start_p starting point
	 * \param stencil_pnt stencil points in the grid
	 *
	 *
	 */
	template<unsigned int dim2,typename ginfo>
	inline void calc_offsets(const ginfo & g,
							 const grid_key_dx<dim2> & start_p,
			                 const grid_key_dx<dim2> (& stencil_pnt)[Np])
	{
		for (size_t i = 0 ; i < Np ; i++)
		{
			grid_key_dx<dim2> offset_point;
			this->stencil_pnt[i] = stencil_pnt[i];
			offset_point = start_p + stencil_pnt[i];
			stencil_offset[i] = g.LinId(offset_point);
		}
	}

	/*! \brief Increment the offsets by one
	 *
	 *
	 */
	inline void increment()
	{
		// update the offsets
		for (size_t i = 0 ; i < Np ; i++)
			stencil_offset[i] += 1;
	}

	/*! \brief Adjust the offset
	 *
	 * \param i component
	 * \param idr component
	 * \param grid_base obbect containing the grid informations
	 *
	 */
	template<typename ginfo> inline void adjust_offset(size_t i, size_t idr, const ginfo & grid_base)
	{
		size_t str_dw = (i == 0)?1:grid_base.size_s(i-1);

		// update the offsets
		for (size_t k = 0 ; k < Np ; k++)
		{stencil_offset[k] += -str_dw*idr + grid_base.size_s(i);}
	}
};

/*! \brief no stencil
 *
 */
struct no_stencil
{
	//! dimensions of space, should be zero but is one 1
	//! otherwise can produce zero length arrays
	static const unsigned int nsp = 1;

	/*! \brief get stencil point
	 *
	 * do nothing
	 *
	 * \return 0
	 *
	 */
	template<unsigned int id> size_t getStencil() const
	{return 0;}


	/*! \brief Calculate the stencil offsets
	 *
	 * \tparam ginfo grid information object
	 * \tparam dim2 dimensionality of the starting point
	 *
	 * do nothing
	 *
	 * \param g grid information
	 * \param start_p starting point
	 *
	 */
	template<unsigned int dim2, typename ginfo>
	inline void calc_offsets(const ginfo & g,
							 const grid_key_dx<dim2> & start_p)
	{}

	/*! \brief Increment do nothing
	 *
	 * do nothing
	 *
	 */
	inline void increment()
	{}

	/*! \brief Set the stencil points
	 *
	 * \param stencil_pnt stencil points
	 *
	 */
	template<unsigned int dim2> void set_stencil(const grid_key_dx<dim2> (& stencil_pnt)[1])
	{
	}

	/*! \brief Adjust the offset
	 *
	 * do nothing
	 *
	 * \param i component
	 * \param idr previous component value
	 * \param grid_base information of the grid
	 *
	 */
	template<typename ginfo> inline void adjust_offset(size_t i, size_t idr, const ginfo & grid_base)
	{}
};



#endif /* OPENFPM_DATA_SRC_GRID_ITERATORS_STENCIL_TYPE_HPP_ */
