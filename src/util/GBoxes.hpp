/*
 * Gboxes.hpp
 *
 *  Created on: May 2, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_UTIL_GBOXES_HPP_
#define OPENFPM_IO_SRC_UTIL_GBOXES_HPP_


/*! \brief This structure store the Box that define the domain inside the Ghost + domain box
 *
	\verbatim

                          (Ghost + Domain)
     +------------------+
     |                  |
     |  +------------+ <---------- (Domain)
     |  |            |  |
     |  |  Domain    |  |
     |  |  Box       |  |
     |  |            |  |
     |  |            |  |
     |  +------------+  |
     |                  |
     +------------------+
(0,0) local coordinate ---> ( x, y )

	\endverbatim

 *
 *  * Domain
 *
 * \tparam dim dimensionality
 *
 */
template<unsigned int dim>
struct GBoxes
{
	//! Ghost + Domain ghost
	Box<dim,long int> GDbox;
	//! Domain box
	Box<dim,long int> Dbox;
	//! origin of GDbox in global grid coordinates
	Point<dim,long int> origin;

	//! In case the grid is not defined everywhere but is defined by a set of boxes
	//! indicate from which box it come from
	size_t k;

	/*! \brief Indicate that this structure has no pointers inside
	 *
	 * \return true
	 *
	 */
	static bool noPointers()
	{
		return true;
	}
};


#endif /* OPENFPM_IO_SRC_UTIL_GBOXES_HPP_ */
