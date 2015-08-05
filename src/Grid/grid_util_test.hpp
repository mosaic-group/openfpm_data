/*
 * grid_util_test.hpp
 *
 *  Created on: Jul 18, 2015
 *      Author: Pietro Incardona
 */

#ifndef SRC_GRID_GRID_UTIL_TEST_HPP_
#define SRC_GRID_GRID_UTIL_TEST_HPP_

/*! \brief Fill the grid with some data
 *
 * \param grid to fill
 *
 */
template<unsigned int dim, typename T> void fill_grid(T & grid)
{
	auto key_it = grid.getIterator();

	while (key_it.isNext())
	{
		grid_key_dx<dim> kk = key_it.get();

		grid.template get<P::x>(kk) = grid.getGrid().LinId(kk);
		grid.template get<P::y>(kk) = grid.getGrid().LinId(kk)+1;
		grid.template get<P::z>(kk) = grid.getGrid().LinId(kk)+2;
		grid.template get<P::s>(kk) = grid.getGrid().LinId(kk)+3;

		grid.template get<P::v>(kk)[0] = grid.getGrid().LinId(kk)+123;
		grid.template get<P::v>(kk)[1] = grid.getGrid().LinId(kk)+124;
		grid.template get<P::v>(kk)[2] = grid.getGrid().LinId(kk)+125;

		grid.template get<P::t>(kk)[0][0] = grid.getGrid().LinId(kk)+567;
		grid.template get<P::t>(kk)[0][1] = grid.getGrid().LinId(kk)+568;
		grid.template get<P::t>(kk)[0][2] = grid.getGrid().LinId(kk)+569;
		grid.template get<P::t>(kk)[1][0] = grid.getGrid().LinId(kk)+570;
		grid.template get<P::t>(kk)[1][1] = grid.getGrid().LinId(kk)+571;
		grid.template get<P::t>(kk)[1][2] = grid.getGrid().LinId(kk)+572;
		grid.template get<P::t>(kk)[2][0] = grid.getGrid().LinId(kk)+573;
		grid.template get<P::t>(kk)[2][1] = grid.getGrid().LinId(kk)+574;
		grid.template get<P::t>(kk)[2][2] = grid.getGrid().LinId(kk)+575;

		++key_it;
	}
}



#endif /* SRC_GRID_GRID_UTIL_TEST_HPP_ */
