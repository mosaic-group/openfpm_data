/*
 * grid_util_test.hpp
 *
 *  Created on: Jul 18, 2015
 *      Author: Pietro Incardona
 */

#ifndef SRC_GRID_GRID_UTIL_TEST_HPP_
#define SRC_GRID_GRID_UTIL_TEST_HPP_

#include "map_grid.hpp"
#include "Point_test.hpp"
#include "grid_key.hpp"
#include "Space/Shape/HyperCube.hpp"

/*! \brief Fill the grid with some data
 *
 * \param grid to fill
 *
 */
template<unsigned int dim, typename T> void fill_grid(T & grid)
{
	typedef Point_test<float> P;

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

template<unsigned int dim, typename g> void test_layout_gridNd(g & c3, size_t sz)
{
#ifdef VERBOSE_TEST
	std::cout << dim << "D Array with grid_key (without redundant dimension): " << "\n";

	timer t;
	t.start();
#endif

	//! [Access to an N-dimensional grid with an iterator]
	typedef Point_test<float> P;

	grid_key_dx_iterator<dim> key_it = c3.getIterator();

	while (key_it.isNext())
	{
		grid_key_dx<dim> kk = key_it.get();

		c3.template get<P::x>(kk) = 1.1f;
		c3.template get<P::y>(kk) = 1.2f;
		c3.template get<P::z>(kk) = 1.3f;
		c3.template get<P::s>(kk) = 1.0f;

		c3.template get<P::v>(kk)[0] = 1.0f;
		c3.template get<P::v>(kk)[1] = 2.0f;
		c3.template get<P::v>(kk)[2] = 3.0f;

		c3.template get<P::t>(kk)[0][0] = 1.0f;
		c3.template get<P::t>(kk)[0][1] = 2.0f;
		c3.template get<P::t>(kk)[0][2] = 3.0f;
		c3.template get<P::t>(kk)[1][0] = 4.0f;
		c3.template get<P::t>(kk)[1][1] = 5.0f;
		c3.template get<P::t>(kk)[1][2] = 6.0f;
		c3.template get<P::t>(kk)[2][0] = 7.0f;
		c3.template get<P::t>(kk)[2][1] = 8.0f;
		c3.template get<P::t>(kk)[2][2] = 9.0f;

		++key_it;

	}
	//! [Access to an N-dimensional grid with an iterator]

#ifdef VERBOSE_TEST
	t.stop();

	std::cout << "End : " << pow(sz,dim)*16*4/1024/1024 << " MB " << "  Bandwidth: " << pow(sz,dim)*16*4/1024/1024/t.getwct() << " MB/s  " << "\n";
#endif

	/////////////////////////////////// MEM CHECK ////////////////////////////////////////////////////////

	bool passed = true;

	key_it = c3.getIterator();

	while (key_it.isNext())
	{
		grid_key_dx<dim> kk = key_it.get();

		c3.template get<P::x>(kk) = c3.getGrid().LinId(kk);
		c3.template get<P::y>(kk) = c3.getGrid().LinId(kk)+1;
		c3.template get<P::z>(kk) = c3.getGrid().LinId(kk)+2;
		c3.template get<P::s>(kk) = c3.getGrid().LinId(kk)+3;

		c3.template get<P::v>(kk)[0] = c3.getGrid().LinId(kk)+123;
		c3.template get<P::v>(kk)[1] = c3.getGrid().LinId(kk)+124;
		c3.template get<P::v>(kk)[2] = c3.getGrid().LinId(kk)+125;

		c3.template get<P::t>(kk)[0][0] = c3.getGrid().LinId(kk)+567;
		c3.template get<P::t>(kk)[0][1] = c3.getGrid().LinId(kk)+568;
		c3.template get<P::t>(kk)[0][2] = c3.getGrid().LinId(kk)+569;
		c3.template get<P::t>(kk)[1][0] = c3.getGrid().LinId(kk)+570;
		c3.template get<P::t>(kk)[1][1] = c3.getGrid().LinId(kk)+571;
		c3.template get<P::t>(kk)[1][2] = c3.getGrid().LinId(kk)+572;
		c3.template get<P::t>(kk)[2][0] = c3.getGrid().LinId(kk)+573;
		c3.template get<P::t>(kk)[2][1] = c3.getGrid().LinId(kk)+574;
		c3.template get<P::t>(kk)[2][2] = c3.getGrid().LinId(kk)+575;

		++key_it;
	}


	key_it = c3.getIterator();

	while (key_it.isNext())
	{
		grid_key_dx<dim> kk = key_it.get();

		if (c3.template get<P::x>(kk) != c3.getGrid().LinId(kk)) passed = false;
		if (c3.template get<P::y>(kk) != c3.getGrid().LinId(kk)+1) passed = false;
		if (c3.template get<P::z>(kk) != c3.getGrid().LinId(kk)+2) passed = false;
		if (c3.template get<P::s>(kk) != c3.getGrid().LinId(kk)+3) passed = false;

		if (c3.template get<P::v>(kk)[0] != c3.getGrid().LinId(kk)+123) passed = false;
		if (c3.template get<P::v>(kk)[1] != c3.getGrid().LinId(kk)+124) passed = false;
		if (c3.template get<P::v>(kk)[2] != c3.getGrid().LinId(kk)+125) passed = false;

		if (c3.template get<P::t>(kk)[0][0] != c3.getGrid().LinId(kk)+567) passed = false;
		if (c3.template get<P::t>(kk)[0][1] != c3.getGrid().LinId(kk)+568) passed = false;
		if (c3.template get<P::t>(kk)[0][2] != c3.getGrid().LinId(kk)+569) passed = false;
		if (c3.template get<P::t>(kk)[1][0] != c3.getGrid().LinId(kk)+570) passed = false;
		if (c3.template get<P::t>(kk)[1][1] != c3.getGrid().LinId(kk)+571) passed = false;
		if (c3.template get<P::t>(kk)[1][2] != c3.getGrid().LinId(kk)+572) passed = false;
		if (c3.template get<P::t>(kk)[2][0] != c3.getGrid().LinId(kk)+573) passed = false;
		if (c3.template get<P::t>(kk)[2][1] != c3.getGrid().LinId(kk)+574) passed = false;
		if (c3.template get<P::t>(kk)[2][2] != c3.getGrid().LinId(kk)+575) passed = false;

		++key_it;
	}

	BOOST_REQUIRE_EQUAL(passed,true);

	// Check sub iterator

	/*
	 * Basically we first fill the interior part of the grid than the borders
	 * creating sub iterator always of smaller dimension
	 *
	 * Example:
	 *
	 * 2D
	 *
	 * if we have a square grid 16x16 we first will with 1 the interior 14x14 grid
	 *
	 * than the 4 line borders 14x1 with one
	 * than the 4 point borders 1x1
	 *
	 * We check that
	 *
	 * 1) The number of points for each sub-grid correspond
	 * 2) No point is filled more than one time
	 * 3) All point are filled
	 *
	 * we use the property x of c3
	 *
	 */

	// Erase the property x

	key_it = c3.getIterator();

	while (key_it.isNext())
	{
		grid_key_dx<dim> kk = key_it.get();

		c3.template get<P::x>(kk) = 0.0;

		++key_it;
	}

	for(size_t i = 0 ; i <= dim ; i++)
	{
		// get the combination of dimension dim-i
		std::vector<comb<dim>> combs = HyperCube<dim>::getCombinations_R(dim-i);

		// For each combination create a sub iterator

		for (size_t j = 0 ; j < combs.size() ; j++)
		{
			// Grid key of the sub-iterator

			grid_key_dx<dim> start;
			grid_key_dx<dim> stop;

			// sub iterator

			for (size_t k = 0 ; k < dim ; k++)
			{
				// if combination is 0 the hyper-cube

				if (combs[j].c[k] == -1)
				{
					start.set_d(k,0);
					stop.set_d(k,0);
				}
				else if (combs[j].c[k] == 1)
				{
					start.set_d(k,c3.getGrid().size(k)-1);
					stop.set_d(k,c3.getGrid().size(k)-1);
				}
				else
				{
					start.set_d(k,1);
					stop.set_d(k,c3.getGrid().size(k)-2);
				}
			}

			bool make_test = true;

#ifdef SE_CLASS1

			if (c3.size() == 0)
			{make_test = false;}

#endif

			if (make_test == true)
			{
				auto key_it = c3.getSubIterator(start,stop);

				while (key_it.isNext())
				{
					grid_key_dx<dim> kk = key_it.get();

					BOOST_REQUIRE_EQUAL(c3.template get<P::x>(kk),0.0);

					c3.template get<P::x>(kk) = 1.0;

					++key_it;
				}
			}
		}
	}

	// Check that everything is 1.0

	key_it = c3.getIterator();

	while (key_it.isNext())
	{
		grid_key_dx<dim> kk = key_it.get();

		BOOST_REQUIRE_EQUAL(c3.template get<P::x>(kk),1.0);

		++key_it;
	}
}

#endif /* SRC_GRID_GRID_UTIL_TEST_HPP_ */
