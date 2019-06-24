/*
 * grid_test_utils.hpp
 *
 *  Created on: Jun 16, 2019
 *      Author: i-bird
 */

#ifndef GRID_TEST_UTILS_HPP_
#define GRID_TEST_UTILS_HPP_


template<typename grid_type>
void copy_test(grid_type & g_src, grid_type & g_dst,
			   Box<grid_type::dims,size_t> & box_src, Box<grid_type::dims,size_t> & box_dst)
{
	g_dst.setMemory();
	g_src.setMemory();

	auto itd = g_dst.getIterator();

	while (itd.isNext())
	{
		auto k = itd.get();

		g_dst.template get<0>(k) = 0;

		g_dst.template get<1>(k)[0] = 0;
		g_dst.template get<1>(k)[1] = 0;
		g_dst.template get<1>(k)[2] = 0;

		g_dst.template get<2>(k)[0][0] = 0;
		g_dst.template get<2>(k)[0][1] = 0;
		g_dst.template get<2>(k)[0][2] = 0;
		g_dst.template get<2>(k)[1][0] = 0;
		g_dst.template get<2>(k)[1][1] = 0;
		g_dst.template get<2>(k)[1][2] = 0;
		g_dst.template get<2>(k)[2][0] = 0;
		g_dst.template get<2>(k)[2][1] = 0;
		g_dst.template get<2>(k)[2][2] = 0;

		++itd;
	}

	auto & gs = g_src.getGrid();

	auto its = g_src.getIterator();

	while (its.isNext())
	{
		auto k = its.get();

		g_src.template get<0>(k) = gs.LinId(k);

		g_src.template get<1>(k)[0] = gs.LinId(k) + 100;
		g_src.template get<1>(k)[1] = gs.LinId(k) + 200;
		g_src.template get<1>(k)[2] = gs.LinId(k) + 300;

		g_src.template get<2>(k)[0][0] = gs.LinId(k) + 1000;
		g_src.template get<2>(k)[0][1] = gs.LinId(k) + 2000;
		g_src.template get<2>(k)[0][2] = gs.LinId(k) + 3000;
		g_src.template get<2>(k)[1][0] = gs.LinId(k) + 4000;
		g_src.template get<2>(k)[1][1] = gs.LinId(k) + 5000;
		g_src.template get<2>(k)[1][2] = gs.LinId(k) + 6000;
		g_src.template get<2>(k)[2][0] = gs.LinId(k) + 7000;
		g_src.template get<2>(k)[2][1] = gs.LinId(k) + 8000;
		g_src.template get<2>(k)[2][2] = gs.LinId(k) + 9000;

		++its;
	}

	// copy
	g_dst.copy_to(g_src,box_src,box_dst);


	// Check

	itd = g_dst.getIterator();

	while (itd.isNext())
	{
		auto k = itd.get();
		Point<grid_type::dims,size_t> p;

		for (size_t i = 0 ; i < grid_type::dims ; i++)
		{p.get(i) = k.get(i);}

		if (box_dst.isInside(p) == true)
		{
			grid_key_dx<grid_type::dims> ks = k + box_src.getKP1() - box_dst.getKP1();

			BOOST_REQUIRE_EQUAL(g_dst.template get<0>(k),gs.LinId(ks));

			BOOST_REQUIRE_EQUAL(g_dst.template get<1>(k)[0],gs.LinId(ks) + 100);
			BOOST_REQUIRE_EQUAL(g_dst.template get<1>(k)[1],gs.LinId(ks) + 200);
			BOOST_REQUIRE_EQUAL(g_dst.template get<1>(k)[2],gs.LinId(ks) + 300);

			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[0][0],gs.LinId(ks) + 1000);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[0][1],gs.LinId(ks) + 2000);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[0][2],gs.LinId(ks) + 3000);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[1][0],gs.LinId(ks) + 4000);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[1][1],gs.LinId(ks) + 5000);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[1][2],gs.LinId(ks) + 6000);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[2][0],gs.LinId(ks) + 7000);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[2][1],gs.LinId(ks) + 8000);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[2][2],gs.LinId(ks) + 9000);
		}
		else
		{
			BOOST_REQUIRE_EQUAL(g_dst.template get<0>(k),0);

			BOOST_REQUIRE_EQUAL(g_dst.template get<1>(k)[0],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<1>(k)[1],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<1>(k)[2],0);

			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[0][0],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[0][1],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[0][2],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[1][0],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[1][1],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[1][2],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[2][0],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[2][1],0);
			BOOST_REQUIRE_EQUAL(g_dst.template get<2>(k)[2][2],0);
		}

		++itd;
	}
}


#endif /* GRID_TEST_UTILS_HPP_ */
