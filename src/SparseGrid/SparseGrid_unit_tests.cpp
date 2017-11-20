/*
 * SparseGrid_unit_tests.cpp
 *
 *  Created on: Oct 22, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_UNIT_TESTS_CPP_
#define OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_UNIT_TESTS_CPP_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "SparseGrid/SparseGrid.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include <math.h>

BOOST_AUTO_TEST_SUITE( sparse_grid_test )

BOOST_AUTO_TEST_CASE( sparse_grid_use_test)
{
	size_t sz[3] = {10000,10000,10000};

	sgrid_cpu<3,aggregate<float>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	// We fill a sphere with a band

	grid_key_dx<3> key1({5000,5000,5000});
	grid_key_dx<3> key2({5001,5001,5001});
	grid_key_dx<3> key3({5002,5003,5003});

	grid.template insert<0>(key1) = 1.0;
	grid.template insert<0>(key2) = 2.0;
	grid.template insert<0>(key3) = 3.0;

	BOOST_REQUIRE_EQUAL(grid.template get<0>(key1),1.0);
	BOOST_REQUIRE_EQUAL(grid.template get<0>(key2),2.0);
	BOOST_REQUIRE_EQUAL(grid.template get<0>(key3),3.0);

	auto it = grid.getDomainIterator();

	size_t count = 0;

	while (it.isNext())
	{
		count++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(count,(size_t)3);
}

BOOST_AUTO_TEST_CASE( sparse_grid_fill_all_test)
{
	size_t sz[3] = {171,171,171};

	sgrid_cpu<3,aggregate<float>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	grid_sm<3,void> g_sm(sz);

	grid_key_dx_iterator<3> kit(g_sm);

	while (kit.isNext())
	{
		auto key = kit.get();

		grid.template insert<0>(key) = g_sm.LinId(key);

		++kit;
	}

	auto it = grid.getDomainIterator();

	size_t count = 0;

	bool match = true;

	while (it.isNext())
	{
		auto key = it.getKey();

		// return a grid_key_dx
		auto key_pos = it.getKeyF();

		match &= (grid.template get<0>(key) == g_sm.LinId(key_pos));

		count++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(count,(size_t)171*171*171);
	BOOST_REQUIRE_EQUAL(grid.size(),(size_t)171*171*171);
	BOOST_REQUIRE_EQUAL(match,true);

	// remove all points

	grid_key_dx_iterator<3> kit2(g_sm);

	while (kit2.isNext())
	{
		auto key = kit2.get();

		grid.remove(key);

		++kit2;
	}

	size_t tot = grid.size();
	BOOST_REQUIRE_EQUAL(tot,0ul);
}


BOOST_AUTO_TEST_CASE( sparse_grid_fill_sparse_test)
{
	size_t sz[3] = {500,500,500};

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	CellDecomposer_sm<3, float, shift<3,float>> cdsm;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	cdsm.setDimensions(domain, sz, 0);

	double r = 0.3;
	double omega = 0.0;
	double phi = 0.0;

	// 3D sphere

	for (r = 0.3 ; r < 0.35 ;r += 0.001)
	{
		for (omega = 0.0; omega < M_PI ; omega += 0.006)
		{
			for (phi = 0.0; phi < 2.0*M_PI ; phi += 0.006)
			{
				Point<3,float> p;

				p.get(0) = r*sin(omega)*sin(phi) + 0.5;
				p.get(1) = r*sin(omega)*cos(phi) + 0.5;
				p.get(2) = r*cos(omega) + 0.5;

				// convert point into grid point

				grid_key_dx<3> kd = cdsm.getCellGrid(p);

				grid.template insert<0>(kd) = sin(omega)*sin(omega)*sin(2*phi);
				grid.template insert<1>(kd) = 0;
			}
		}
	}

	for (r = 0.3 ; r < 0.35 ;r += 0.001)
	{
		for (omega = 0.0; omega < M_PI ; omega += 0.006)
		{
			for (phi = 0.0; phi < 2.0*M_PI ; phi += 0.006)
			{

				Point<3,float> p;

				p.get(0) = r*sin(omega)*sin(phi) + 0.5;
				p.get(1) = r*sin(omega)*cos(phi) + 0.5;
				p.get(2) = r*cos(omega) + 0.5;

				// convert point into grid point

				grid_key_dx<3> kd = cdsm.getCellGrid(p);


				if (grid.template get<0>(kd) == sin(omega)*sin(omega)*sin(2*phi))
				{grid.template insert<1>(kd) = 1;}

			}
		}
	}

	auto it = grid.getDomainIterator();

	bool match = true;

	while(it.isNext())
	{
		auto key = it.getKey();

		if (grid.template get<1>(key) == 0)
		{match = false;}

		++it;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	// remove the points

	for (r = 0.3 ; r < 0.35 ;r += 0.001)
	{
		for (omega = 0.0; omega < M_PI ; omega += 0.006)
		{
			for (phi = 0.0; phi < 2.0*M_PI ; phi += 0.006)
			{

				Point<3,float> p;

				p.get(0) = r*sin(omega)*sin(phi) + 0.5;
				p.get(1) = r*sin(omega)*cos(phi) + 0.5;
				p.get(2) = r*cos(omega) + 0.5;

				// convert point into grid point

				grid_key_dx<3> kd = cdsm.getCellGrid(p);


				grid.remove(kd);

			}
		}
	}

	size_t tot;
	tot = grid.size();

	BOOST_REQUIRE_EQUAL(tot,0ul);
}


BOOST_AUTO_TEST_CASE( sparse_grid_resize_test)
{
	size_t sz[3] = {500,500,500};

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid(sz);
	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid2(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	CellDecomposer_sm<3, float, shift<3,float>> cdsm;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	cdsm.setDimensions(domain, sz, 0);

	double r = 0.3;
	double omega = 0.0;
	double phi = 0.0;

	// 3D sphere

	for (r = 0.3 ; r < 0.35 ;r += 0.001)
	{
		for (omega = 0.0; omega < M_PI ; omega += 0.006)
		{
			for (phi = 0.0; phi < 2.0*M_PI ; phi += 0.006)
			{
				Point<3,float> p;

				p.get(0) = r*sin(omega)*sin(phi) + 0.5;
				p.get(1) = r*sin(omega)*cos(phi) + 0.5;
				p.get(2) = r*cos(omega) + 0.5;

				// convert point into grid point

				grid_key_dx<3> kd = cdsm.getCellGrid(p);

				grid.template insert<0>(kd) = sin(omega)*sin(omega)*sin(2*phi);
				grid.template insert<1>(kd) = 0;
				grid2.template insert<0>(kd) = sin(omega)*sin(omega)*sin(2*phi);
				grid2.template insert<1>(kd) = 0;
			}
		}
	}

	size_t sz_b[3] = {1024,1024,1024};
	grid2.resize(sz_b);

	// Check that both grid contain the same information

	auto it = grid2.getDomainIterator();

	bool match = true;

	while(it.isNext())
	{
		auto key = it.getKey();

		if (grid.template get<0>(key) != grid2.template get<0>(key))
		{match = false;}

		++it;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	// now we resize smalle

	size_t sz_s[3] = {250,250,250};
	grid2.resize(sz_s);

	//

	auto it2 = grid.getDomainIterator();

	match = true;

	while(it2.isNext())
	{
		auto key = it2.getKeyF();

		// we check if the key is inside

		bool cin = true;

		cin &= (size_t)key.get(0) < sz_s[0];
		cin &= (size_t)key.get(1) < sz_s[1];
		cin &= (size_t)key.get(2) < sz_s[2];


		if (cin == true)
		{
			if (grid.template get<0>(key) != grid2.template get<0>(key))
			{match = false;}
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}



BOOST_AUTO_TEST_CASE( sparse_grid_fill_all_with_resize_test)
{
	size_t sz[3] = {10,10,171};

	sgrid_cpu<3,aggregate<float>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	grid_sm<3,void> g_sm(sz);

	grid_key_dx_iterator<3> kit(g_sm);

	while (kit.isNext())
	{
		auto key = kit.get();

		grid.template insert<0>(key) = g_sm.LinId(key);

		++kit;
	}

	size_t sz_b[3] = {20,20,200};

	// now we increase the size
	grid.resize(sz_b);

	auto it = grid.getDomainIterator();

	size_t count = 0;

	bool match = true;

	while (it.isNext())
	{
		auto key = it.getKey();

		// return a grid_key_dx
		auto key_pos = it.getKeyF();

		match &= (grid.template get<0>(key) == g_sm.LinId(key_pos));

		count++;

		++it;
	}

	BOOST_REQUIRE_EQUAL(count,(size_t)10*10*171);
	BOOST_REQUIRE_EQUAL(grid.size(),(size_t)10*10*171);
	BOOST_REQUIRE_EQUAL(match,true);

	// refill with the full set of point

	grid_sm<3,void> g_sm2(sz_b);

	grid_key_dx_iterator<3> kit2(g_sm2);

	while (kit2.isNext())
	{
		auto key = kit2.get();

		grid.template insert<0>(key) = g_sm2.LinId(key);

		++kit2;
	}

	auto it2 = grid.getDomainIterator();

	count = 0;

	match = true;

	while (it2.isNext())
	{
		auto key = it2.getKey();

		// return a grid_key_dx
		auto key_pos = it2.getKeyF();

		match &= (grid.template get<0>(key) == g_sm2.LinId(key_pos));

		count++;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(count,(size_t)20*20*200);
	BOOST_REQUIRE_EQUAL(grid.size(),(size_t)20*20*200);
	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( sparse_grid_insert_o_test)
{
	size_t sz[3] = {10,10,171};

	sgrid_cpu<3,aggregate<float>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	size_t ele;
	grid_key_dx<3> key({5,5,90});

	auto & flt = grid.insert_o(key,ele).template get<0>();
	flt[ele] = 117.0;

	BOOST_REQUIRE_EQUAL(grid.template get<0>(key),117.0);
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRID_UNIT_TESTS_CPP_ */
