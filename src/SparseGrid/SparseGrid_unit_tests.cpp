/*
 * SparseGrid_unit_tests.cpp
 *
 *  Created on: Oct 22, 2017
 *      Author: i-bird
 */

#define DISABLE_MPI_WRITTERS

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "SparseGrid/SparseGrid.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include <math.h>
#include "util/debug.hpp"

BOOST_AUTO_TEST_SUITE( sparse_grid_test )

template <typename grid_type, typename cell_decomposer>
size_t fill_sphere(grid_type & grid, cell_decomposer & cdsm)
{
	size_t tot_count = 0;
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

	auto it = grid.getIterator();

	while (it.isNext())
	{
		tot_count++;

		++it;
	}

	return tot_count;
}

template <typename grid_type, typename cell_decomposer>
size_t fill_sphere_quad(grid_type & grid, cell_decomposer & cdsm)
{
	size_t tot_count = 0;
	double r = 0.3;
	double omega = 0.0;
	double phi = 0.0;

	// 3D sphere

	for (r = 0.3 ; r < 0.4 ;r += 0.001)
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

				grid.template insert<0>(kd) = kd.get(0)*kd.get(0) + kd.get(1)*kd.get(1) + kd.get(2)*kd.get(2);
				grid.template insert<1>(kd) = 0;
			}
		}
	}

	auto it = grid.getIterator();

	while (it.isNext())
	{
		tot_count++;

		++it;
	}

	return tot_count;
}

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

	auto it = grid.getIterator();

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

	auto it = grid.getIterator();

	size_t count = 0;

	bool match = true;

	while (it.isNext())
	{
		auto key = it.get();

		// return a grid_key_dx
		auto key_pos = it.getKeyF();

		match &= (grid.template get<0>(key_pos) == g_sm.LinId(key));

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

	fill_sphere(grid,cdsm);

	double r = 0.3;
	double omega = 0.0;
	double phi = 0.0;

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

	auto it = grid.getIterator();

	bool match = true;

	while(it.isNext())
	{
		auto key = it.get();

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

	auto it = grid2.getIterator();

	bool match = true;

	while(it.isNext())
	{
		auto key = it.get();

		if (grid.template get<0>(key) != grid2.template get<0>(key))
		{match = false;}

		++it;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	// now we resize smalle

	size_t sz_s[3] = {250,250,250};
	grid2.resize(sz_s);

	//

	auto it2 = grid.getIterator();

	match = true;

	while(it2.isNext())
	{
		auto key = it2.get();

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

	auto it = grid.getIterator();

	size_t count = 0;

	bool match = true;

	while (it.isNext())
	{
		auto key = it.get();

		// return a grid_key_dx
		auto key_pos = it.getKeyF();

		match &= (grid.template get<0>(key_pos) == g_sm.LinId(key));

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

	auto it2 = grid.getIterator();

	count = 0;

	match = true;

	while (it2.isNext())
	{
		auto key = it2.get();

		// return a grid_key_dx
		auto key_pos = it2.getKeyF();

		match &= (grid.template get<0>(key_pos) == g_sm2.LinId(key));

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


BOOST_AUTO_TEST_CASE( sparse_grid_sub_grid_it)
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

	grid_key_dx<3> start({21,21,21});
	grid_key_dx<3> stop({90,90,90});

	bool error = false;
	size_t count = 0;
	auto it_sub = grid.getIterator(start,stop);

	while (it_sub.isNext())
	{
		auto gkey = it_sub.get();

		if (gkey.get(0) < start.get(0) ||
			gkey.get(1) < start.get(1) ||
			gkey.get(2) < start.get(2) ||
			gkey.get(0) > stop.get(0) ||
			gkey.get(1) > stop.get(1) ||
			gkey.get(2) > stop.get(2))
		{
			error = true;
		}

		count++;

		++it_sub;
	}

	size_t tot = (stop.get(2) - start.get(2) + 1)*(stop.get(1) - start.get(1) + 1)*(stop.get(0) - start.get(0) + 1);
	BOOST_REQUIRE_EQUAL(error,false);
	BOOST_REQUIRE_EQUAL(count,tot);
}


BOOST_AUTO_TEST_CASE( sparse_grid_sub_grid_it_quarter_sphere)
{
	size_t sz[3] = {501,501,501};
	size_t sz_cell[3] = {500,500,500};

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	CellDecomposer_sm<3, float, shift<3,float>> cdsm;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	cdsm.setDimensions(domain, sz_cell, 0);

	fill_sphere(grid,cdsm);

	grid_key_dx<3> start({0,0,0});
	grid_key_dx<3> stop({250,250,250});

	bool error = false;
	size_t count = 0;
	auto it_sub = grid.getIterator(start,stop);

	while (it_sub.isNext())
	{
		auto gkey = it_sub.get();

		if (gkey.get(0) < start.get(0) ||
			gkey.get(1) < start.get(1) ||
			gkey.get(2) < start.get(2) ||
			gkey.get(0) > stop.get(0) ||
			gkey.get(1) > stop.get(1) ||
			gkey.get(2) > stop.get(2))
		{
			error = true;
		}

		// Check that the point is in the sphere

		double radius = (gkey.get(0) - 250)*(gkey.get(0) - 250) +
						(gkey.get(1) - 250)*(gkey.get(1) - 250) +
						(gkey.get(2) - 250)*(gkey.get(2) - 250);

		radius = sqrt(radius);

		if (radius < 150 || radius >= 175)
		{
			// if is not in the radius remove it
			grid.remove(gkey);
		}

		count++;

		++it_sub;
	}

	BOOST_REQUIRE_EQUAL(error,false);

	// We go again across the point now every point out the sphere is an error

	count = 0;
	auto it_sub2 = grid.getIterator(start,stop);

	while (it_sub2.isNext())
	{
		auto gkey = it_sub2.get();

		if (gkey.get(0) < start.get(0) ||
			gkey.get(1) < start.get(1) ||
			gkey.get(2) < start.get(2) ||
			gkey.get(0) > stop.get(0) ||
			gkey.get(1) > stop.get(1) ||
			gkey.get(2) > stop.get(2))
		{
			error = true;
		}

		// Check that the point is in the sphere

		double radius = (gkey.get(0) - 250)*(gkey.get(0) - 250) +
						(gkey.get(1) - 250)*(gkey.get(1) - 250) +
						(gkey.get(2) - 250)*(gkey.get(2) - 250);

		radius = sqrt(radius);

		if (radius < 150 || radius >= 175)
		{
			error = true;
		}

		count++;

		++it_sub2;
	}

	BOOST_REQUIRE_EQUAL(error,false);
}

BOOST_AUTO_TEST_CASE( sparse_grid_fast_stencil)
{
	size_t sz[3] = {501,501,501};
	size_t sz_cell[3] = {500,500,500};

	sgrid_cpu<3,aggregate<double,double,int>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	CellDecomposer_sm<3, float, shift<3,float>> cdsm;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	cdsm.setDimensions(domain, sz_cell, 0);

	fill_sphere_quad(grid,cdsm);

	grid_key_dx<3> start({1,1,1});
	grid_key_dx<3> stop({499,499,499});

	/////// Check

	timer t;
	t.start();


	for (int i = 0 ; i < 100 ; i++)
	{

	auto it = grid.getBlockIterator<1>(start,stop);

	unsigned char mask[decltype(it)::sizeBlockBord];
	double block_bord_src[decltype(it)::sizeBlockBord];
	double block_bord_dst[decltype(it)::sizeBlock];

	while (it.isNext())
	{
		it.loadBlockBorder<0>(block_bord_src,mask);

		for (int k = it.start_b(2) ; k < it.stop_b(2) ; k++)
		{
			for (int j = it.start_b(1) ; j < it.stop_b(1) ; j++)
			{
				for (int i = it.start_b(0) ; i < it.stop_b(0) ; i++)
				{
					int c = it.LinB(i,j,k);

					int xp = it.LinB(i+1,j,k);
					int xm = it.LinB(i-1,j,k);

					int yp = it.LinB(i,j+1,k);
					int ym = it.LinB(i,j-1,k);

					int zp = it.LinB(i,j,k+1);
					int zm = it.LinB(i,j,k-1);

					// we do only id exist the point
					if (mask[c] == false) {continue;}

					bool surround = mask[xp] & mask[xm] & mask[ym] & mask[yp] & mask[zp] & mask[zm];

					double Lap = block_bord_src[xp] + block_bord_src[xm] +
					block_bord_src[yp] + block_bord_src[ym] +
					block_bord_src[zp] + block_bord_src[zm] - 6.0*block_bord_src[c];

					block_bord_dst[it.LinB_off(i,j,k)] = (surround)?Lap:6.0;
				}
			}
		}

		it.storeBlock<1>(block_bord_dst);

		++it;
	}

	}

	t.stop();

	std::cout << "Laplacian " << t.getwct() << " points: " << grid.size() << "   " << grid.size_all() << std::endl;

	bool check = true;
	auto it2 = grid.getIterator();
	while (it2.isNext())
	{
		auto p = it2.get();

		check &= grid.template get<1>(p) == 6;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(check,true);
	// Check correct-ness

//	print_grid("debug_out",grid);
}

BOOST_AUTO_TEST_CASE( sparse_grid_slow_stencil)
{
	size_t sz[3] = {501,501,501};
	size_t sz_cell[3] = {500,500,500};

	sgrid_cpu<3,aggregate<double,double,int>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	CellDecomposer_sm<3, float, shift<3,float>> cdsm;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	cdsm.setDimensions(domain, sz_cell, 0);

	fill_sphere_quad(grid,cdsm);

	grid_key_dx<3> start({1,1,1});
	grid_key_dx<3> stop({499,499,499});

	/////// Check

	timer t;
	t.start();

	auto it = grid.getIterator();

	while (it.isNext())
	{
		// center point
		auto p = it.get();

		// plus,minus X,Y,Z
		auto mx = p.move(0,-1);
		auto px = p.move(0,1);
		auto my = p.move(1,-1);
		auto py = p.move(1,1);
		auto mz = p.move(2,-1);
		auto pz = p.move(2,1);

		bool surround = grid.existPoint(mx) & grid.existPoint(px) & grid.existPoint(py) & grid.existPoint(my) & grid.existPoint(mz) & grid.existPoint(pz);

		double Lap = grid.template get<0>(mz) + grid.template get<0>(pz) +
				     grid.template get<0>(my) + grid.template get<0>(py) +
				     grid.template get<0>(mx) + grid.template get<0>(px) - 6.0*grid.template get<0>(p);

		grid.template insert<1>(p) = (surround)?Lap:6.0;

		++it;
	}

	t.stop();

	std::cout << "Laplacian " << t.getwct() << " points: " << grid.size() << std::endl;

	bool check = true;
	auto it2 = grid.getIterator();
	while (it2.isNext())
	{
		auto p = it2.get();

		check &= grid.template get<1>(p) == 6;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(check,true);
	// Check correct-ness
	//grid.write("debug_out");
}

template<typename sgrid> void Test_unpack_and_check_full(sgrid & grid)
{
	grid_key_dx<3> end;
	grid_key_dx<3> zero({0,0,0});
	size_t req2 = 0;
	size_t req3 = 0;
	size_t sz[3];

	for (size_t i = 0 ; i < 3 ; i++)
	{
		end.set_d(i,grid.getGrid().size(i) - 1);
		sz[i] = 0;
	}

	grid.template packRequest<0>(req2);
	auto sub_it = grid.getIterator(zero,end);
	grid.template packRequest<0>(sub_it,req3);

	BOOST_REQUIRE_EQUAL(req2,req3);

	Pack_stat sts2;
	Pack_stat sts3;

	// allocate the memory 2
	HeapMemory pmem2;
	pmem2.allocate(req2);
	ExtPreAlloc<HeapMemory> & mem2 = *(new ExtPreAlloc<HeapMemory>(req2,pmem2));
	mem2.incRef();

	// allocate the memory 3
	HeapMemory pmem3;
	pmem3.allocate(req3);
	ExtPreAlloc<HeapMemory> & mem3 = *(new ExtPreAlloc<HeapMemory>(req3,pmem3));
	mem3.incRef();

	grid.template pack<0>(mem2,sts2);
	grid.template pack<0>(mem3,sub_it,sts3);

	BOOST_REQUIRE_EQUAL(mem2.size(),mem3.size());

	bool check = true;

	char * p2 = (char *)pmem2.getPointer();
	char * p3 = (char *)pmem3.getPointer();

	for (size_t i = 0 ; i < mem2.size(); i++)
	{
		check &= (p2[i] == p3[i]);
	}

	BOOST_REQUIRE_EQUAL(check,true);

	// Unpack on a Sparse grid without information

	Unpack_stat ps;
	sgrid empty(sz);

	empty.template unpack<0>(mem3,ps);

	BOOST_REQUIRE_EQUAL(empty.getGrid().size(0),grid.getGrid().size(0));
	BOOST_REQUIRE_EQUAL(empty.getGrid().size(1),grid.getGrid().size(1));
	BOOST_REQUIRE_EQUAL(empty.getGrid().size(2),grid.getGrid().size(2));

	auto it = empty.getIterator();

	while (it.isNext())
	{
		auto p = it.get();

		check &= (grid.template get<0>(p) == empty.template get<0>(p));

		++it;
	}

	BOOST_REQUIRE_EQUAL(check,true);
}


template<typename sgrid> void Test_unpack_and_check_full_noprp(sgrid & grid)
{
	grid_key_dx<3> end;
	grid_key_dx<3> zero({0,0,0});
	size_t req2 = 0;
	size_t req3 = 0;
	size_t sz[3];

	for (size_t i = 0 ; i < 3 ; i++)
	{
		end.set_d(i,grid.getGrid().size(i) - 1);
		sz[i] = 0;
	}

	grid.template packRequest(req2);
	auto sub_it = grid.getIterator(zero,end);
	grid.template packRequest<0,1>(sub_it,req3);

	BOOST_REQUIRE_EQUAL(req2,req3);

	Pack_stat sts2;
	Pack_stat sts3;

	// allocate the memory 2
	HeapMemory pmem2;
	pmem2.allocate(req2);
	ExtPreAlloc<HeapMemory> & mem2 = *(new ExtPreAlloc<HeapMemory>(req2,pmem2));
	mem2.incRef();

	// allocate the memory 3
	HeapMemory pmem3;
	pmem3.allocate(req3);
	ExtPreAlloc<HeapMemory> & mem3 = *(new ExtPreAlloc<HeapMemory>(req3,pmem3));
	mem3.incRef();

	grid.template pack(mem2,sts2);
	grid.template pack<0,1>(mem3,sub_it,sts3);

	BOOST_REQUIRE_EQUAL(mem2.size(),mem3.size());

	bool check = true;

	char * p2 = (char *)pmem2.getPointer();
	char * p3 = (char *)pmem3.getPointer();

	for (size_t i = 0 ; i < mem2.size(); i++)
	{
		check &= (p2[i] == p3[i]);
	}

	BOOST_REQUIRE_EQUAL(check,true);

	// Unpack on a Sparse grid without information

	Unpack_stat ps;
	sgrid empty(sz);

	empty.template unpack(mem2,ps);

	BOOST_REQUIRE_EQUAL(empty.getGrid().size(0),grid.getGrid().size(0));
	BOOST_REQUIRE_EQUAL(empty.getGrid().size(1),grid.getGrid().size(1));
	BOOST_REQUIRE_EQUAL(empty.getGrid().size(2),grid.getGrid().size(2));

	auto it = empty.getIterator();

	while (it.isNext())
	{
		auto p = it.get();

		check &= (grid.template get<0>(p) == empty.template get<0>(p));

		++it;
	}

	BOOST_REQUIRE_EQUAL(check,true);
}

template<typename sgrid> void Test_unpack_and_check(sgrid & grid, sgrid & grid2, size_t (& sz_cell)[3])
{
	grid.getBackgroundValue().template get<0>() = 0.0;

	grid_key_dx<3> start({0,0,0});
	grid_key_dx<3> stop({250,250,250});

	Box<3,size_t> bx({0,0,0},{250,250,250});

	auto sub_it = grid.getIterator(start,stop);
	size_t req = 0;
	grid.template packRequest<0>(sub_it,req);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(req,pmem));
	mem.incRef();

	Pack_stat sts;

	grid.template pack<0>(mem,sub_it,sts);

	// now we unpack on another grid

	int ctx = 0;
	Unpack_stat usts;
	grid2.template unpack<0>(mem,sub_it,usts,ctx);

	bool match = true;
	auto it = grid.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		if (bx.isInside(key.toPoint()) == true)
		{match &= grid.template get<0>(key) == grid2.template get<0>(key);}

		++it;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( sparse_grid_sub_grid_it_packing)
{
	size_t sz[3] = {501,501,501};
	size_t sz_cell[3] = {500,500,500};

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid(sz);
	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid2(sz);
	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid3(sz);

	CellDecomposer_sm<3, float, shift<3,float>> cdsm;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	cdsm.setDimensions(domain, sz_cell, 0);

	fill_sphere(grid,cdsm);

	Test_unpack_and_check(grid,grid2,sz_cell);

	grid_key_dx<3> start({251,251,251});
	grid_key_dx<3> stop({490,490,490});


	size_t cnt2 = 0;

	auto sub_it = grid.getIterator(start,stop);

	while (sub_it.isNext())
	{
		auto p = sub_it.get();

		grid3.template insert<0>(p) = grid.template get<0>(p);
		grid3.template insert<1>(p) = grid.template get<1>(p);

		cnt2++;

		++sub_it;
	}

	Test_unpack_and_check(grid,grid3,sz_cell);

	grid_key_dx<3> start2({251,251,251});
	grid_key_dx<3> stop2({490,490,490});

	size_t cnt = 0;

	auto sub_it2 = grid.getIterator(start2,stop2);

	bool match = true;

	while (sub_it2.isNext())
	{
		auto p = sub_it2.get();

		match &= grid3.template insert<0>(p) == grid.template get<0>(p);
		match &= grid3.template insert<1>(p) == grid.template get<1>(p);

		cnt++;

		++sub_it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
	BOOST_REQUIRE_EQUAL(cnt,cnt2);
}


BOOST_AUTO_TEST_CASE( sparse_grid_remove_area)
{
	size_t sz[3] = {501,501,501};

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid(sz);

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});


	Box<3,size_t> bx_create({100,100,100},{400,400,400});
	Box<3,size_t> bx_delete({150,150,150},{350,350,350});

	grid_sm<3,void> gs(sz);
	grid_key_dx_iterator_sub<3> sub(gs,bx_create.getKP1(),bx_create.getKP2());

	while (sub.isNext())
	{
		auto p = sub.get();

		grid.template insert<0>(p) = 1.0;

		++sub;
	}

	grid.remove(bx_delete);

	bool check = true;

	size_t cnt = 0;
	auto it2 = grid.getIterator(bx_create.getKP1(),bx_create.getKP2());
	while (it2.isNext())
	{
		auto p = it2.get();

		check &= bx_delete.isInside(p.toPoint()) == false;

		cnt++;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(check,true);

	BOOST_REQUIRE_EQUAL(cnt,bx_create.getVolumeKey() - bx_delete.getVolumeKey());
}

BOOST_AUTO_TEST_CASE( sparse_grid_copy_to)
{
	size_t sz[3] = {501,501,501};

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	size_t sz_g2[3] = {259,27,27};
	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid2(sz_g2);

	grid_key_dx_iterator<3> key_it(grid2.getGrid());
	auto gs = grid2.getGrid();

	while (key_it.isNext())
	{
		auto key = key_it.get();

		grid2.insert<0>(key) = gs.LinId(key);

		++key_it;
	}

	Box<3,size_t> bx_src({1,1,1},{255,23,23});
	Box<3,size_t> bx_dst({5,5,5},{259,27,27});

	grid.copy_to(grid2,bx_src,bx_dst);

	BOOST_REQUIRE(grid.size() != 0);
	BOOST_REQUIRE_EQUAL(grid.size(),bx_dst.getVolumeKey());

	bool match = true;
	auto it_check = grid.getIterator(bx_dst.getKP1(),bx_dst.getKP2());

	while (it_check.isNext())
	{
		auto key = it_check.get();

		grid_key_dx<3> key2 = key - bx_dst.getKP1();
		key2 += bx_src.getKP1();

		if (grid.template get<0>(key) != gs.LinId(key2))
		{match = false;}

		++it_check;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	bx_dst += Point<3,size_t>({0,40,40});
	grid.copy_to(grid2,bx_src,bx_dst);

	BOOST_REQUIRE(grid.size() != 0);
	BOOST_REQUIRE_EQUAL(grid.size(),2*bx_dst.getVolumeKey());

	match = true;
	auto it_check2 = grid.getIterator(bx_dst.getKP1(),bx_dst.getKP2());

	while (it_check2.isNext())
	{
		auto key = it_check2.get();

		grid_key_dx<3> key2 = key - bx_dst.getKP1();
		key2 += bx_src.getKP1();

		if (grid.template get<0>(key) != gs.LinId(key2))
		{match = false;}

		++it_check2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( sparse_pack_full )
{
	size_t sz[3] = {501,501,501};
	size_t sz_cell[3] = {500,500,500};

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid(sz);

	CellDecomposer_sm<3, float, shift<3,float>> cdsm;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	cdsm.setDimensions(domain, sz_cell, 0);

	fill_sphere(grid,cdsm);

	Test_unpack_and_check_full(grid);
}

BOOST_AUTO_TEST_CASE( sparse_pack_full_noprp )
{
	size_t sz[3] = {501,501,501};
	size_t sz_cell[3] = {500,500,500};

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid(sz);

	CellDecomposer_sm<3, float, shift<3,float>> cdsm;

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	cdsm.setDimensions(domain, sz_cell, 0);

	fill_sphere(grid,cdsm);

	Test_unpack_and_check_full_noprp(grid);
}

BOOST_AUTO_TEST_CASE( sparse_operator_equal )
{
	size_t sz[3] = {270,270,270};

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid(sz);

	grid.getBackgroundValue().template get<0>() = 0.0;

	sgrid_cpu<3,aggregate<double,int>,HeapMemory> grid2;

	grid_key_dx_iterator<3> key_it(grid.getGrid());
	auto gs = grid.getGrid();

	while (key_it.isNext())
	{
		auto key = key_it.get();

		grid.insert<0>(key) = gs.LinId(key);

		++key_it;
	}

	grid2 = grid;

	BOOST_REQUIRE(grid.size() == grid2.size());

	bool match = true;
	auto it_check = grid.getIterator();

	while (it_check.isNext())
	{
		auto key = it_check.get();

		if (grid.template get<0>(key) != grid2.template get<0>(key))
		{match = false;}

		++it_check;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_SUITE_END()

