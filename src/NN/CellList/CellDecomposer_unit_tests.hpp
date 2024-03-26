/*
 * CellDecomposer_unit_tests.hpp
 *
 *  Created on: Feb 12, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLDECOMPOSER_UNIT_TESTS_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLDECOMPOSER_UNIT_TESTS_HPP_

BOOST_AUTO_TEST_SUITE( CellDecomposer_test )

BOOST_AUTO_TEST_CASE( CellDecomposer_get_grid_points )
{
	{
	// Cell decomposer
	CellDecomposer_sm<3,float> cd;

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Divisions
	size_t div[] = {10,10,10};

	// padding
	size_t padding = 1;

	// Set the dimensions of the decomposer
	cd.setDimensions(box,div,padding);

	Box<3,float> box_small({0.2,0.3,0.4},{0.5,0.5,0.6});

	// Get the grid points box
	Box<3,size_t> gp = cd.getGridPoints(box_small);

	BOOST_REQUIRE_EQUAL(gp.getLow(0),2+padding);
	BOOST_REQUIRE_EQUAL(gp.getLow(1),3+padding);
	BOOST_REQUIRE_EQUAL(gp.getLow(2),4+padding);

	BOOST_REQUIRE_EQUAL(gp.getHigh(0),5+padding-1);
	BOOST_REQUIRE_EQUAL(gp.getHigh(1),5+padding-1);
	BOOST_REQUIRE_EQUAL(gp.getHigh(2),6+padding-1);

	// Get the volume of the box
	size_t vol = gp.getVolumeKey();

	BOOST_REQUIRE_EQUAL(vol,12ul);
	}

	//////////////////////////
	// 2D test
	//////////////////////////

	{
	// Cell decomposer
	CellDecomposer_sm<2,float> cd;

	// Box
	Box<2,float> box({0.0,0.0},{1.0,1.0});

	// Divisions
	size_t div[] = {5,5};

	// padding
	size_t padding = 1;

	// Set the dimensions of the decomposer
	cd.setDimensions(box,div,padding);

	Box<2,float> box_small({0.4,0.4},{0.8,0.8});

	// Get the grid points box
	Box<2,size_t> gp = cd.getGridPoints(box_small);

	BOOST_REQUIRE_EQUAL(gp.getLow(0),2+padding);
	BOOST_REQUIRE_EQUAL(gp.getLow(1),2+padding);

	BOOST_REQUIRE_EQUAL(gp.getHigh(0),3+padding);
	BOOST_REQUIRE_EQUAL(gp.getHigh(1),3+padding);

	// Get the volume of the box
	size_t vol = gp.getVolumeKey();

	BOOST_REQUIRE_EQUAL(vol,4ul);
	}
}


BOOST_AUTO_TEST_CASE( CellDecomposer_use )
{
	//! [Cell decomposer use without shift]
	//Space where is living the Cell list
	Box<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Point<3,double> p({0.5,0.5,0.5});
	double pp[3] = {0.5,0.5,0.5};

	// Number of cell on each dimension
	size_t div[3] = {16,16,16};

	// grid info
	grid_sm<3,void> g_info(div);

	// Origin
	Point<3,double> org({0.0,0.0,0.0});

	// Test Cell decomposer
	{
	CellDecomposer_sm<3,double> cd(box,div,0);
	size_t cell = cd.getCell(p);
	BOOST_REQUIRE_EQUAL(cell,(size_t)(8*16*16 + 8*16 + 8));
	auto key = cd.getCellGrid(p);
	BOOST_REQUIRE_EQUAL(key.get(0),8);
	BOOST_REQUIRE_EQUAL(key.get(1),8);
	BOOST_REQUIRE_EQUAL(key.get(2),8);

	cell = cd.getCell(pp);
	BOOST_REQUIRE_EQUAL(cell,(size_t)(8*16*16 + 8*16 + 8));
	key = cd.getCellGrid(pp);
	BOOST_REQUIRE_EQUAL(key.get(0),8);
	BOOST_REQUIRE_EQUAL(key.get(1),8);
	BOOST_REQUIRE_EQUAL(key.get(2),8);
	}

	//! [Cell decomposer use without shift]

	//! [Test Cell decomposer with padding]
	{
	CellDecomposer_sm<3,double> cd(box,div,1);
	size_t cell = cd.getCell(p);
	BOOST_REQUIRE_EQUAL(cell,(size_t)(9*18*18 + 9*18 + 9));
	auto key = cd.getCellGrid(p);
	BOOST_REQUIRE_EQUAL(key.get(0),9);
	BOOST_REQUIRE_EQUAL(key.get(1),9);
	BOOST_REQUIRE_EQUAL(key.get(2),9);

	cell = cd.getCell(pp);
	BOOST_REQUIRE_EQUAL(cell,(size_t)(9*18*18 + 9*18 + 9));
	key = cd.getCellGrid(pp);
	BOOST_REQUIRE_EQUAL(key.get(0),9);
	BOOST_REQUIRE_EQUAL(key.get(1),9);
	BOOST_REQUIRE_EQUAL(key.get(2),9);
	}

	//! [Test Cell decomposer with padding]

	//! [Test Cell decomposer with shift]
	{
	Point<3,double> sht({1.0,2.0,3.0});
	Box<3,double> box2 = box;
	box2+= sht;
	double psht[3] = {1.5,2.5,3.5};

	CellDecomposer_sm< 3,double,shift<3,double> > cd(box2,div,1);
	size_t cell = cd.getCell(p + sht);
	BOOST_REQUIRE_EQUAL(cell,(size_t)(9*18*18 + 9*18 + 9));
	auto key = cd.getCellGrid(p + sht);
	BOOST_REQUIRE_EQUAL(key.get(0),9);
	BOOST_REQUIRE_EQUAL(key.get(1),9);
	BOOST_REQUIRE_EQUAL(key.get(2),9);

	cell = cd.getCell(psht);
	BOOST_REQUIRE_EQUAL(cell,(size_t)(9*18*18 + 9*18 + 9));
	key = cd.getCellGrid(psht);
	BOOST_REQUIRE_EQUAL(key.get(0),9);
	BOOST_REQUIRE_EQUAL(key.get(1),9);
	BOOST_REQUIRE_EQUAL(key.get(2),9);
	}

	//! [Test Cell decomposer with shift]

	//! [Test Cell decomposer getCellDom getCellPad]
	{
	CellDecomposer_sm<3,double> cd(box,div,1);
	size_t cell_cor_dom = cd.getCellDom(Point<3,double>({-0.000001,-0.000001,-0.000001}));
	size_t cell_cor_pad = cd.getCellPad(Point<3,double>({0.000001,0.000001,0.000001}));
	size_t cell_not_cor_pad1 = cd.getCell(Point<3,double>({-0.000001,-0.000001,-0.000001}));
	size_t cell_not_cor_pad2 = cd.getCell(Point<3,double>({0.000001,0.000001,0.000001}));
	BOOST_REQUIRE_EQUAL(cell_cor_pad,0ul);
	BOOST_REQUIRE_EQUAL(cell_cor_dom,(size_t)(1*18*18 + 1*18 + 1));
	BOOST_REQUIRE_EQUAL(cell_not_cor_pad1,0ul);
	BOOST_REQUIRE_EQUAL(cell_not_cor_pad2,(size_t)(1*18*18 + 1*18 + 1));
	}
}

#define N_POINTS 1000

BOOST_AUTO_TEST_CASE( CellDecomposer_consistent_use )
{
	{

	//! [Test Cell decomposer consistent construction]

	// number of divisions in each directions
	// random engine
	size_t div[3] = {16,16,16};
	std::default_random_engine g;
	std::uniform_real_distribution<double> d(0.0,1.0);

	Box<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Point<3,double> sht({1.1,2.1,3.1});
	Box<3,double> box2 = box;
	box2 += sht;

	CellDecomposer_sm< 3,double,shift<3,double> > cd(box2,div,1);

	// Here we create another extended Cell decomposer consistent with the
	// previous

	Box<3,size_t> ext({1,2,3},{1,1,1});
	CellDecomposer_sm< 3,double,shift<3,double> > cd2(cd,ext);

	//! [Test Cell decomposer consistent construction]

	// we create 2 grid_sm for the 2 CellDecomposer (used to check the consistency)
	size_t div_ext[3] = {6+16,6+16,6+16};
	// The old one has a padding 1 the new padding 3, result is padding 2
	grid_key_dx<3> key_base({2,2,2});

	grid_sm<3,void> cd_gr(div);
	grid_sm<3,void> cd2_gr(div_ext);

	// we randomly choose N_POINTS points inside box

	for (size_t i = 0 ; i < N_POINTS ; i++)
	{
		Point<3,double> p;
		p.get(0) = d(g);
		p.get(1) = d(g);
		p.get(2) = d(g);

		p += sht;

		grid_key_dx<3> key = cd.getCellGrid(p);
		grid_key_dx<3> key2 = cd2.getCellGrid(p);

		grid_key_dx<3> key_res = key2 - key_base;

		BOOST_REQUIRE(key == key_res);
	}

	grid_key_dx<3> key = cd.getCellGrid(sht);
	grid_key_dx<3> key2 = cd2.getCellGrid(sht);

	grid_key_dx<3> key_res = key2 - key_base;

	BOOST_REQUIRE(key == key_res);

	Point<3,double> p({1.0,1.0,1.0});
	p += sht;
	key = cd.getCellGrid(p);
	key2 = cd2.getCellGrid(p);

	key_res = key2 - key_base;

	BOOST_REQUIRE(key == key_res);

	}


}

BOOST_AUTO_TEST_CASE( CellDecomposer_swap_use )
{
	// number of divisions in each directions
	// random engine
	size_t div[3] = {16,16,16};

	Box<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	CellDecomposer_sm< 3,double,shift<3,double> > cd1(box,div,1);

	// another cell Decomposer

	// number of divisions in each directions
	// random engine
	size_t div2[3] = {15,15,15};

	Box<3,double> box2({-1.0f,-1.0f,-1.0f},{1.0f,1.0f,1.0f});

	CellDecomposer_sm< 3,double,shift<3,double> > cd2(box2,div2,2);

	auto cd1_old = cd1;
	auto cd2_old = cd2;

	cd1.swap(cd2);

	BOOST_REQUIRE(cd2 == cd1_old);
	BOOST_REQUIRE(cd1 == cd2_old);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLDECOMPOSER_UNIT_TESTS_HPP_ */
