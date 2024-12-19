/*
 * CellList_test.hpp
 *
 *  Created on: Mar 23, 2015
 *      Author: Pietro Incardona
 */

#include "NN/CellList/CellList.hpp"
#include "NN/CellList/multiphase/CellListM.hpp"
#include "Grid/grid_sm.hpp"

#ifndef CELLLIST_TEST_HPP_
#define CELLLIST_TEST_HPP_


/*! \brief Test cell structure
 *
 * \tparam CellS
 *
 */
template<unsigned int dim, typename T, typename CellS> void Test_cell_s(Box<dim,T> & box)
{
	//! [Declare a cell list]
	//Space where is living the Cell list
	//Box<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Subdivisions
	size_t div[dim] = {16,16,16};

	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// id Cell list
	CellS cl2(box,div);
	//! [Declare a cell list]

	// grid info
	grid_sm<dim,void> g_info(div);

	// Test force reallocation in case of Cell list fast
	for (size_t i = 0 ; i < CELL_REALLOC * 3 ; i++)
	{
		cl2.add(org,i);
	}

	// Check the elements
	BOOST_REQUIRE_EQUAL(cl2.getNelements(cl2.getCell(org)),CELL_REALLOC * 3ul);
	for (size_t i = 0 ; i < CELL_REALLOC * 3 ; i++)
	{
		BOOST_REQUIRE_EQUAL(cl2.get(cl2.getCell(org),i),i);
	}

	//! [Usage of cell list]

	// id Cell list
	CellS cl1(box,div);

	// Create a grid iterator
	grid_key_dx_iterator<dim> g_it(g_info);

	// Iterate through each element
	// Add 1 element for each cell

	// Usefull definition of points
	Point<dim,T> end = box.getP2() - box.getP1();
	Point<dim,T> middle = end / div / 2.0;
	Point<dim,T> spacing = end / div;

	Point<dim,T> offset[dim] = {middle,middle,middle};

	// Create offset shift vectors
	for (size_t i = 0 ; i < dim ; i++)
	{
		offset[i].get(i) += (1.0 / div[i]) / 8.0;
	}

	openfpm::vector<Point<dim,T>> pos;
	size_t id = 0;

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<dim,T> key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[0] + box.getP1();
		pos.add(key);

		cl1.add(key,id);
		++id;

		key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[1] + box.getP1();
		pos.add(key);

		cl1.add(key,id);
		++id;

		++g_it;
	}

	//! [Usage of cell list]

	// check the cell are correctly filled

	// reset iterator
	g_it.reset();

	while (g_it.isNext())
	{
		// Check that there are 2 particles on each cell

		Point<dim,T> key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[2] + box.getP1();

		size_t cell = cl1.getCell(key);
		size_t n_ele = cl1.getNelements(cell);

		BOOST_REQUIRE_EQUAL(n_ele,2ul);
		BOOST_REQUIRE_EQUAL((long int)(cl1.get(cell,1) - cl1.get(cell,0)),1);

		++g_it;
	}

	// reset itarator
	g_it.reset();

	//! [remove one particle from each cell]

	while (g_it.isNext())
	{
		// remove 1 particle on each cell

		Point<dim,T> key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[0] + box.getP1();

		auto cell = cl1.getCell(key);

		// Remove the first particle in the cell
		cl1.remove(cell,0);
		++g_it;
	}

	//! [remove one particle from each cell]

	// Check we have 1 object per cell
	g_it.reset();

	while (g_it.isNext())
	{
		// remove 1 particle on each cell

		Point<dim,T> key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[0] + box.getP1();

		auto cell = cl1.getCell(key);
		size_t n_ele = cl1.getNelements(cell);

		BOOST_REQUIRE_EQUAL(n_ele,1ul);
		++g_it;
	}


	// Check we have 1 object per cell

	// Create a grid iterator
	grid_key_dx<dim> p1(1,1,1);
	grid_key_dx<dim> p2(div[0]-2,div[1]-2,div[2]-2);
	grid_key_dx_iterator_sub<dim> g_it_s(g_info,p1,p2);
	id = 0;

	while (g_it_s.isNext())
	{
		//! [Usage of the neighborhood iterator]

		Point<dim,T> key = Point<dim,T>(g_it_s.get().toPoint());
		key = pmul(key,spacing) + offset[0] + box.getP1();

		auto NN = cl1.getNNIteratorBox(cl1.getCell(key));
		size_t total = 0;

		while(NN.isNext())
		{
			// total

			total++;

			++NN;
		}

		//! [Usage of the neighborhood iterator]

		BOOST_REQUIRE_EQUAL(total,(size_t)openfpm::math::pow(3,dim));

		// in SE1_CLASS the cell list consider this construction as an hack
		// disable the test

#ifndef SE_CLASS1

		id = cl1.get(cl1.getCell(key),0);
		auto NNSym = cl1.getNNIteratorBoxSym(cl1.getCell(key),id,pos);
		total = 0;

		while(NNSym.isNext())
		{
			// total

			total++;

			++NNSym;
		}

		BOOST_REQUIRE_EQUAL(total,(size_t)openfpm::math::pow(3,dim) / 2 + 1);

#endif

		++g_it_s;
	}

}


/*! \brief Test cell structure
 *
 * \tparam CellS
 *
 */
template<unsigned int dim, typename T, typename CellS> void Test_cell_sM(Box<dim,T> & box)
{
	//Space where is living the Cell list
	//Box<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Subdivisions
	size_t div[dim] = {16,16,16};

	// Origin
	Point<dim,T> org({0.0,0.0,0.0});

	// grid info
	grid_sm<dim,void> g_info(div);

	//! [Usage of cell list multi]

	// CellS = CellListM<dim,T,8>
	CellS cl1(box,div);

	// Create a grid iterator
	grid_key_dx_iterator<dim> g_it(g_info);

	// Iterate through each element
	// Add 1 element for each cell

	// Usefull definition of points
	Point<dim,T> end = box.getP2() - box.getP1();
	Point<dim,T> middle = end / div / 2.0;
	Point<dim,T> spacing = end / div;

	Point<dim,T> offset[dim] = {middle,middle,middle};

	// Create offset shift vectors
	for (size_t i = 0 ; i < dim ; i++)
	{
		offset[i].get(i) += (1.0 / div[i]) / 8.0;
	}

	openfpm::vector<Point<dim,T>> phase1;
	openfpm::vector<Point<dim,T>> phase2;

	openfpm::vector<pos_v<openfpm::vector<Point<dim,T>>>> phases;
	phases.add(pos_v<openfpm::vector<Point<dim,T>>>(phase1));
	phases.add(pos_v<openfpm::vector<Point<dim,T>>>(phase2));

	size_t id = 0;

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<dim,T> key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[0] + box.getP1();

		phase1.add(key);

		cl1.add(key,id,0);

		key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[1] + box.getP1();

		phase2.add(key);
		cl1.add(key,id,1);
		++id;

		++g_it;
	}

	// check the cell are correctly filled

	// reset iterator
	g_it.reset();

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<dim,T> key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[2] + box.getP1();

		size_t cell = cl1.getCell(key);
		size_t n_ele = cl1.getNelements(cell);

		size_t p1 = cl1.getP(cell,1);
		size_t p2 = cl1.getP(cell,0);

		size_t v1 = cl1.getV(cell,1);
		size_t v2 = cl1.getV(cell,0);

		BOOST_REQUIRE_EQUAL(n_ele,2ul);
		BOOST_REQUIRE_EQUAL((long int)(p1 - p2),0);
		BOOST_REQUIRE_EQUAL((long int)(v1 - v2),1);
		++g_it;
	}
}



template<typename CellList> void Test_CellDecomposer_consistent()
{
	Box<2,float> bx({-1.0/3.0,-1.0/3.0},{1.0/3.0,1.0/3.0});

	size_t div[2] = {36,36};

	CellDecomposer_sm<2,float,shift<2,float>> cd(bx,div,1);

	Box<2,float> bx_sub({-1.0/5.0,-1.0/5.0},{1.0/5.0,1.0/5.0});

	size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
	CellList cl(cd,bx_sub);
	Box<2,long int> bx_int = cd.convertDomainSpaceIntoGridUnits(bx_sub,bc);

	BOOST_REQUIRE_EQUAL(bx_int.getLow(0),8);
	BOOST_REQUIRE_EQUAL(bx_int.getLow(1),8);

	BOOST_REQUIRE_EQUAL(bx_int.getHigh(0),28);
	BOOST_REQUIRE_EQUAL(bx_int.getHigh(1),28);

	cd.convertCellUnitsIntoDomainSpace(bx_sub);

	BOOST_REQUIRE_EQUAL(bx_sub.getLow(0),-1.0f/5.0f);
	BOOST_REQUIRE_EQUAL(bx_sub.getLow(1),-1.0f/5.0f);

	BOOST_REQUIRE_EQUAL(bx_sub.getHigh(0),1.0f/5.0f);
	BOOST_REQUIRE_EQUAL(bx_sub.getHigh(1),1.0f/5.0f);
}

/*! \brief Test cell structure
 *
 * \tparam CellS
 *
 */
template<unsigned int dim, typename T, typename CellS> void Test_NN_iterator_radius(Box<dim,T> & box)
{
	//Space where is living the Cell list
	//Box<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

	// Subdivisions
	size_t div[dim];
	size_t div2[dim];

	for (size_t i = 0 ; i < dim ; i++)
	{
		div[i] = 17;
		div2[i] = 34;
	}

	// grid info
	grid_sm<dim,void> g_info(div);
	grid_sm<dim,void> g_info2(div2);

	//! [Usage of cell list multi]

	// CellS = CellListM<dim,T,8>
	CellS cl1(box,div,2);
	CellS cl2(box,div2,3);

	T radius = (box.getHigh(0) - box.getLow(0))/div[0];

	cl2.setRadius( radius );

	// create a vector of random point

	openfpm::vector<Point<dim,float>> vrp;

	for (size_t j = 0 ; j < 10000 ; j++)
	{
		vrp.add();

		for (size_t i = 0 ; i < dim ; i++)
		{
			vrp.template get<0>(j)[i] = ((float)rand() / (float)RAND_MAX)*(box.getHigh(i) - box.getLow(i)) + box.getLow(i);
		}
	}

	auto g_it = vrp.getIterator();

	while (g_it.isNext())
	{
		Point<dim,T> xp = vrp.get(g_it.get());

		size_t debug = cl1.getCell(xp);

		cl1.add(xp,g_it.get());
		cl2.add(xp,g_it.get());

		++g_it;
	}

	// Get the neighborhood of each particle and compare the numbers of cell

	bool match = true;
	auto g_it2 = vrp.getIterator();

	size_t number_of_nn = 0;
	size_t number_of_nn2 = 0;

	while (g_it2.isNext())
	{
		Point<dim,T> xp = vrp.get(g_it2.get());

		openfpm::vector<size_t> ids1;
		openfpm::vector<size_t> ids2;

		size_t local1 = 0;
		size_t local2 = 0;

		auto NNit = cl1.getNNIteratorBox(cl1.getCell(xp));

		while (NNit.isNext())
		{
			auto q = NNit.get();

			// calculate distance

			Point<dim,T> xq = vrp.get(q);
			Point<dim,T> r = xq - xp;

			if (r.norm() <= radius)
			{ids1.add(q);}

			number_of_nn++;
			local1++;

			++NNit;
		}

		auto NN2it = cl2.getNNIteratorRadius(cl2.getCell(xp));

		while (NN2it.isNext())
		{
			auto q = NN2it.get();

			Point<dim,T> xq = vrp.get(q);
			Point<dim,T> r = xq - xp;

			if (r.norm() <= radius)
			{ids2.add(q);}

			number_of_nn2++;
			local2++;

			++NN2it;
		}

		// Sort ids1
		ids1.sort();
		ids2.sort();

		match &= ids1.size() == ids2.size();

		for (size_t i = 0 ; i < ids1.size() ; i++)
		{match &= ids1.get(i) == ids2.get(i);}

		if (match == false)
		{break;}

		++g_it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
	BOOST_REQUIRE(number_of_nn2 < number_of_nn);
}

BOOST_AUTO_TEST_SUITE( CellList_test )

BOOST_AUTO_TEST_CASE ( NN_radius_check )
{
	Box<2,float> box1({0.1,0.1},{0.3,0.5});
	Box<3,float> box2({0.0,0.1,0.2},{0.3,0.7,0.5});
	Box<3,float> box3({0.0,0.0,0.0},{1.0,1.0,1.0});

	std::cout << "Test cell list radius" << "\n";

	Test_NN_iterator_radius<2,float,CellList<2,float,Mem_fast<>,shift<2,float>>>(box1);
	Test_NN_iterator_radius<3,float,CellList<3,float,Mem_fast<>,shift<3,float>>>(box2);
	Test_NN_iterator_radius<3,float,CellList<3,float,Mem_fast<>,shift<3,float>>>(box3);

	std::cout << "End cell list" << "\n";
}

BOOST_AUTO_TEST_CASE( CellList_use)
{
	std::cout << "Test cell list" << "\n";

	Box<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	Box<3,double> box2({-1.0f,-1.0f,-1.0f},{1.0f,1.0f,1.0f});
	Test_cell_s<3,double,CellList<3,double,Mem_fast<>>>(box);
	Test_cell_s<3,double,CellList<3,double,Mem_fast<>,shift<3,double>> >(box2);
	Test_cell_sM<3,double,CellListM<3,double,8>>(box);
	Test_cell_sM<3,double,CellListM<3,double,8>>(box2);


	Test_cell_s<3,double,CellList<3,double,Mem_bal<>>>(box);
	Test_cell_s<3,double,CellList<3,double,Mem_mw<>>>(box);

	std::cout << "End cell list" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( CellList_consistent )
{
	Test_CellDecomposer_consistent<CellList<2,float,Mem_fast<>,shift<2,float>>>();
}

BOOST_AUTO_TEST_CASE( CellList_NNc_csr_calc )
{
	openfpm::vector<std::pair<grid_key_dx<3>,grid_key_dx<3>>> cNN;

	NNcalc_csr(cNN);

	BOOST_REQUIRE_EQUAL(cNN.size(),14ul);

	BOOST_REQUIRE(cNN.get(0).first == grid_key_dx<3>(0,0,0));
	BOOST_REQUIRE(cNN.get(0).second == grid_key_dx<3>(0,0,0));

	BOOST_REQUIRE(cNN.get(1).first == grid_key_dx<3>(0,0,0));
	BOOST_REQUIRE(cNN.get(1).second == grid_key_dx<3>(0,0,1));

	BOOST_REQUIRE(cNN.get(2).first == grid_key_dx<3>(0,0,1));
	BOOST_REQUIRE(cNN.get(2).second == grid_key_dx<3>(0,1,0));

	BOOST_REQUIRE(cNN.get(3).first == grid_key_dx<3>(0,0,0));
	BOOST_REQUIRE(cNN.get(3).second == grid_key_dx<3>(0,1,0));

	BOOST_REQUIRE(cNN.get(4).first == grid_key_dx<3>(0,0,0));
	BOOST_REQUIRE(cNN.get(4).second == grid_key_dx<3>(0,1,1));

	BOOST_REQUIRE(cNN.get(5).first == grid_key_dx<3>(0,1,1));
	BOOST_REQUIRE(cNN.get(5).second == grid_key_dx<3>(1,0,0));

	BOOST_REQUIRE(cNN.get(6).first == grid_key_dx<3>(0,1,0));
	BOOST_REQUIRE(cNN.get(6).second == grid_key_dx<3>(1,0,0));

	BOOST_REQUIRE(cNN.get(7).first == grid_key_dx<3>(0,1,0));
	BOOST_REQUIRE(cNN.get(7).second == grid_key_dx<3>(1,0,1));

	BOOST_REQUIRE(cNN.get(8).first == grid_key_dx<3>(0,0,1));
	BOOST_REQUIRE(cNN.get(8).second == grid_key_dx<3>(1,0,0));

	BOOST_REQUIRE(cNN.get(9).first == grid_key_dx<3>(0,0,0));
	BOOST_REQUIRE(cNN.get(9).second == grid_key_dx<3>(1,0,0));

	BOOST_REQUIRE(cNN.get(10).first == grid_key_dx<3>(0,0,0));
	BOOST_REQUIRE(cNN.get(10).second == grid_key_dx<3>(1,0,1));

	BOOST_REQUIRE(cNN.get(11).first == grid_key_dx<3>(0,0,1));
	BOOST_REQUIRE(cNN.get(11).second == grid_key_dx<3>(1,1,0));

	BOOST_REQUIRE(cNN.get(12).first == grid_key_dx<3>(0,0,0));
	BOOST_REQUIRE(cNN.get(12).second == grid_key_dx<3>(1,1,0));

	BOOST_REQUIRE(cNN.get(13).first == grid_key_dx<3>(0,0,0));
	BOOST_REQUIRE(cNN.get(13).second == grid_key_dx<3>(1,1,1));

}



BOOST_AUTO_TEST_SUITE_END()

#endif /* CELLLIST_TEST_HPP_ */
