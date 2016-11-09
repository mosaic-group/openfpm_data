/*
 * CellList_test.hpp
 *
 *  Created on: Mar 23, 2015
 *      Author: Pietro Incardona
 */

#include "CellList.hpp"
#include "CellListM.hpp"
#include "Grid/grid_sm.hpp"

#ifndef CELLLIST_TEST_HPP_
#define CELLLIST_TEST_HPP_


/*! \brief Test cell structure
 *
 * \tparam CellS
 *
 */
template<unsigned int dim, typename T, typename CellS> void Test_cell_s(SpaceBox<dim,T> & box)
{
	//! [Declare a cell list]
	//Space where is living the Cell list
	//SpaceBox<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

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
		// remove 1 particle on each cell

		//! [Usage of the neighborhood iterator]

		Point<dim,T> key = Point<dim,T>(g_it_s.get().toPoint());
		key = pmul(key,spacing) + offset[0] + box.getP1();

		auto NN = cl1.template getNNIterator<NO_CHECK>(cl1.getCell(key));
		size_t total = 0;

		while(NN.isNext())
		{
			// total

			total++;

			++NN;
		}

		//! [Usage of the neighborhood iterator]

		BOOST_REQUIRE_EQUAL(total,(size_t)openfpm::math::pow(3,dim));

		auto NNSym = cl1.template getNNIteratorSym<NO_CHECK>(cl1.getCell(key),id,pos);
		total = 0;

		while(NNSym.isNext())
		{
			// total

			total++;

			++NNSym;
		}

		BOOST_REQUIRE_EQUAL(total,(size_t)openfpm::math::pow(3,dim) / 2 + 1);

		++id;
		++g_it_s;
	}

}


/*! \brief Test cell structure
 *
 * \tparam CellS
 *
 */
template<unsigned int dim, typename T, typename CellS> void Test_cell_sM(SpaceBox<dim,T> & box)
{
	//Space where is living the Cell list
	//SpaceBox<dim,T> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});

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

	size_t id = 0;

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<dim,T> key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[0] + box.getP1();

		cl1.add(key,id,1);
		++id;

		key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + offset[1] + box.getP1();

		cl1.add(key,id,2);
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
		BOOST_REQUIRE_EQUAL((long int)(p1 - p2),1);
		BOOST_REQUIRE_EQUAL((long int)(v1 - v2),1);
		++g_it;
	}


	// Create a grid iterator
	grid_key_dx<dim> p1(1,1,1);
	grid_key_dx<dim> p2(div[0]-2,div[1]-2,div[2]-2);
	grid_key_dx_iterator_sub<dim> g_it_s(g_info,p1,p2);

	while (g_it_s.isNext())
	{
		Point<dim,T> key = Point<dim,T>(g_it_s.get().toPoint());
		key = pmul(key,spacing) + offset[0] + box.getP1();

		auto NN = cl1.template getNNIterator<NO_CHECK>(cl1.getCell(key));
		size_t total1 = 0;
		size_t total2 = 0;

		while(NN.isNext())
		{
			// total

			if (NN.getV() == 1)
				total1++;
			else
				total2++;

			++NN;
		}

		BOOST_REQUIRE_EQUAL(total1,(size_t)openfpm::math::pow(3,dim));
		BOOST_REQUIRE_EQUAL(total2,(size_t)openfpm::math::pow(3,dim));


		auto NNSym = cl1.template getNNIteratorSym<NO_CHECK>(cl1.getCell(key));
		total1 = 0;
		total2 = 0;

		while(NNSym.isNext())
		{
			// total

			if (NNSym.getV() == 1)
				total1++;
			else
				total2++;

			++NNSym;
		}

		BOOST_REQUIRE_EQUAL(total1,(size_t)openfpm::math::pow(3,dim) / 2 + 1);
		BOOST_REQUIRE_EQUAL(total2,(size_t)openfpm::math::pow(3,dim) / 2 + 1);

		++g_it_s;
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

BOOST_AUTO_TEST_SUITE( CellList_test )

BOOST_AUTO_TEST_CASE( CellList_use)
{
	std::cout << "Test cell list" << "\n";

	SpaceBox<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	SpaceBox<3,double> box2({-1.0f,-1.0f,-1.0f},{1.0f,1.0f,1.0f});
	Test_cell_s<3,double,CellList<3,double,Mem_fast<3,double>>>(box);
	Test_cell_s<3,double,CellList<3,double,Mem_fast<3,double>,shift<3,double>> >(box2);
	Test_cell_sM<3,double,CellListM<3,double,8>>(box);
	Test_cell_sM<3,double,CellListM<3,double,8>>(box2);


	Test_cell_s<3,double,CellList<3,double,Mem_bal<3,double>>>(box);
//	Test_cell_s<3,double,CellList<3,double,MEMORY>>();

	std::cout << "End cell list" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_CASE( CellList_consistent )
{
	Test_CellDecomposer_consistent<CellList<2,float,Mem_fast<3,double>,shift<2,float>>>();
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* CELLLIST_TEST_HPP_ */
