/*
 * VerletList_test.hpp
 *
 *  Created on: Aug 16, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_TEST_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_TEST_HPP_

#include "NN/VerletList/VerletList.hpp"

/*! \brief create a vector of particles on a grid between 0.0 and 1.0
 *
 * \param box where the box is defined
 * \param v vector of positions
 *
 */
template<unsigned int dim, typename T> void create_particles_on_grid(grid_sm<dim,void> & g_info, SpaceBox<dim,T> & box, openfpm::vector<Point<dim,T>> & v)
{
	// Create a grid iterator
	grid_key_dx_iterator<dim> g_it(g_info);

	// Usefull definition of points
	Point<dim,T> end = box.getP2() - box.getP1();

	Point<dim,T> middle;
	Point<dim,T> spacing;

	for (size_t i = 0 ; i < dim ; i++)
	{
		middle.get(i) = end.get(i) / g_info.size(i) / 2.0;
		spacing.get(i) = end.get(i) / g_info.size(i);
	}

	while (g_it.isNext())
	{
		// Add 2 particles on each cell

		Point<dim,T> key = Point<dim,T>(g_it.get().toPoint());
		key = pmul(key,spacing) + box.getP1() + middle;

		v.add(key);

		++g_it;
	}
}

/*! \brief Test cell structure
 *
 * \tparam CellS
 *
 */
template<unsigned int dim, typename T, typename VerS> void Verlet_list_s(SpaceBox<dim,T> & box)
{
	T r_cut = 0.0753;

	for (size_t k = 0 ; k < 5; k++)
	{
		// id Cell list

		size_t div[dim];
		for (size_t i = 0 ; i < dim ; i++)
			div[i] = 20;

		grid_sm<3,void> ginfo(div);

		openfpm::vector<Point<dim,T>> pos;
		create_particles_on_grid(ginfo,box,pos);

		VerS vl1;

		Box<dim,T> innerBox;
		for (size_t i = 0 ; i < dim ; i++)
		{
			innerBox.setLow(i,r_cut);
			innerBox.setHigh(i,1.0-r_cut);
		}

		vl1.Initialize(box,box,r_cut,pos,pos.size());

		// Check that the verlet is consistent
		for (size_t i = 0 ; i < pos.size() ; i++)
		{
			if (innerBox.isInside(pos.get(i)) == false)
				continue;

			Point<dim,T> p = pos.get(i);

			if (k == 0)
			{BOOST_REQUIRE_EQUAL(vl1.getNNPart(i),19ul);}

			auto NN = vl1.getNNIterator(i);

			bool ret = true;

			while (NN.isNext())
			{
				auto k = NN.get();

				T dist = p.distance(Point<dim,T>(pos.get(k)));

				ret &= (dist < r_cut);

				if (ret == false)
				{
					std::cout << "Point1 " << p.toString() << "    " << "Point2 " << Point<dim,T>(pos.get(k)).toString() << "  " << r_cut << std::endl;
					exit(1);
				}

				++NN;
			}

			BOOST_REQUIRE_EQUAL(ret,true);
		}

		r_cut += 0.075;
	}

	r_cut = 0.0753;

	size_t div[dim];
	for (size_t i = 0 ; i < dim ; i++)
		div[i] = 20-1;

	grid_sm<3,void> ginfo(div);

	openfpm::vector<Point<dim,T>> pos;
	create_particles_on_grid(ginfo,box,pos);

	//! [Fill external cell list]

	// Initialize an external cell-list
	CellList<dim,T,Mem_fast<3,double>,shift<3,double>> cli;
	Box<dim,T> bt = box;

	// Calculate the divisions for the Cell-lists
	cl_param_calculate(bt,div,r_cut,Ghost<dim,T>(0.0));

	// Initialize a cell-list
	cli.Initialize(bt,div,5);

	for (size_t i = 0; i < pos.size() ; i++)
		cli.add(pos.get(i), i);

	//! [Fill external cell list]

	// Create cell-list

	for (size_t k = 0 ; k < 4; k++)
	{
		//! [create verlet cell]

		VerS vl1;
		vl1.Initialize(cli,r_cut,pos,pos.size());

		//! [create verlet cell]

		//! [create verlet]

		VerS vl2;
		vl2.Initialize(box,box,r_cut,pos,pos.size());

		//! [create verlet]

		Box<dim,T> innerBox;
		for (size_t i = 0 ; i < dim ; i++)
		{
			innerBox.setLow(i,r_cut);
			innerBox.setHigh(i,1.0-r_cut);
		}

		// Check that the verlet is consistent
		for (size_t i = 0 ; i < pos.size() ; i++)
		{
			if (innerBox.isInside(pos.get(i)) == false)
				continue;

			Point<dim,T> p = pos.get(i);

			if (k == 0)
			{BOOST_REQUIRE_EQUAL(vl1.getNNPart(i),19ul);}
			BOOST_REQUIRE_EQUAL(vl1.getNNPart(i),vl2.getNNPart(i));

			openfpm::vector<size_t> v1;
			openfpm::vector<size_t> v2;

			//! [usage of verlet]

			bool ret = true;
			auto NN = vl1.getNNIterator(i);
			while (NN.isNext())
			{
				auto k = NN.get();

				T dist = p.distance(Point<dim,T>(pos.get(k)));

				ret &= (dist < r_cut);

				v1.add(k);

				++NN;
			}

			//! [usage of verlet]

			auto NN2 = vl2.getNNIterator(i);
			while (NN2.isNext())
			{
				auto k = NN2.get();

				v2.add(k);

				++NN2;
			}

			ret = true;
			v1.sort();
			v2.sort();
			for (size_t i = 0 ; i < v1.size(); i++)
			{
				ret &= v1.get(i) == v2.get(i);
			}

			BOOST_REQUIRE_EQUAL(ret,true);
		}

		r_cut += 0.075;
	}
}



BOOST_AUTO_TEST_SUITE( VerletList_test )

BOOST_AUTO_TEST_CASE( VerletList_use)
{
	std::cout << "Test verlet list" << "\n";

	SpaceBox<3,double> box({0.0f,0.0f,0.0f},{1.0f,1.0f,1.0f});
	SpaceBox<3,double> box2({-1.0f,-1.0f,-1.0f},{1.0f,1.0f,1.0f});
	Verlet_list_s<3,double,VerletList<3,double,FAST,shift<3,double>>>(box);
//	Verlet_list_s<3,double,CellList<3,double,FAST,shift<3,double>> >(box2);
//	Verlet_list_sM<3,double,CellListM<3,double,8>>(box);
//	Verlet_list_sM<3,double,CellListM<3,double,8>>(box2);

	std::cout << "End verlet list" << "\n";

	// Test the cell list
}

BOOST_AUTO_TEST_SUITE_END()



#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_TEST_HPP_ */
