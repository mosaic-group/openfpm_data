/*
 * CellListIterator_test.hpp
 *
 *  Created on: May 7, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_TEST_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_TEST_HPP_

#include "NN/CellList/CellListIterator.hpp"
#include "NN/CellList/ParticleIt_Cells.hpp"
#include "NN/CellList/ParticleItCRS_Cells.hpp"

/*! \brief Fill the cell-list with particles in the box 0.0,1.0
 *
 * \param k Number of particles
 * \param NN Cell-list
 *
 */
template<unsigned int dim, typename CellList> void FillCellList(size_t k, CellList & NN)
{
	float pos[dim];

	//Fill with particles
	for (size_t i = 0; i < k; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			pos[j] = rand()/double(RAND_MAX);
		}
		NN.add(pos,i);
	}
}

#include "CellListFast_gen.hpp"

BOOST_AUTO_TEST_SUITE( celllist_gen_and_iterator_tests )

BOOST_AUTO_TEST_CASE( celllist_lin_and_iterator_test )
{
	///////// INPUT DATA //////////

	const size_t dim = 3;

	size_t div[dim] = {4,5,6};

	//Number of particles
	size_t k = 300;

	///////////////////////////////

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Initialize a cell list
	CellList_gen<dim,float,Process_keys_lin> NN;

	NN.Initialize(box,div,k*0.9,1);

	float pos[dim];

	//Fill with particles
	for (size_t i = 0; i < k*0.9; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			pos[j] = rand()/double(RAND_MAX);
		}
		NN.add(pos,i);
	}

	//Test the iterator
	auto it_cl = NN.getIterator();

	size_t count = 0;

	while (it_cl.isNext())
	{
		auto p_key = it_cl.get();

		BOOST_REQUIRE(p_key < NN.get_gm());

		count++;
		++it_cl;
	}

	BOOST_REQUIRE_EQUAL(count,NN.get_gm());
}

BOOST_AUTO_TEST_CASE( celllist_hilb_and_iterator_test )
{
	///////// INPUT DATA //////////

	const size_t dim = 3;

	size_t div[dim] = {4,5,6};

	//Number of particles
	size_t k = 300;

	///////////////////////////////

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Initialize a cell list
	CellList_gen<dim,float,Process_keys_hilb> NN;

	NN.Initialize(box,div,k*0.9,1);

	FillCellList<dim>((size_t)k*0.9,NN);


	//Test the iterator
	auto it_cl = NN.getIterator();

	size_t count = 0;

	while (it_cl.isNext())
	{
		auto p_key = it_cl.get();

		BOOST_REQUIRE(p_key < NN.get_gm());

		count++;
		++it_cl;
	}

	BOOST_REQUIRE_EQUAL(count,NN.get_gm());

	// Load previous results and check equality

	openfpm::vector<size_t> keys_old;

	keys_old.load("NN_hilb_keys");

	for (size_t i = 0; i < keys_old.size(); i++)
	{
		size_t a1 = keys_old.get(i);
		size_t a2 = NN.getCellSFC().getKeys().get(i);

		BOOST_REQUIRE_EQUAL(a1,a2);
	}

	size_t s1 = keys_old.size();
	size_t s2 = NN.getCellSFC().getKeys().size();

	BOOST_REQUIRE_EQUAL(s1,s2);
}

BOOST_AUTO_TEST_CASE( ParticleItCRS_Cells_iterator )
{
	///////// INPUT DATA //////////

	const size_t dim = 3;

	size_t div[dim] = {4,5,6};
	size_t div_p[dim] = {6,7,8};

	grid_sm<3,void> gs(div_p);

	//Number of particles
	size_t k = 300;

	///////////////////////////////

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Initialize a cell list
	CellList<dim,float,Mem_fast<dim,float>,shift<dim,float>> NN;

	NN.Initialize(box,div,1);

	FillCellList<dim>(k,NN);

	float pos[dim];

	grid_key_dx<3> start(0,0,0);
	grid_key_dx<3> stop(div[2]-1+2*NN.getPadding(2),div[1]-1+2*NN.getPadding(1),div[0]-1+2*NN.getPadding(0));

	openfpm::vector<size_t> dom;
	openfpm::vector<subsub_lin<dim>> anom;

	// No anomalous cells + all domain cells
	grid_key_dx_iterator_sub<dim> it(gs,start,stop);

	while (it.isNext())
	{
		auto key = it.get();

		dom.add(gs.LinId(key));

		++it;
	}

	//Test the iterator
	ParticleItCRS_Cells<dim,CellList<dim,float,Mem_fast<dim,float>,shift<dim,float>>> it_cl(NN,dom,anom,NN.getNNc_sym());

	size_t count = 0;

	while (it_cl.isNext())
	{
		count++;
		++it_cl;
	}

	BOOST_REQUIRE_EQUAL(count,k);

	//////////////////////////////// KEY ///////////////////////////////

	NN.clear();

	grid_key_dx<3> start2(1,1,1);
	grid_key_dx<3> stop2(div[2]-2+2*NN.getPadding(2),div[1]-2+2*NN.getPadding(1),div[0]-2+2*NN.getPadding(0));

	//Fill with particles
	for (size_t i = 0; i < k; i++)
	{
		pos[0] = 0.999*rand()/double(RAND_MAX) + 0.0001;
		pos[1] = 0.999*rand()/double(RAND_MAX) + 0.0001;
		pos[2] = 0.999*rand()/double(RAND_MAX) + 0.0001;

		NN.add(pos,i);
	}

	dom.clear();

	bool alternate = false;

	// No anomalous cells + all domain cells
	grid_key_dx_iterator_sub<dim> it2(gs,start,stop);

	while (it2.isNext())
	{
		auto key = it2.get();

		if (alternate == false)
		{
			dom.add(gs.LinId(key));
			alternate = true;
		}
		else
		{
			anom.add();
			anom.last().subsub = gs.LinId(key);
			alternate = false;
		}

		++it2;
	}

	ParticleItCRS_Cells<dim,CellList<dim,float,Mem_fast<dim,float>,shift<dim,float>>> it_cl2(NN,dom,anom,NN.getNNc_sym());

	count = 0;

	while (it_cl2.isNext())
	{
		count++;
		++it_cl2;
	}

	BOOST_REQUIRE_EQUAL(count,k);
}

BOOST_AUTO_TEST_CASE( ParticleIt_Cells_NN_iterator )
{
	///////// INPUT DATA //////////

	const size_t dim = 3;

	size_t div[dim] = {4,5,6};
	size_t div_p[dim] = {6,7,8};

	grid_sm<3,void> gs2(div);
	grid_sm<3,void> gs(div_p);

	///////////////////////////////

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Initialize a cell list
	CellList<dim,float,Mem_fast<dim,float>,shift<dim,float>> NN;

	NN.Initialize(box,div,1);

	grid_key_dx_iterator<3> it(gs2);

	float spacing[dim];
	float middle[dim];
	for (size_t i = 0 ; i < dim ; i++)
	{
		spacing[i] = box.getHigh(i) / div[i];
		middle[i] = spacing[i]/2;
	}

	size_t cnt = 0;

	//Fill with particles
	openfpm::vector<Point<dim,float>> vp;
	while (it.isNext())
	{
		auto key = it.get();

		Point<3,float> p;

		p.get(0) = key.get(0)*spacing[0] + middle[0];
		p.get(1) = key.get(1)*spacing[1] + middle[1];
		p.get(2) = key.get(2)*spacing[2] + middle[2];

		vp.add(p);

		NN.add(p,cnt);
		cnt++;

		++it;
	}

	grid_key_dx<3> start2(2,2,2);
	grid_key_dx<3> stop2(div[2]-3+2*NN.getPadding(2),div[1]-3+2*NN.getPadding(1),div[0]-3+2*NN.getPadding(0));

	openfpm::vector<size_t> dom;
	openfpm::vector<subsub_lin<dim>> anom;

	// No anomalous cells + all domain cells
	grid_key_dx_iterator_sub<dim> it2(gs,start2,stop2);

	bool alternate = false;
	while (it2.isNext())
	{
		auto key = it2.get();

		if (alternate == false)
		{
			dom.add(gs.LinId(key));
			alternate = true;
		}
		else
		{
			anom.add();
			anom.last().subsub = gs.LinId(key);

			for(size_t j = 0 ; j < openfpm::math::pow(3,dim)/2+1 ; j++)
				anom.last().NN_subsub.add(NN.getNNc_sym()[j]);

			alternate = false;
		}

		++it2;
	}

	//Test the iterator
	ParticleItCRS_Cells<dim,CellList<dim,float,Mem_fast<dim,float>,shift<dim,float>>> it_cl(NN,dom,anom,NN.getNNc_sym());

	size_t count = 0;

	while (it_cl.isNext())
	{
		auto NN_it = it_cl.getNNIteratorCSR(vp);

		size_t size_NN = 0;

		while (NN_it.isNext())
		{
			size_NN++;

			++NN_it;
		}

		BOOST_REQUIRE_EQUAL(size_NN,14ul);

		count++;
		++it_cl;
	}

	BOOST_REQUIRE_EQUAL(count,(div[0]-2)*(div[1]-2)*(div[2]-2));
}

BOOST_AUTO_TEST_CASE( ParticleIt_Cells_iterator )
{
	///////// INPUT DATA //////////

	const size_t dim = 3;

	size_t div[dim] = {4,5,6};
	size_t div_p[dim] = {6,7,8};

	grid_sm<3,void> gs(div_p);

	//Number of particles
	size_t k = 300;

	///////////////////////////////

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Initialize a cell list
	CellList<dim,float,Mem_fast<dim,float>,shift<dim,float>> NN;

	NN.Initialize(box,div,1);

	FillCellList<dim>(k,NN);

	grid_key_dx<3> start(0,0,0);
	grid_key_dx<3> stop(div[2]-1+2*NN.getPadding(2),div[1]-1+2*NN.getPadding(1),div[0]-1+2*NN.getPadding(0));


	// No anomalous cells + all domain cells
	grid_key_dx_iterator_sub<dim> it(gs,start,stop);

	openfpm::vector<size_t> dom;

	while (it.isNext())
	{
		auto key = it.get();

		dom.add(gs.LinId(key));

		++it;
	}

	//Test the iterator
	ParticleIt_Cells<dim,CellList<dim,float,Mem_fast<dim,float>,shift<dim,float>>> it_cl(NN,dom);

	size_t count = 0;

	while (it_cl.isNext())
	{
		count++;
		++it_cl;
	}

	BOOST_REQUIRE_EQUAL(count,k);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLISTITERATOR_TEST_HPP_ */
