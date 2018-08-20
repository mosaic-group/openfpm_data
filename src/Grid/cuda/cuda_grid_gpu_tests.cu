/*
 * cuda_gpu_func.cpp
 *
 *  Created on: Jun 3, 2018
 *      Author: i-bird
 */

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Grid/map_grid.hpp"
#include "Point_test.hpp"
#include "Grid/grid_util_test.hpp"
#include "cuda_grid_unit_tests_func.cuh"

BOOST_AUTO_TEST_SUITE( grid_gpu_func_test )


BOOST_AUTO_TEST_CASE (gpu_computation_func)
{
#ifdef CUDA_GPU

	size_t sz[3] = {64,64,64};
	grid_gpu<3, Point_test<float> > c3(sz);

	grid_key_dx<3> k1({1,1,1});
	grid_key_dx<3> k2({62,62,62});

	c3.setMemory();

	auto gcf = c3.getGPUIterator(k1,k2);

	BOOST_REQUIRE_EQUAL(gcf.thr.x,16ul);
	BOOST_REQUIRE_EQUAL(gcf.thr.y,8ul);
	BOOST_REQUIRE_EQUAL(gcf.thr.z,8ul);

	BOOST_REQUIRE_EQUAL(gcf.wthr.x,4ul);
	BOOST_REQUIRE_EQUAL(gcf.wthr.y,8ul);
	BOOST_REQUIRE_EQUAL(gcf.wthr.z,8ul);

	grid_key_dx<3> k3({50,50,50});
	grid_key_dx<3> k4({62,62,62});
	grid_key_dx<3> k5({60,61,62});

	auto gcf2 = c3.getGPUIterator(k3,k4);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,13ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,8ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,8ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,2ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,2ul);

	gcf2 = c3.getGPUIterator(k3,k4,511);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,8ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,8ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,4ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,2ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,2ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,4ul);

	gcf2 = c3.getGPUIterator(k3,k4,1);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,1ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,13ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,13ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,13ul);

	gcf2 = c3.getGPUIterator(k3,k5,32);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,4ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,4ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,2ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,3ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,3ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,7ul);

	gcf2 = c3.getGPUIterator(k3,k5,1);

	BOOST_REQUIRE_EQUAL(gcf2.thr.x,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.y,1ul);
	BOOST_REQUIRE_EQUAL(gcf2.thr.z,1ul);

	BOOST_REQUIRE_EQUAL(gcf2.wthr.x,11ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.y,12ul);
	BOOST_REQUIRE_EQUAL(gcf2.wthr.z,13ul);

#endif
}


BOOST_AUTO_TEST_CASE (gpu_computation)
{
	#ifdef CUDA_GPU

	{
	size_t sz[3] = {64,64,64};
	grid_gpu<3, Point_test<float> > c3(sz);

	c3.setMemory();
	test_layout_gridNd<3>(c3,sz[0]);

	gpu_grid_3D_compute(c3);

	c3.deviceToHost<0>();

	auto it = c3.getIterator();

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();

		good &= c3.getGrid().LinId(key) == c3.template get<0>(key);

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);

	}

	#endif
}

BOOST_AUTO_TEST_CASE (gpu_computation_stencil)
{
	#ifdef CUDA_GPU

	{
	size_t sz[3] = {64,64,64};
	grid_gpu<3, Point_test<float> > c3(sz);
	grid_gpu<3, Point_test<float> > c2(sz);
	grid_key_dx<3> key1({1,1,1});
	grid_key_dx<3> key2({62,62,62});


	c3.setMemory();
	c2.setMemory();
	test_layout_gridNd<3>(c3,sz[0]);
	test_layout_gridNd<3>(c2,sz[0]);

	gpu_grid_3D_one(c2);

	// Check property 1 is 1.0
	c2.deviceToHost<0>();

	{
	auto it = c2.getIterator();

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();

		good &= c2.get<0>(key) == 1.0;

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);
	}

	gpu_grid_3D_compute(c3);
	c3.deviceToHost<0>();

	{
	auto it = c3.getIterator();

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();

		good &= c3.getGrid().LinId(key) == c3.get<0>(key);

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);
	}

	gpu_grid_3D_compute_stencil(c3,c2,key1,key2);

	c2.deviceToHost<0>();

	auto it = c2.getIterator(key1,key2);

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();

		good &= c2.get<0>(key) == 0;

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);

	}

	#endif
}

BOOST_AUTO_TEST_CASE (gpu_computation_grid_stencil)
{
	#ifdef CUDA_GPU

	{
	size_t sz[3] = {64,64,64};
	grid_gpu<3, Point_test<float> > c3(sz);
	grid_gpu<3, Point_test<float> > c2(sz);
	grid_key_dx<3> key1({1,1,1});
	grid_key_dx<3> zero({0,0,0});
	grid_key_dx<3> key2({62,62,62});
	grid_key_dx<3> keyl({63,63,63});


	c3.setMemory();
	c2.setMemory();
	test_layout_gridNd<3>(c3,sz[0]);
	test_layout_gridNd<3>(c2,sz[0]);

	gpu_grid_3D_one(c2);

	// Check property 1 is 1.0
	c2.deviceToHost<0>();

	{
	auto it = c2.getIterator();

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();

		good &= c2.get<0>(key) == 1.0;

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);
	}

	gpu_grid_3D_compute(c3);
	c3.deviceToHost<0>();

	{
	auto it = c3.getIterator();

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();

		good &= c3.getGrid().LinId(key) == c3.get<0>(key);

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);
	}

	gpu_grid_3D_compute_grid_stencil(c3,c2,key1,key2);

	c2.deviceToHost<0>();

	auto it = c2.getIterator(key1,key2);

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();
		good &= c2.get<0>(key) == 0;

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);

	// We also try to fill a vectorial quantity

	gpu_grid_fill_vector(c3,zero,keyl);

	}

	#endif
}

BOOST_AUTO_TEST_CASE (gpu_computation_grid_stencil_vector)
{
	#ifdef CUDA_GPU

	{
	size_t sz[3] = {64,64,64};
	grid_gpu<3, Point_test<float> > c3(sz);
	grid_gpu<3, Point_test<float> > c2(sz);
	grid_key_dx<3> key1({1,1,1});
	grid_key_dx<3> zero({0,0,0});
	grid_key_dx<3> key2({62,62,62});
	grid_key_dx<3> keyl({63,63,63});


	c3.setMemory();
	c2.setMemory();

	gpu_grid_fill_vector(c3,zero,keyl);

	// Check property 1 is 1.0
	c3.deviceToHost<4>();

	{
	auto it = c3.getIterator(key1,key2);

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();

		good &= c3.get<4>(key)[0] == 1.0;
		good &= c3.get<4>(key)[1] == 2.0;
		good &= c3.get<4>(key)[2] == 3.0;

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);
	}

	// Fill c3

	gpu_grid_3D_compute(c3);
	gpu_grid_gradient_vector(c3,c2,key1,key2);

	// Check property 1 is 1.0
	c2.deviceToHost<4>();

	{
	auto it = c2.getIterator(key1,key2);

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();

		good &= c2.get<4>(key)[0] == 1.0;
		good &= c2.get<4>(key)[1] == 64.0;
		good &= c2.get<4>(key)[2] == 4096.0;

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);
	}

	}

	#endif
}

BOOST_AUTO_TEST_CASE (gpu_swap_vector)
{
	#ifdef CUDA_GPU

	{
	size_t sz[3] = {64,64,64};
	grid_gpu<3, Point_test<float> > c3(sz);
	grid_gpu<3, Point_test<float> > c2(sz);
	grid_key_dx<3> key1({1,1,1});
	grid_key_dx<3> zero({0,0,0});
	grid_key_dx<3> key2({62,62,62});
	grid_key_dx<3> keyl({63,63,63});


	c3.setMemory();
	c2.setMemory();

	gpu_grid_fill_vector(c2,zero,keyl);
	gpu_grid_fill_vector2(c3,zero,keyl);

	auto it4 = c3.getIterator(zero,keyl);

	// fill CPU
	while(it4.isNext())
	{
		auto key = it4.get();

		c2.get<4>(key)[0] = 1.0;
		c2.get<4>(key)[1] = 2.0;
		c2.get<4>(key)[2] = 3.0;

		c3.get<4>(key)[0] = 1001.0;
		c3.get<4>(key)[1] = 1002.0;
		c3.get<4>(key)[2] = 1003.0;

		++it4;
	}

	// now we swap

	// Check property 1 is 1.0
	c3.swap(c2);

	{
	auto it = c3.getIterator(zero,keyl);

	bool good = true;
	while(it.isNext())
	{
		auto key = it.get();

		good &= c3.get<4>(key)[0] == 1.0;
		good &= c3.get<4>(key)[1] == 2.0;
		good &= c3.get<4>(key)[2] == 3.0;

		good &= c2.get<4>(key)[0] == 1001.0;
		good &= c2.get<4>(key)[1] == 1002.0;
		good &= c2.get<4>(key)[2] == 1003.0;

		if (good == false)	{break;}

		// Set to zero

		c3.get<4>(key)[0] = 0.0;
		c3.get<4>(key)[1] = 0.0;
		c3.get<4>(key)[2] = 0.0;

		c2.get<4>(key)[0] = 0.0;
		c2.get<4>(key)[1] = 0.0;
		c2.get<4>(key)[2] = 0.0;

		++it;
	}

	BOOST_REQUIRE_EQUAL(good,true);

	c2.template deviceToHost<4>();
	c3.template deviceToHost<4>();

	auto it2 = c3.getIterator(zero,keyl);

	good = true;
	while(it2.isNext())
	{
		auto key = it2.get();

		good &= c3.get<4>(key)[0] == 1.0;
		good &= c3.get<4>(key)[1] == 2.0;
		good &= c3.get<4>(key)[2] == 3.0;

		good &= c2.get<4>(key)[0] == 1001.0;
		good &= c2.get<4>(key)[1] == 1002.0;
		good &= c2.get<4>(key)[2] == 1003.0;

		if (good == false)	{break;}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(good,true);
	}


	}

	#endif
}

BOOST_AUTO_TEST_CASE (gpu_copy_device)
{
	size_t sz[3] = {64,64,64};
	grid_gpu<3, Point_test<float> > c3(sz);

	grid_sm<3,void> g(sz);
	c3.setMemory();

	auto it4 = c3.getIterator();
	while (it4.isNext())
	{
		auto key = it4.get();

		c3.template get<0>(key) = g.LinId(key);

		++it4;
	}

	c3.hostToDevice<0>();

	size_t sz2[3] = {117,117,117};
	c3.resize(sz2);

	auto it = c3.getIterator();

	bool match = true;
	while (it.isNext())
	{
		auto key = it.get();

		if (key.get(0) < (unsigned int)sz[0] && key.get(1) < (unsigned int)sz[2] && key.get(2) < (unsigned int)sz[2])
		{
			match = c3.template get<0>(key) == g.LinId(key);
		}

		++it;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	// reset the memory

	auto it2 = c3.getIterator();

	match = true;
	while (it2.isNext())
	{
		auto key = it2.get();

		if (key.get(0) < (unsigned int)sz[0] && key.get(1) < (unsigned int)sz[2] && key.get(2) < (unsigned int)sz[2])
		{c3.template get<0>(key) = 0;}

		++it2;
	}

	// brint to CPU

	c3.deviceToHost<0>();

	auto it3 = c3.getIterator();

	match = true;
	while (it3.isNext())
	{
		auto key = it3.get();

		if (key.get(0) < (unsigned int)sz[0] && key.get(1) < (unsigned int)sz[2] && key.get(2) < (unsigned int)sz[2])
		{
			match = c3.template get<0>(key) == g.LinId(key);
		}

		++it3;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_SUITE_END()
