/*
 * grid_performance_tests.hpp
 *
 *  Created on: Nov 1, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_PERFORMANCE_TESTS_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_PERFORMANCE_TESTS_HPP_

#include "grid_util_test.hpp"
#include "timer.hpp"

#define N_STAT 256
#define N_STAT_SMALL 32
#define N_TRY 8

#ifdef PERFORMANCE_TEST

BOOST_AUTO_TEST_SUITE( grid_performance )

BOOST_AUTO_TEST_CASE(grid_performance_set_obj)
{
	size_t sz[] = {128,128,128};

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);

	Point_test<float> f;
	f.fill();

	std::vector<double> times(N_STAT + 1);
	times[0] = 1000;

	for (size_t j = 0 ; j < 8 ; j++)
	{
		for (size_t i = 1 ; i < N_STAT+1 ; i++)
		{
			timer t;
			t.start();

			auto it = c3.getIterator();

			while (it.isNext())
			{
				c3.set(it.get(),f);

				++it;
			}

			t.stop();

			times[i] = t.getcputime();
		}
		std::sort(times.begin(),times.end());
		sleep(5);
	}

	std::cout << "Time : " <<  times[0] << "s ";

}

BOOST_AUTO_TEST_CASE(grid_performance_set_other_grid)
{
	size_t sz[] = {128,128,128};

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);
	grid_cpu<3, Point_test<float> > c1 = c3.duplicate();

	std::vector<double> times(N_STAT + 1);
	times[0] = 1000;

	for (size_t j = 0 ; j < 8 ; j++)
	{
		for (size_t i = 1 ; i < N_STAT+1 ; i++)
		{
			timer t;
			t.start();

			auto it = c3.getIterator();

			while (it.isNext())
			{
				c3.set(it.get(),c1,it.get());

				++it;
			}

			t.stop();

			times[i] = t.getcputime();
		}
		std::sort(times.begin(),times.end());
		sleep(5);
	}

	std::cout << "Time : " <<  times[0] << "s ";

}

BOOST_AUTO_TEST_CASE(grid_performance_set_other_grid_encap)
{
	size_t sz[] = {128,128,128};

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);
	grid_cpu<3, Point_test<float> > c1 = c3.duplicate();

	std::vector<double> times(N_STAT + 1);
	times[0] = 1000;

	for (size_t j = 0 ; j < 8 ; j++)
	{
		for (size_t i = 1 ; i < N_STAT+1 ; i++)
		{
			timer t;
			t.start();

			auto it = c3.getIterator();

			while (it.isNext())
			{
				c3.set(it.get(),c1.get_o(it.get()));

				++it;
			}

			t.stop();

			times[i] = t.getcputime();
		}
		std::sort(times.begin(),times.end());
		sleep(5);
	}

	std::cout << "Time : " <<  times[0] << "s ";
}

BOOST_AUTO_TEST_CASE(grid_performance_duplicate)
{
	size_t sz[] = {128,128,128};

	grid_cpu<3, Point_test<float> > c3(sz);
	c3.setMemory();

	fill_grid<3>(c3);
	grid_cpu<3, Point_test<float> > c1;

	std::vector<double> times(N_STAT_SMALL + 1);
	times[0] = 1000;

	for (size_t j = 0 ; j < 8 ; j++)
	{
		for (size_t i = 1 ; i < N_STAT_SMALL+1 ; i++)
		{
			timer t;
			t.start();

			c1 = c3.duplicate();

			t.stop();

			times[i] = t.getcputime();
		}
		std::sort(times.begin(),times.end());
		sleep(5);
	}

	std::cout << "Time : " <<  times[0] << "s ";
}

BOOST_AUTO_TEST_SUITE_END()

#endif

#endif /* OPENFPM_DATA_SRC_GRID_GRID_PERFORMANCE_TESTS_HPP_ */
