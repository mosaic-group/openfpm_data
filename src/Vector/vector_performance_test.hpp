/*
 * vector_performance_test.hpp
 *
 *  Created on: Jan 11, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_VECTOR_VECTOR_PERFORMANCE_TEST_HPP_
#define OPENFPM_DATA_SRC_VECTOR_VECTOR_PERFORMANCE_TEST_HPP_

#define NADD 128*128*128

openfpm::vector<std::string> testsv;
openfpm::vector<float> per_timesv;

BOOST_AUTO_TEST_CASE(vector_add_performance)
{
	std::vector<double> times(N_STAT + 1);
	times[0] = 1000;

	for (size_t j = 0 ; j < 8 ; j++)
	{
		for (size_t i = 1 ; i < N_STAT+1 ; i++)
		{
			timer t;
			t.start();

			// create a vector
			openfpm::vector<Point_test<float>> v1;

			// Point
			Point_test<float> p;
			p.setx(1.0);
			p.sety(2.0);
			p.setz(3.0);
			p.sets(4.0);

			p.get<P::v>()[0] = 1.0;
			p.get<P::v>()[1] = 2.0;
			p.get<P::v>()[2] = 7.0;

			p.get<P::t>()[0][0] = 10.0;
			p.get<P::t>()[0][1] = 13.0;
			p.get<P::t>()[0][2] = 8.0;
			p.get<P::t>()[1][0] = 19.0;
			p.get<P::t>()[1][1] = 23.0;
			p.get<P::t>()[1][2] = 5.0;
			p.get<P::t>()[2][0] = 4.0;
			p.get<P::t>()[2][1] = 3.0;
			p.get<P::t>()[2][2] = 11.0;

			// push objects

			for (size_t i = 0 ; i < NADD ; i++)
			{
				v1.add(p);
			}

			t.stop();
			times[i] = t.getwct();
		}
		std::sort(times.begin(),times.end());

		sleep(5);
	}

	testsv.add("Vector add");
	per_timesv.add(times[0]);
}

/////// THIS IS NOT A TEST IT WRITE THE PERFORMANCE RESULT ///////

BOOST_AUTO_TEST_CASE(vector_performance_write_report)
{
	openfpm::vector<std::string> yn;
	openfpm::vector<openfpm::vector<float>> y;

	std::cout << "Write vector report\n";

	// Get the directory of the performance test files
	std::string per_dir(test_dir);

	load_and_combine(per_dir + std::string("/previous_measurev"),y,per_timesv);

	// Adding the dataset names
	if (y.size() != 0)
	{
		for (size_t j = 0; j < y.get(0).size(); j++)
			yn.add("config " + std::to_string(j));
	}

	// Google charts options
	GCoptions options;

	options.title = std::string("Vector Performances");
	options.yAxis = std::string("Time (seconds)");
	options.xAxis = std::string("Benchmark");
	options.stype = std::string("bars");

	cg.addHTML("<h2>Vector performance test</h2>");
	cg.AddColumsGraph(testsv,y,yn,options);
}

#endif /* OPENFPM_DATA_SRC_VECTOR_VECTOR_PERFORMANCE_TEST_HPP_ */
