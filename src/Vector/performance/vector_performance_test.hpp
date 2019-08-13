/*
 * vector_performance_test.hpp
 *
 *  Created on: Jul 15, 2019
 *      Author: i-bird
 */

#ifndef VECTOR_PERFORMANCE_TEST_HPP_
#define VECTOR_PERFORMANCE_TEST_HPP_

/*
 * vector_performance_test.hpp
 *
 *  Created on: Jan 11, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_VECTOR_VECTOR_PERFORMANCE_TEST_HPP_
#define OPENFPM_DATA_SRC_VECTOR_VECTOR_PERFORMANCE_TEST_HPP_

#define NADD 128*128*128

// Property tree
struct report_vector_func_tests
{
	boost::property_tree::ptree graphs;
};

report_vector_func_tests report_vector_funcs;

BOOST_AUTO_TEST_SUITE( vector_performance )

BOOST_AUTO_TEST_CASE(vector_performance)
{
	report_vector_funcs.graphs.put("performance.vector(0).funcs.nele",NADD);
	report_vector_funcs.graphs.put("performance.vector(0).funcs.name","add");

	report_vector_funcs.graphs.put("performance.vector(1).funcs.nele",NADD);
	report_vector_funcs.graphs.put("performance.vector(1).funcs.name","get");

	std::vector<double> times(N_STAT + 1);
	std::vector<double> times_g(N_STAT + 1);

	// get test
	double tot_accu = 0.0;

	for (size_t i = 0 ; i < N_STAT+1 ; i++)
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

		// Add test

		for (size_t i = 0 ; i < NADD ; i++)
		{
			v1.add(p);
		}

		t.stop();
		times[i] = t.getwct();

		timer tg;
		tg.start();

		for (size_t i = 0 ; i < NADD ; i++)
		{
			double accu1 = v1.template get<P::x>(i);
			double accu2 = v1.template get<P::y>(i);
			double accu3 = v1.template get<P::z>(i);
			double accu4 = v1.template get<P::s>(i);

			double accu5 = v1.template get<P::v>(i)[0];
			double accu6 = v1.template get<P::v>(i)[1];
			double accu7 = v1.template get<P::v>(i)[2];

			double accu8 = v1.template get<P::t>(i)[0][0];
			double accu9 = v1.template get<P::t>(i)[0][1];
			double accu10 = v1.template get<P::t>(i)[0][2];
			double accu11 = v1.template get<P::t>(i)[1][0];
			double accu12 = v1.template get<P::t>(i)[1][1];
			double accu13 = v1.template get<P::t>(i)[1][2];
			double accu14 = v1.template get<P::t>(i)[2][0];
			double accu15 = v1.template get<P::t>(i)[2][1];
			double accu16 = v1.template get<P::t>(i)[2][2];

			tot_accu += accu1 + accu2 + accu3 + accu4 + accu5 + accu6 + accu7 + accu8 + accu9 + accu10 + accu11 + accu12 +
					   accu13 + accu14 + accu15 + accu16;
		}

		tg.stop();

		times_g[i] = tg.getwct();
	}

	std::cout << "Tot: " << tot_accu / 1000000.0 << std::endl;

	double mean;
	double dev;
	standard_deviation(times,mean,dev);

	report_vector_funcs.graphs.put("performance.vector(0).y.data.mean",mean);
	report_vector_funcs.graphs.put("performance.vector(0).y.data.dev",dev);

	standard_deviation(times_g,mean,dev);

	report_vector_funcs.graphs.put("performance.vector(1).y.data.mean",mean);
	report_vector_funcs.graphs.put("performance.vector(1).y.data.dev",dev);
}

/////// THIS IS NOT A TEST IT WRITE THE PERFORMANCE RESULT ///////

BOOST_AUTO_TEST_CASE(vector_performance_write_report)
{
	// Create a graphs

	report_vector_funcs.graphs.put("graphs.graph(0).type","line");
	report_vector_funcs.graphs.add("graphs.graph(0).title","Vector add and get");
	report_vector_funcs.graphs.add("graphs.graph(0).x.title","Tests");
	report_vector_funcs.graphs.add("graphs.graph(0).y.title","Time seconds");
	report_vector_funcs.graphs.add("graphs.graph(0).y.data(0).source","performance.vector(#).y.data.mean");
	report_vector_funcs.graphs.add("graphs.graph(0).x.data(0).source","performance.vector(#).funcs.name");
	report_vector_funcs.graphs.add("graphs.graph(0).y.data(0).title","Actual");
	report_vector_funcs.graphs.add("graphs.graph(0).interpolation","lines");

	boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
	boost::property_tree::write_xml("vector_performance_funcs.xml", report_vector_funcs.graphs,std::locale(),settings);

	GoogleChart cg;

	std::string file_xml_ref(test_dir);
	file_xml_ref += std::string("/openfpm_data/vector_performance_funcs_ref.xml");

	StandardXMLPerformanceGraph("vector_performance_funcs.xml",file_xml_ref,cg);

	addUpdtateTime(cg,1);

	cg.write("vector_performance_funcs.html");
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_VECTOR_VECTOR_PERFORMANCE_TEST_HPP_ */



#endif /* VECTOR_PERFORMANCE_TEST_HPP_ */
