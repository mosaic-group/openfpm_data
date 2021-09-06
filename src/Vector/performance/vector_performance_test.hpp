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

	double mean;
	double dev;
	standard_deviation(times,mean,dev);

	report_vector_funcs.graphs.put("performance.vector(0).y.data.mean",mean);
	report_vector_funcs.graphs.put("performance.vector(0).y.data.dev",dev);

	standard_deviation(times_g,mean,dev);

	report_vector_funcs.graphs.put("performance.vector(1).y.data.mean",mean);
	report_vector_funcs.graphs.put("performance.vector(1).y.data.dev",dev);
}

template<typename vector_prop_type, typename vector_pos_type>
__device__ __host__ void read_write(vector_prop_type & vd_prop, vector_pos_type & vd_pos, unsigned int p)
{
	vd_prop.template get<0>(p) = vd_pos.template get<0>(p)[0] + vd_pos.template get<0>(p)[1];

	vd_prop.template get<1>(p)[0] = vd_pos.template get<0>(p)[0];
	vd_prop.template get<1>(p)[1] = vd_pos.template get<0>(p)[1];

	vd_prop.template get<2>(p)[0][0] = vd_pos.template get<0>(p)[0];
	vd_prop.template get<2>(p)[0][1] = vd_pos.template get<0>(p)[1];
	vd_prop.template get<2>(p)[1][0] = vd_pos.template get<0>(p)[0] + 
                                      	   vd_pos.template get<0>(p)[1];
	vd_prop.template get<2>(p)[1][1] = vd_pos.template get<0>(p)[1] - 
                                      	   vd_pos.template get<0>(p)[0];

	vd_pos.template get<0>(p)[0] += 0.01f;
	vd_pos.template get<0>(p)[1] += 0.01f;
}


struct ele
{
	double s;
	double v[2];
	double t[2][2];
};

__device__ __host__ void read_write_lin(double * pos, ele * prp, unsigned int p)
{
    prp[p].s = pos[2*p] + pos[2*p+1];

    prp[p].v[0] = pos[2*p];
    prp[p].v[1] = pos[2*p+1];

    prp[p].t[0][0] = pos[2*p];
    prp[p].t[0][1] = pos[2*p+1];
    prp[p].t[1][0] = pos[2*p] + 
                                  pos[2*p+1];
    prp[p].t[1][1] = pos[2*p+1] - 
                                  pos[2*p];

    pos[2*p] += 0.01f;
    pos[2*p+1] += 0.01f;
}

__device__ __host__ void read_write_inte(double * pos, double * prp0, double * prp1, double * prp2, unsigned int p, unsigned int n_pos)
{
    prp0[0*n_pos + p] = pos[0*n_pos + p] + pos[1*n_pos+p];

    prp1[0*n_pos + p] = pos[0*n_pos + p];
    prp1[1*n_pos + p] = pos[1*n_pos + p];

    prp2[0*n_pos*2+0*n_pos + p] = pos[0*n_pos + p];
    prp2[0*n_pos*2+1*n_pos + p] = pos[1*n_pos + p];
    prp2[1*n_pos*2+0*n_pos + p] = pos[0*n_pos + p] + 
                                  pos[1*n_pos + p];
    prp2[1*n_pos*2+1*n_pos + p] = pos[1*n_pos + p] - 
                                  pos[0*n_pos + p];

    pos[0*n_pos + p] += 0.01f;
    pos[1*n_pos + p] += 0.01f;
}

BOOST_AUTO_TEST_CASE(vector_performance_layout_vs_plain_array)
{
	std::vector<double> times(N_STAT + 1);
	std::vector<double> times_g(N_STAT + 1);

	std::vector<double> times2(N_STAT + 1);
	std::vector<double> times2_g(N_STAT + 1);

	// get test
	double tot_accu = 0.0;

	report_vector_funcs.graphs.put("performance.vector_layout(0).funcs.nele",NADD);
	report_vector_funcs.graphs.put("performance.vector_layout(0).funcs.name","read_write_lin");

	for (size_t i = 0 ; i < N_STAT+1 ; i++)
	{
		// create a vector
		openfpm::vector<aggregate<double,double[2],double[2][2]>> v1;
		openfpm::vector<aggregate<double[2]>> v2;

		// Point
		aggregate<double[2]> p;
		p.get<0>()[0] = 1.0;
		p.get<0>()[1] = 2.0;

		aggregate<double,double[2],double[2][2]> pa;
		pa.get<0>() = 1.0;

		pa.get<1>()[0] = 1.0;
		pa.get<1>()[1] = 1.0;

		pa.get<2>()[0][0] = 1.0;
		pa.get<2>()[0][1] = 1.0;
		pa.get<2>()[1][0] = 1.0;
		pa.get<2>()[1][1] = 1.0;

		// Add test

		for (size_t i = 0 ; i < NADD ; i++)
		{
			v1.add(pa);
			v2.add(p);
		}

		timer tg;
		tg.start();

		for (size_t i = 0 ; i < NADD ; i++)
		{
			read_write(v1,v2,i);
		}

		tg.stop();

		times_g[i] = tg.getwct();

		timer tga;
		tga.start();

		double * prp = (double *)v1.getPointer<0>();
		double * pos = (double *)v2.getPointer<0>();

		for (size_t i = 0 ; i < NADD ; i++)
		{
			read_write_lin(pos,(struct ele *)prp,i);
		}

		tga.stop();

		times[i] = tga.getwct();
	}

	double mean;
	double dev;
	standard_deviation(times_g,mean,dev);

	double mean_;
	double dev_;
	standard_deviation(times,mean_,dev_);

	report_vector_funcs.graphs.put("performance.vector_layout(0).y.data.mean",mean_/mean);

	// Deviation od x/y = x/y^2 dy + 1/y dx

	report_vector_funcs.graphs.put("performance.vector_layout(0).y.data.dev",mean_/(mean*mean)*dev + dev_ / mean );

	report_vector_funcs.graphs.put("performance.vector_layout(1).funcs.nele",NADD);
	report_vector_funcs.graphs.put("performance.vector_layout(1).funcs.name","read_write_inte");

	for (size_t i = 0 ; i < N_STAT+1 ; i++)
	{
		// create a vector
		openfpm::vector<aggregate<double,double[2],double[2][2]>,HeapMemory,memory_traits_inte> v1;
		openfpm::vector<aggregate<double[2]>,HeapMemory,memory_traits_inte> v2;

		// Point
		aggregate<double[2]> p;
		p.get<0>()[0] = 1.0;
		p.get<0>()[1] = 2.0;

		aggregate<double,double[2],double[2][2]> pa;
		pa.get<0>() = 1.0;

		pa.get<1>()[0] = 1.0;
		pa.get<1>()[1] = 1.0;

		pa.get<2>()[0][0] = 1.0;
		pa.get<2>()[0][1] = 1.0;
		pa.get<2>()[1][0] = 1.0;
		pa.get<2>()[1][1] = 1.0;

		// Add test

		for (size_t i = 0 ; i < NADD ; i++)
		{
			v1.add(pa);
			v2.add(p);
		}

		timer tg;
		tg.start();

		for (size_t i = 0 ; i < NADD ; i++)
		{
			read_write(v1,v2,i);
		}

		tg.stop();

		times2_g[i] = tg.getwct();
		int sz = v1.size();

		timer tga;
		tga.start();

		double * prp0 = (double *)v1.getPointer<0>();
		double * prp1 = (double *)v1.getPointer<1>();
		double * prp2 = (double *)v1.getPointer<2>();

		double * pos = (double *)v2.getPointer<0>();

		for (size_t i = 0 ; i < NADD ; i++)
		{
			read_write_inte(pos,prp0,prp1,prp2,i,sz);
		}

		tga.stop();

		times2[i] = tga.getwct();
	}

	double mean2;
	double dev2;
	standard_deviation(times2_g,mean2,dev2);

	double mean2_;
	double dev2_;
	standard_deviation(times2,mean2_,dev2_);

	report_vector_funcs.graphs.put("performance.vector_layout(1).y.data.mean",mean2_/mean2);

	// Deviation od x/y = x/y^2 dy + 1/y dx

	report_vector_funcs.graphs.put("performance.vector_layout(1).y.data.dev",mean2_/(mean2*mean2)*dev2 + dev2_ / mean2 );
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

	report_vector_funcs.graphs.put("graphs.graph(1).type","line");
	report_vector_funcs.graphs.add("graphs.graph(1).title","Vector read write");
	report_vector_funcs.graphs.add("graphs.graph(1).x.title","Layout");
	report_vector_funcs.graphs.add("graphs.graph(1).y.title","Time seconds");
	report_vector_funcs.graphs.add("graphs.graph(1).y.data(0).source","performance.vector_layout(#).y.data.mean");
	report_vector_funcs.graphs.add("graphs.graph(1).x.data(0).source","performance.vector_layout(#).funcs.name");
	report_vector_funcs.graphs.add("graphs.graph(1).y.data(0).title","Actual");
	report_vector_funcs.graphs.add("graphs.graph(1).interpolation","lines");

	boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
	boost::property_tree::write_xml("vector_performance_funcs.xml", report_vector_funcs.graphs,std::locale(),settings);

	GoogleChart cg;

	std::string file_xml_ref(test_dir);
	file_xml_ref += std::string("/openfpm_data/vector_performance_funcs_ref.xml");

	StandardXMLPerformanceGraph("vector_performance_funcs.xml",file_xml_ref,cg);

	addUpdtateTime(cg,1);
	createCommitFile("pdata");

	cg.write("vector_performance_funcs.html");
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_DATA_SRC_VECTOR_VECTOR_PERFORMANCE_TEST_HPP_ */



#endif /* VECTOR_PERFORMANCE_TEST_HPP_ */
