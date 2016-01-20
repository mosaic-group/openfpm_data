/*
 * performance.hpp
 *
 *  Created on: Jan 11, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_PERFORMANCE_HPP_
#define OPENFPM_DATA_SRC_PERFORMANCE_HPP_

#include "Plot/GoogleChart.hpp"
#include "timer.hpp"

//#define N_STAT 256
//#define N_STAT_SMALL 32
//#define N_TRY 8

#define N_STAT 1
#define N_STAT_SMALL 1
#define N_TRY 1

#ifdef PERFORMANCE_TEST

GoogleChart cg;
const char * test_dir;

// Declaration of functions


void load_and_combine(std::string file, openfpm::vector<openfpm::vector<float>> & y, openfpm::vector<float> & per_times);
void speedup_calculate(openfpm::vector<openfpm::vector<float>> & y_ref_sup, openfpm::vector<openfpm::vector<float>> & y, openfpm::vector<openfpm::vector<float>> & y_ref ,openfpm::vector<std::string> & yn);

BOOST_AUTO_TEST_SUITE( performance )

//// Include tests ////////

#include "Grid/grid_performance_tests.hpp"
#include "Vector/vector_performance_test.hpp"

BOOST_AUTO_TEST_SUITE_END()

#define MEASURE_SET 5

/*! \brief Execute a command getting its output
 *
 * \param cmd command to execute
 *
 *
 * \return the output string
 *
 */
std::string exec(const char* cmd)
{
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe)
    	return "ERROR";

    char buffer[128];
    std::string result = "";
    while (!feof(pipe.get()))
    {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    return result;
}

/*! \brief It store the sha-1 git commit
 *
 *
 */
struct configuration
{
	char data_commit[41];
	char pdata_commit[41];
	char devices_commit[41];
	char io_commit[41];
	char vcluster_commit[41];
	char numeric_commit[41];

	/*! \brief Convert to HTML
	 *
	 * \param i configuration number
	 *
	 */
	std::string toHTML(size_t i)
	{
		std::stringstream html;

		html << "<div style=\" border-radius: 25px; padding: 20px 20px 20px 20px ;\" >\n\
	  	<div style=\"background-color:white; border-radius: 25px;border: 2px solid #000000; padding: 20px;\" >\n\
	  	<h4>config " << i << "</h4>\n";
	    html << "<strong>OpenFPM numeric:</strong> " << numeric_commit << "<br>\n";
	    html << "<strong>OpenFPM pdata:</strong> " << pdata_commit << "<br>\n";
	    html << "<strong>OpenFPM vcluster:</strong> " << vcluster_commit << "<br>\n";
	    html << "<strong>OpenFPM io:</strong> " << io_commit << "<br>\n";
	    html << "<strong>OpenFPM pdata:</strong> " << pdata_commit << "<br>\n";
	    html << "<strong>OpenFPM data:</strong> " << data_commit << "<br>\n";
	    html << "<strong>OpenFPM devices:</strong> " << devices_commit << "<br>\n";
	    html << "</div>\n";
	    html << "</div>\n";

	    return html.str();
	}
};

/*! \brief Get the actual configuration
 *
 *
 *
 */
void getConfig(configuration & conf)
{
	// take the head of the actual version of data
	std::string out = exec("git rev-parse HEAD");

	// copy the actual commit version
	for (size_t i = 0 ; i < 40 ; i++)
		conf.data_commit[i] = out[i];
	out[40] = 0;

	// take the head of the actual version of devices
	std::string out2 = exec("cd ../../openfpm_devices ; git rev-parse HEAD");

	// copy the actual commit version
	for (size_t i = 0 ; i < 40 ; i++)
		conf.devices_commit[i] = out2[i];
	out2[40] = 0;

	// take the head of the actual version of io
	std::string out3 = exec("cd ../../openfpm_io ; git rev-parse HEAD");

	// copy the actual commit version
	for (size_t i = 0 ; i < 40 ; i++)
		conf.io_commit[i] = out3[i];
	out[40] = 0;

	// take the head of the actual version of vcluster
	std::string out4 = exec("cd ../../openfpm_vcluster ; git rev-parse HEAD");

	// copy the actual commit version
	for (size_t i = 0 ; i < 40 ; i++)
		conf.vcluster_commit[i] = out4[i];
	out[40] = 0;

	// take the head of the actual version of pdata
	std::string out5 = exec("cd ../../ ; git rev-parse HEAD");

	// copy the actual commit version
	for (size_t i = 0 ; i < 40 ; i++)
		conf.pdata_commit[i] = out5[i];
	out[40] = 0;

	// take the head of the actual version of numerics
	std::string out6 = exec("cd ../../openfpm_numerics ; git rev-parse HEAD");

	// copy the actual commit version
	for (size_t i = 0 ; i < 40 ; i++)
		conf.numeric_commit[i] = out6[i];
	out[40] = 0;
}

/*! \brief Load the previous test result combine with the actual result and save
 *
 * \param file file that contain the previous result
 * \param commit vector with the previous commit SHA-1
 *
 */
void load_and_combine_commit(std::string file, openfpm::vector<configuration> & commits)
{
	// Load the previous measure and combine the previous measure with the actual measure
	commits.clear();
	commits.load(file);

	// remove the oldest commit
	if (commits.size() >= MEASURE_SET)
		commits.remove(0);

	commits.add();

	getConfig(commits.last());

	bool saved = commits.save(file);
	if (saved == false)
		std::cerr << __FILE__ << ":" << __LINE__ << " error saving previous commit \n";
}

/*! \brief Load the previous test result combine with the actual result and save
 *
 * load the previous result combining the the actual results eliminating the oldest
 * dataset
 *
 * \param file file that contain the previous performance test results
 * \param y vector loaded with the previous performance results
 * \param per_times for each test the result
 *
 */
void load_and_combine(std::string file, openfpm::vector<openfpm::vector<float>> & y, openfpm::vector<float> & per_times)
{
	// Load the previous measure and combine the previous measure with the actual measure
	y.clear();

	y.load(file);

	if (per_times.size() == 0)
		return;

	if (y.size() == 0)
		y.resize(per_times.size());

	for(size_t i = 0 ; i < y.size() ; i++)
	{
		if (y.get(i).size() >= MEASURE_SET)
			y.get(i).remove(0);
		y.get(i).add(per_times.get(i));
	}

	y.save(file);
}

/*! \brief Calculate the speedup
 *
 * Given the references times and the actual results, calculate the relative
 * speedup
 *
 * \param y_ref_sup vector with the calculated speed-up
 * \param y vector with the actual performance results
 * \param y_ref vector with the reference times
 * \param Output vector usefull for plotting r_ref_sup with GoogleChart
 *
 */
void speedup_calculate(openfpm::vector<openfpm::vector<float>> & y_ref_sup, openfpm::vector<openfpm::vector<float>> & y, openfpm::vector<openfpm::vector<float>> & y_ref ,openfpm::vector<std::string> & yn)
{
	yn.clear();

	for (size_t i = 0 ; i < y_ref.size() ; i++)
	{
		// Get the minimum and maximum of the reference time
		float min = 100.0;
		float max = 0.0;
		float average = 0.0;

		for (size_t j = 0 ; j < y_ref.get(i).size() ; j++)
		{
			if (y_ref.get(i).get(j) < min)
				min = y_ref.get(i).get(j);

			if (y_ref.get(i).get(j) > max)
				max = y_ref.get(i).get(j);

			average += y_ref.get(i).get(j);
		}
		average /= y_ref.get(i).size();

		// calculate speedup percentage

		y_ref_sup.add();

		for (size_t j = 0 ; j < y.get(i).size() ; j++)
		{
			y_ref_sup.last().add((average/y.get(i).get(j) - 1.0) * 100.0);
			yn.add("config " + std::to_string(j));
		}

		y_ref_sup.last().add((average/min - 1.0) * 100.0);
		yn.add("interval");
		y_ref_sup.last().add((average/max - 1.0) * 100.0);
		yn.add("interval");
	}
}

void write_performance_report()
{

	openfpm::vector<configuration> pcommit;

	// Get the directory of the performance test files
	std::string per_dir(test_dir);

	// Load the previous git commit SHA1 of the previous tests
	load_and_combine_commit(per_dir + std::string("/prev_commit"),pcommit);

	std::string config_list;

	config_list += std::string("<h3>Highest configuration number is the latest version</h3>\n");

	for (size_t i = 0 ; i < pcommit.size() ; i++)
		config_list += pcommit.get(i).toHTML(i);

	cg.addHTML(config_list);
	cg.write(test_dir + std::string("/openfpm_data_performance.html"));
}

#endif

#endif /* OPENFPM_DATA_SRC_PERFORMANCE_HPP_ */
