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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "util/performance/performance_util.hpp"

constexpr int N_STAT = 32;
constexpr int N_STAT_SMALL = 32;
constexpr int N_TRY = 8;

#ifdef PERFORMANCE_TEST

GoogleChart cg;
const char * test_dir;

// Declaration of functions


void load_and_combine(std::string file, openfpm::vector<openfpm::vector<float>> & y, openfpm::vector<float> & per_times);
void speedup_calculate(openfpm::vector<openfpm::vector<float>> & y_ref_sup, openfpm::vector<openfpm::vector<float>> & y, openfpm::vector<openfpm::vector<float>> & y_ref ,openfpm::vector<std::string> & yn);

BOOST_AUTO_TEST_SUITE( performance )

//// Include tests ////////

#include "Grid/performance/grid_performance_tests.hpp"
//#include "Vector/performance/vector_performance_test.hpp"

BOOST_AUTO_TEST_SUITE_END()


#endif

#endif /* OPENFPM_DATA_SRC_PERFORMANCE_HPP_ */
