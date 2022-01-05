//
// Created by jstark on 28.12.21.
//
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "CSVReader/CSVReader.hpp"

BOOST_AUTO_TEST_SUITE(CSVReaderTestSuite)
BOOST_AUTO_TEST_CASE(csv_reader_int)
		{
				std::cout << "CWD = " << get_cwd() << std::endl;
				std::string csv_file = std::string("../../../openfpm_io/test_data/integer.csv");
				// Read csv file into vector while linearizing
				openfpm::vector<int> v_lin; // Vector to which csv file will be read to
				size_t m, n; // Number of rows m and columns n
				read_csv_to_vector(csv_file, v_lin, m, n);
				
				BOOST_CHECK(m == 4);
				BOOST_CHECK(n == 3);
				BOOST_CHECK(m * n == v_lin.size());
				
				for(int i = 0; i < v_lin.size() / n; ++i)
				{
					BOOST_CHECK( v_lin.get(i * n) == i + 1);
					BOOST_CHECK( v_lin.get(i * n + 1) == (i + 1) * 2);
					BOOST_CHECK( v_lin.get(i * n + 2) == v_lin.get(i * n) * v_lin.get(i * n + 1));
				}
			
		}


BOOST_AUTO_TEST_CASE(csv_reader_char)
		{
				// Read csv file into vector while linearizing
				const std::string csv_file = "../../../openfpm_io/test_data/char.csv"; // csv file to be read
				openfpm::vector<std::string> v_lin; // Vector to which csv file will be read to
				size_t m, n; // Number of rows m and columns n
				
				read_csv_to_vector(csv_file, v_lin, m, n);
				
				BOOST_CHECK(m == 5);
				BOOST_CHECK(n == 2);
				BOOST_CHECK(m * n == v_lin.size());
				
				openfpm::vector<std::string> col1 = {"a", "b", "c", "d", "e"};
				openfpm::vector<std::string> col2 = {"antilope", "ballena", "camel", "dackel", "elefant"};
				
				for(int i = 0; i < v_lin.size() / n; ++i)
				{
					BOOST_CHECK(col1.get(i) == v_lin.get(i * n));
					BOOST_CHECK(col2.get(i) == v_lin.get(i * n + 1));
				}
			
		}
BOOST_AUTO_TEST_SUITE_END()
