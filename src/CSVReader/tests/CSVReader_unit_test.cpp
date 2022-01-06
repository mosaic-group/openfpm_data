//
// Created by jstark on 28.12.21.
//
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "CSVReader/CSVReader.hpp"

BOOST_AUTO_TEST_SUITE(CSVReaderTestSuite)
BOOST_AUTO_TEST_CASE(csv_reader_int_test)
		{
#ifdef OPENFPM_PDATA
			std::string csv_file = std::string("openfpm_io/test_data/integer.csv");
#else
			std::string csv_file = std::string("test_data/integer.csv");
#endif
			std::cout << "CWD = " << get_cwd() << std::endl;
			// Read csv file into vector while linearizing
			openfpm::vector<int> v_lin; // Vector to which csv file will be read to
			size_t m, n; // Number of rows m and columns n
			read_csv_to_vector(csv_file, v_lin, m, n);
		
			auto & v_cl = create_vcluster();
			
			std::cout << "My rank is " << v_cl.rank() << std::endl;
			auto v_iter = v_lin.getIterator();
			while(v_iter.isNext())
			{
				auto key = v_iter.get();
				std::cout << "Element number " << key << " of rank " <<  v_cl.rank() << " is " << v_lin.get(key) <<
				std::endl;
				++v_iter;
			}
			
			BOOST_CHECK(m == 4);
			BOOST_CHECK(n == 3);
			BOOST_CHECK(m * n == v_lin.size());

			for(int i = 0; i < v_lin.size() / n; ++i)
			{
				BOOST_CHECK( v_lin.get(i * n) == i + 1);
				BOOST_CHECK( v_lin.get(i * n + 1) == (i + 1) * 2);
				BOOST_CHECK( v_lin.get(i * n + 2) == v_lin.get(i * n) * v_lin.get(i * n + 1));
			}
		
			std::cout << "Rank " << v_cl.rank() << " reaches end of unit test." << std::endl;
		}


BOOST_AUTO_TEST_CASE(csv_reader_char_test)
		{
#ifdef OPENFPM_PDATA
			std::string csv_file = std::string("openfpm_io/test_data/char.csv");
#else
			std::string csv_file = std::string("test_data/char.csv");
#endif
			// Read csv file into vector while linearizing
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
