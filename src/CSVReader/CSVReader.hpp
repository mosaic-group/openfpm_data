//
// Created by jstark on 23.12.21.
//

#ifndef OPENFPM_IO_CSVREADER_HPP
#define OPENFPM_IO_CSVREADER_HPP

#include <iostream>
#include <sstream>
#include <string>

#include "Vector/map_vector.hpp"
#include "util/PathsAndFiles.hpp"
#include "VCluster/VCluster.hpp"

/**@brief Converts string into template type T.
 *
 * @tparam T Destination type.
 * @param string_to_convert String that will be streamed into a variable of type T
 * @return T-type variable storing the content of string_to_convert
 */
template <typename T>
T string_to_type(const std::string & string_to_convert)
{
	T x_new_type;
	std::istringstream istream_string_to_convert(string_to_convert);
	istream_string_to_convert >> x_new_type;
	return x_new_type;
}


/**@brief Csv file is read by process rank 0 and the elements are saved into a linearized vector. The vector is
 * broadcasted to the other processes.
 *
 * @tparam T Type that elements from csv file shall be converted to.
 * @param input_file std::string containing the path to the input csv file.
 * @param output_vector Empty openfpm::vector<T> to which linearized matrix from csv file will be written.
 * @param m Size_t variable to which number of rows of the csv file will be returned.
 * @param n Size_t variable to which number of columns of the csv file will be returned.
 */
template <typename T>
void read_csv_to_vector(const std::string & input_file, openfpm::vector<T> & output_vector, size_t & m, size_t & n)
{
	auto & v_cl = create_vcluster();
	if (v_cl.rank() == 0)
	{
		if(!check_if_file_exists(input_file))
		{
			std::cout << "Cannot find < " <<  input_file << " >. Please check path specification. Aborting..."
					<< std::endl;
			abort();
		}
		std::ifstream file(input_file); // Create file pointer and open file
		std::string line, element;
		m = 0;
		n = 0;
		size_t elems = 0;
		while (getline(file, line)) // Read entire row and store as one string in line
		{
			std::stringstream stream_line(line); // Needed to break line into elements later one
			while (getline(stream_line, element, ',')) // Read every column of line and store content as string in element
			{
				output_vector.add(string_to_type<T>(element)); // Convert element from string to type T and append to
				// output_vector
				elems += 1; // Increase n by one for each element read
			}
			m += 1; // Increase m by one for each row read
		}
		if(m >= 1) n = elems / m; // If at least one row present, divide number of elements read by number of rows to get
		// the number of columns
	}
	v_cl.Bcast(output_vector, 0);
	v_cl.execute();
}

#endif //OPENFPM_IO_CSVREADER_HPP
