//
// Created by jstark on 23.12.21.
//

#ifndef OPENFPM_IO_CSVREADER_HPP
#define OPENFPM_IO_CSVREADER_HPP

#include <iostream>
#include <sstream>
#include <string>

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


/**@brief Reads csv file and saves elements as linearized vector.
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
	std::ifstream file(input_file); // Create file pointer and open file
	std::string line, element;
	m = 0;
	n = 0;
	
	while (getline(file, line)) // Read entire row and store as one string in line
	{
		std::stringstream stream_line(line); // Needed to break line into elements later one
		
		while (getline(stream_line, element, ',')) // Read every column of line and store content as string in element
		{
			output_vector.add(string_to_type<T>(element)); // Convert element from string to type T and append to
			// output_vector
			n += 1; // Increase n by one for each element read
		}
		m += 1; // Increase m by one for each row read
	}
	n /= m; // Divide number of elements read by number of rows to get number of columns
}

#endif //OPENFPM_IO_CSVREADER_HPP
