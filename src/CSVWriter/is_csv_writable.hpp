/*
 * is_csv_writable.hpp
 *
 *  Created on: Jul 18, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_CSVWRITER_IS_CSV_WRITABLE_HPP_
#define OPENFPM_IO_SRC_CSVWRITER_IS_CSV_WRITABLE_HPP_


//! Indicate if the property T is writable in CSV
template<typename T>
struct is_csv_writable
{
	//! as default if not recognized is not writable
	enum
	{
		value = false
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<float>
{
	//! float is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<double>
{
	//! double is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<char>
{
	//! char is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<unsigned char>
{
	//! unsigned char is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<short>
{
	//! short is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<unsigned short>
{
	//! unsigned is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<int>
{
	//! int is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<unsigned int>
{
	//! unsigned int is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<long int>
{
	//! long int is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<unsigned long int>
{
	//! unsigned long int is writable
	enum
	{
		value = true
	};
};

//! Indicate if the property T is writable in CSV
template<>
struct is_csv_writable<bool>
{
	//! bool is writable
	enum
	{
		value = true
	};
};


#endif /* OPENFPM_IO_SRC_CSVWRITER_IS_CSV_WRITABLE_HPP_ */
