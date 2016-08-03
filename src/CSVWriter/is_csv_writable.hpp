/*
 * is_csv_writable.hpp
 *
 *  Created on: Jul 18, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_CSVWRITER_IS_CSV_WRITABLE_HPP_
#define OPENFPM_IO_SRC_CSVWRITER_IS_CSV_WRITABLE_HPP_



template<typename T>
struct is_csv_writable
{
	enum
	{
		value = false
	};
};

template<>
struct is_csv_writable<float>
{
	enum
	{
		value = true
	};
};

template<>
struct is_csv_writable<double>
{
	enum
	{
		value = true
	};
};

template<>
struct is_csv_writable<char>
{
	enum
	{
		value = true
	};
};

template<>
struct is_csv_writable<unsigned char>
{
	enum
	{
		value = true
	};
};

template<>
struct is_csv_writable<short>
{
	enum
	{
		value = true
	};
};

template<>
struct is_csv_writable<unsigned short>
{
	enum
	{
		value = true
	};
};


template<>
struct is_csv_writable<int>
{
	enum
	{
		value = true
	};
};

template<>
struct is_csv_writable<unsigned int>
{
	enum
	{
		value = true
	};
};

template<>
struct is_csv_writable<long int>
{
	enum
	{
		value = true
	};
};

template<>
struct is_csv_writable<unsigned long int>
{
	enum
	{
		value = true
	};
};


template<>
struct is_csv_writable<bool>
{
	enum
	{
		value = true
	};
};


#endif /* OPENFPM_IO_SRC_CSVWRITER_IS_CSV_WRITABLE_HPP_ */
