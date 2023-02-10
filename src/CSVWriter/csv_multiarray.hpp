/*
 * csv_multiarray_copy.hpp
 *
 *  Created on: Jun 20, 2015
 *      Author: i-bird
 */

#ifndef CSV_MULTIARRAY_COPY_HPP_
#define CSV_MULTIARRAY_COPY_HPP_



/*! \brief This class is an helper to produce csv headers from multi-array
 *
 * Usage:
 *
 * \code{.cpp}
 *
 * float src[3];
 *
 * std::stringstream str;
 * csv_col_str<float[3]> cp(std::string("test"),str);
 *
 * std::cout << str.str() << "\n";
 *
 * \endcode
 *
 * Will produce ",test[0],test[1],test[2]"
 *
 */
template<typename T>
struct csv_col_str
{
	inline csv_col_str(std::string prp, std::stringstream & str)
	{
		str << "," << prp;
	}
};

//! Partial specialization for N=1 1D-Array
template<typename T,size_t N1>
struct csv_col_str<T[N1]>
{
	inline csv_col_str(std::string prp, std::stringstream & str)
	{
		for (size_t i = 0 ; i < N1 ; i++)
			str << "," << prp << "_" << "[" << i << "]";
	}
};

//! Partial specialization for N=2 2D-Array
template<typename T,size_t N1,size_t N2>
struct csv_col_str<T[N1][N2]>
{
	inline csv_col_str(std::string prp, std::stringstream & str)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				str << "," << prp << "_"  << "[" << i1 << "]" << "[" << i2 << "]";
			}
		}
	}
};

//! Partial specialization for N=3
template<typename T,size_t N1,size_t N2,size_t N3>
struct csv_col_str<T[N1][N2][N3]>
{
	inline csv_col_str(std::string prp, std::stringstream & str)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					str << "," << prp << "_"  << "[" << i1 << "]" << "[" << i2 << "]" << "[" << i3 << "]";
				}
			}
		}
	}
};

//! Partial specialization for N=4
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4>
struct csv_col_str<T[N1][N2][N3][N4]>
{
	inline csv_col_str(std::string prp, std::stringstream & str)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						str << "," << prp << "_"  << "[" << i1 << "]" << "[" << i2 << "]" << "[" << i3 << "]" << "[" << i4 << "]";
					}
				}
			}
		}
	}
};


/*! \brief This class is an helper to produce csv data from multi-array
 *
 * Usage:
 *
 * \code{.cpp}
 *
 * float src[3] = {1.0,2.0,3.0};
 *
 * std::stringstream str;
 * csv_value_str<float[3]> cp(src,str);
 *
 * std::cout << str.str() << "\n";
 *
 * \endcode
 *
 * Will produce ",1.0,2.0,3.0"
 *
 */
template<typename T, bool is_writable>
struct csv_value_str
{
	inline csv_value_str(T & v, std::stringstream & str)
	{
		str << "," << v;
	}
};

//! Partial specialization for N=1 1D-Array
template<typename T,size_t N1, bool is_writable>
struct csv_value_str<T[N1], is_writable>
{
	template<typename ArrObject>
	inline csv_value_str(const ArrObject v, std::stringstream & str)
	{
		for (size_t i = 0 ; i < N1 ; i++)
			str << "," << v[i];
	}
};

//! Partial specialization for N=2 2D-Array
template<typename T,size_t N1,size_t N2, bool is_writable>
struct csv_value_str<T[N1][N2], is_writable>
{
	template<typename ArrObject>
	inline csv_value_str(const ArrObject v, std::stringstream & str)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				str << "," << v[i1][i2];
			}
		}
	}
};

//! Partial specialization for N=3
template<typename T,size_t N1,size_t N2,size_t N3, bool is_writable>
struct csv_value_str<T[N1][N2][N3], is_writable>
{
	template<typename ArrObject>
	inline csv_value_str(const  ArrObject v, std::stringstream & str)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					str << "," <<  v[i1][i2][i3];
				}
			}
		}
	}
};

//! Partial specialization for N=4
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4, bool is_writable>
struct csv_value_str<T[N1][N2][N3][N4],is_writable>
{
	template<typename ArrObject>
	inline csv_value_str(const ArrObject v, std::stringstream & str)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						str << "," << v[i1][i2][i3][i4];
					}
				}
			}
		}
	}
};

//! Partial specialization for unknown property
template<typename T>
struct csv_value_str<T,false>
{
	inline csv_value_str(const T v, std::stringstream & str)
	{
		str << "," << 0.0;
	}
};

#endif /* CSV_MULTIARRAY_COPY_HPP_ */
