/*
 * Matrix.hpp
 *
 *  Created on: Jun 3, 2015
 *      Author: Pietro Incardona
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_


#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include "boost/multi_array.hpp"
#include "Grid/grid_key.hpp"

/*! \brief This class implement an NxN (dense) matrix
 *
 * Be carefull when you allocate locally it take memory from the stack
 *
 * \tparam dim dimensionality
 * \tparam T type of the space
 *
 * \warning Matrix is untested so it remain as concept usage until test are not introduced
 *
 */

template<unsigned int dim ,typename T> class Matrix
{
	public:

	typedef T coord_type;

	//! boost fusion that store the point
	typedef boost::fusion::vector<T[dim][dim]> type;

	//! structure that store the data of the point
	type data;

	//! Property id of the point
	static const unsigned int mat = 0;

	/*! \brief Get coordinate
	 *
	 * \param i row
	 * \param j colums
	 * \return the matrix element
	 *
	 */

	inline T get(size_t i,size_t j) const
	{
		return boost::fusion::at_c<mat>(data)[i][j];
	}

	/*! \brief Get coordinate
	 *
	 * \param i row
	 * \param j colum
	 * \return the value
	 *
	 */

	inline T& get(size_t i,size_t j)
	{
		return boost::fusion::at_c<mat>(data)[i][j];
	}

	/*! \brief operator= between Matrix
	 *
	 * \param m Matrix
	 *
	 */
	inline Matrix<dim,T> & operator=(const Matrix<dim,T> & m)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			for (size_t j = 0 ; j < dim ; j++)
				get(i,j) = m.get(i,j);
		}

		return *this;
	}

	/*! \brief Set to zero the point coordinate
	 *
	 *
	 */
	void zero()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			for (size_t j = 0 ; j < dim ; j++)
			{
				get(i,j) = 0;
			}
		}
	}

	/*! \brief Point constructor from point
	 *
	 * \param p the point
	 *
	 */
	Matrix(const Matrix<dim,T> && p)
	{
	    for(size_t i = 0; i < dim ; i++)
	    {
			for (size_t j = 0 ; j < dim ; j++)
			{
				get(i,j) = p.get(i,j);
			}
	    }
	}

	/*! \brief Point constructor from point
	 *
	 * \param p the point
	 *
	 */
	Matrix(const Matrix<dim,T> & p)
	{
	    for(size_t i = 0; i < dim ; i++)
	    {
			for (size_t j = 0 ; j < dim ; j++)
			{
				get(i,j) = p.get(i,j);
			}
	    }
	}

	/*! \brief Constructor from an array
	 *
	 * \param p array with the coordinate of the point
	 *
	 */
	Matrix(const T (&p)[dim][dim])
	{
	    for(size_t i = 0; i < dim ; i++)
	    {
			for (size_t j = 0 ; j < dim ; j++)
			{
				get(i,j) = p[i][j];
			}
	    }
	}

	//! Default contructor
	Matrix()
	{}

	/*! \brief Identity matrix
	 *
	 * \return the identity matrix
	 *
	 */
	inline static Matrix<dim,T> identity()
	{
		Matrix<dim,T> ret;

		for (size_t i = 0 ; i < dim ; i++)
		{
			for (size_t j = 0 ; j < dim ; j++)
			{
				/* coverity[dead_error_line] */
				ret.get(i,j) = (i == j)?1:0;
			}
		}

		return ret;
	}

	//! 1 property
	static const unsigned int max_prop = 1;

	//! dimension of the matrix (it is a square matrix)
	static const unsigned int dims = dim;
};


#endif /* MATRIX_HPP_ */
