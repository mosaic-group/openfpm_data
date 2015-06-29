/*
 * Point_orig.hpp
 *
 *  Created on: Jun 7, 2015
 *      Author: i-bird
 */

#ifndef POINT_ORIG_HPP_
#define POINT_ORIG_HPP_

/*! \brief Definition of a class Point in plain C++ and boost::vector for testing purpose
 *
 * Definition of a class Point in plain C++ and boost::vector for testing purpose
 *
 *	\param T base type of all the fields
 *
 */

template<typename T> class Point_orig
{
public:

  T x;
  T y;
  T z;

  T s;

  T v[3];
  T t[3][3];

  inline void setx(T x_)	{x = x_;};
  inline void sety(T y_)	{y = y_;};
  inline void setz(T z_)	{z = z_;};
  inline void sets(T s_)	{s = s_;};
};



#endif /* POINT_ORIG_HPP_ */
