#ifndef POINT_TEST_HPP
#define POINT_TEST_HPP

#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include "boost/multi_array.hpp"
#include "base_type.hpp"

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
};

// tranformed

template<typename T> class Point_test
{
public:
  
  typedef boost::fusion::vector<T,T,T,T,T[3],T[3][3]> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;
  
  static const unsigned int x = 0;
  static const unsigned int y = 1;
  static const unsigned int z = 2;
  static const unsigned int s = 3;
  static const unsigned int v = 4;
  static const unsigned int t = 5;
  static const unsigned int max_prop = 6;
  
  // Setter method

  inline void setx(T x_)	{boost::fusion::at_c<0>(data) = x_;};
  inline void sety(T y_)	{boost::fusion::at_c<1>(data) = y_;};
  inline void setz(T z_)	{boost::fusion::at_c<2>(data) = z_;};
  inline void sets(T s_)	{boost::fusion::at_c<3>(data) = s_;};
  
  // getter method

  template<unsigned int i> inline typename boost::fusion::result_of::at<type, boost::mpl::int_<i> >::type get()	{return boost::fusion::at_c<i>(data);};
};


#endif
