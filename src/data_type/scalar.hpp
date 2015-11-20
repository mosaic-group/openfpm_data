/*
 * scalar.hpp
 *
 *  Created on: Feb 3, 2015
 *      Author: i-bird
 */

#ifndef SCALAR_HPP_
#define SCALAR_HPP_

#include <boost/fusion/container/vector.hpp>

/*! \brief it define a scalar value compatible with grid_cpu , grid_gpu, vector, graph ...
 *
 * \tparam T scalar type
 *
 */

template <typename  T>
class scalar
{
public:

  typedef boost::fusion::vector<T> type;

  type data;

  static const unsigned int ele = 0;
  static const unsigned int max_prop = 1;

  // Setter method

  inline void set(T s_)	{boost::fusion::at_c<0>(data) = s_;};

  // getter method

  T& get()	{return boost::fusion::at_c<ele>(data);};
};


#endif /* SCALAR_HPP_ */
