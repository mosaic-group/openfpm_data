/*
 * storage_order.hpp
 *
 *  Created on: Jul 1, 2018
 *      Author: i-bird
 */

#ifndef STORAGE_ORDER_HPP_
#define STORAGE_ORDER_HPP_

#include "types.hpp"
#include "array_openfpm.hpp"
#include "boost/multi_array/algorithm.hpp"
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/at.hpp>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <vector>

namespace openfpm
{

template<unsigned int NumDims, unsigned int n, typename T> struct c_storage_order;
template<unsigned int NumDims, unsigned int n, typename T> struct ofp_storage_order;
template<unsigned int NumDims, unsigned int n, typename T> struct fortran_storage_order;

template<unsigned int NumDims, unsigned int n=NumDims-1, typename T=boost::mpl::vector<>>
struct c_storage_order
{
  typedef typename boost::mpl::push_back<typename c_storage_order<NumDims,n-1,T>::value, boost::mpl::int_<NumDims-1-n>>::type value;
};

template <unsigned int NumDims, typename T>
struct c_storage_order<NumDims, 0, T>
{
  typedef typename boost::mpl::push_back<T,boost::mpl::int_<NumDims-1>>::type value;
};

template<unsigned int NumDims, unsigned int n=NumDims-1, typename T=boost::mpl::vector<>>
struct fortran_storage_order
{
  typedef typename boost::mpl::push_back<typename fortran_storage_order<NumDims,n-1,T>::value, boost::mpl::int_<n>>::type value;
};

template <unsigned int NumDims, typename T>
struct fortran_storage_order<NumDims, 0, T>
{
  typedef typename boost::mpl::push_back<T,boost::mpl::int_<0>>::type value;
};

template<unsigned int NumDims, unsigned int n=NumDims-1, typename T=boost::mpl::vector<>>
struct ofp_storage_order
{
  typedef typename boost::mpl::push_back<typename ofp_storage_order<NumDims,n-1,T>::value, boost::mpl::int_<NumDims-n>>::type value;
};

template <unsigned int NumDims, typename T>
struct ofp_storage_order<NumDims, 0, T>
{
  typedef typename boost::mpl::push_back<T,boost::mpl::int_<0>>::type value;
};

} // namespace openfpm


#endif /* STORAGE_ORDER_HPP_ */
