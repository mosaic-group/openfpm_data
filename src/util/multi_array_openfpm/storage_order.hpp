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
#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <vector>

namespace openfpm
{

class c_storage_order;
class fortran_storage_order;
class ofp_storage_order;

template <std::size_t NumDims>
class general_storage_order
{
public:
  typedef detail::multi_array::size_type size_type;
  template <typename OrderingIter, typename AscendingIter>
  general_storage_order(OrderingIter ordering,
                        AscendingIter ascending) {
    boost::detail::multi_array::copy_n(ordering,NumDims,ordering_.begin());
  }

  // RG - ideally these would not be necessary, but some compilers
  // don't like template conversion operators.  I suspect that not
  // too many folk will feel the need to use customized
  // storage_order objects, I sacrifice that feature for compiler support.
  general_storage_order(const c_storage_order&) {
    for (size_type i=0; i != NumDims; ++i) {
      ordering_[i] = NumDims - 1 - i;
    }
  }

  general_storage_order(const fortran_storage_order&) {
    for (size_type i=0; i != NumDims; ++i) {
      ordering_[i] = i;
    }
  }

  general_storage_order(const ofp_storage_order&)
  {
	    ordering_[0] = 0;
	    for (size_type i=1; i != NumDims; ++i)
	    {ordering_[i] = NumDims - i;}
  }

  size_type ordering(size_type dim) const { return ordering_[dim]; }


  bool operator==(general_storage_order const& rhs) const
  {
    return (ordering_ == rhs.ordering_);
  }

protected:
  openfpm::array<size_type,NumDims> ordering_;
};

class ofp_storage_order
{
  typedef detail::multi_array::size_type size_type;
public:
  // This is the idiom for creating your own custom storage orders.
  // Not supported by all compilers though!
#ifndef __MWERKS__ // Metrowerks screams "ambiguity!"
  template <std::size_t NumDims>
  operator general_storage_order<NumDims>() const
  {
    openfpm::array<size_type,NumDims> ordering;

    ordering[0] = 0;
    for (size_type i=1; i != NumDims; ++i)
    {ordering[i] = NumDims - i;}
    return general_storage_order<NumDims>(ordering.begin());
  }
#endif
};

class c_storage_order
{
  typedef detail::multi_array::size_type size_type;
public:
  // This is the idiom for creating your own custom storage orders.
  // Not supported by all compilers though!
#ifndef __MWERKS__ // Metrowerks screams "ambiguity!"
  template <std::size_t NumDims>
  operator general_storage_order<NumDims>() const {
    openfpm::array<size_type,NumDims> ordering;
    openfpm::array<bool,NumDims> ascending;

    for (size_type i=0; i != NumDims; ++i) {
      ordering[i] = NumDims - 1 - i;
    }
    return general_storage_order<NumDims>(ordering.begin());
  }
#endif
};

class fortran_storage_order
{
  typedef detail::multi_array::size_type size_type;
public:
  // This is the idiom for creating your own custom storage orders.
  // Not supported by all compilers though!
#ifndef __MWERKS__ // Metrowerks screams "ambiguity!"
  template <std::size_t NumDims>
  operator general_storage_order<NumDims>() const
  {
    openfpm::array<size_type,NumDims> ordering;

    for (size_type i=0; i != NumDims; ++i) {
      ordering[i] = i;
    }
    return general_storage_order<NumDims>(ordering.begin());
  }
#endif
};

} // namespace openfpm


#endif /* STORAGE_ORDER_HPP_ */
