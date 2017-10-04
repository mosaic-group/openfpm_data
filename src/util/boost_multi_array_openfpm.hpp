// Copyright 2002 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Boost.MultiArray Library
//  Authors: Ronald Garcia
//           Jeremy Siek
//           Andrew Lumsdaine
//
//  Modified by Pietro incardona for openfpm
//
//  See http://www.boost.org/libs/multi_array for documentation.

#ifndef OPENFPM_DATA_SRC_GRID_BOOST_MULTI_ARRAY_OPENFPM_HPP_
#define OPENFPM_DATA_SRC_GRID_BOOST_MULTI_ARRAY_OPENFPM_HPP_


//
// multi_array_ref.hpp - code for creating "views" of array data.
//

#include "boost/multi_array/base.hpp"
#include "boost/multi_array/collection_concept.hpp"
#include "boost/multi_array/concept_checks.hpp"
#include "boost/multi_array/iterator.hpp"
#include "boost/multi_array/storage_order.hpp"
#include "boost/multi_array/subarray.hpp"
#include "boost/multi_array/view.hpp"
#include "boost/multi_array/algorithm.hpp"
#include "boost/type_traits/is_integral.hpp"
#include "boost/utility/enable_if.hpp"
#include "boost/array.hpp"
#include "boost/concept_check.hpp"
#include "boost/functional.hpp"
#include "boost/limits.hpp"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>

namespace boost {


template <typename T, std::size_t NumDims>
class multi_array_ref_openfpm :
  public const_multi_array_ref<T,NumDims,T*>
{
  typedef const_multi_array_ref<T,NumDims,T*> super_type;
public:
  typedef typename super_type::value_type value_type;
  typedef typename super_type::reference reference;
  typedef typename super_type::iterator iterator;
  typedef typename super_type::reverse_iterator reverse_iterator;
  typedef typename super_type::const_reference const_reference;
  typedef typename super_type::const_iterator const_iterator;
  typedef typename super_type::const_reverse_iterator const_reverse_iterator;
  typedef typename super_type::element element;
  typedef typename super_type::size_type size_type;
  typedef typename super_type::difference_type difference_type;
  typedef typename super_type::index index;
  typedef typename super_type::extent_range extent_range;

  typedef typename super_type::storage_order_type storage_order_type;
  typedef typename super_type::index_list index_list;
  typedef typename super_type::size_list size_list;

  template <std::size_t NDims>
  struct const_array_view {
    typedef boost::detail::multi_array::const_multi_array_view<T,NDims> type;
  };

  template <std::size_t NDims>
  struct array_view {
    typedef boost::detail::multi_array::multi_array_view<T,NDims> type;
  };

  template <class ExtentList>
  explicit multi_array_ref_openfpm(T* base, const ExtentList& extents) :
    super_type(base,extents) {
    boost::function_requires<
      CollectionConcept<ExtentList> >();
  }

  template <class ExtentList>
  explicit multi_array_ref_openfpm(T* base, const ExtentList& extents,
                           const general_storage_order<NumDims>& so) :
    super_type(base,extents,so) {
    boost::function_requires<
      CollectionConcept<ExtentList> >();
  }


  explicit multi_array_ref_openfpm(T* base,
                           const detail::multi_array::
                           extent_gen<NumDims>& ranges) :
    super_type(base,ranges) { }


  explicit multi_array_ref_openfpm(T* base,
                           const detail::multi_array::
                           extent_gen<NumDims>&
                             ranges,
                           const general_storage_order<NumDims>& so) :
    super_type(base,ranges,so) { }


  // Assignment from other ConstMultiArray types.
  template <typename ConstMultiArray>
  multi_array_ref_openfpm & operator=(const ConstMultiArray& other) {
    function_requires<
      multi_array_concepts::
      ConstMultiArrayConcept<ConstMultiArray,NumDims> >();

    // make sure the dimensions agree
    BOOST_ASSERT(other.num_dimensions() == this->num_dimensions());
    BOOST_ASSERT(std::equal(other.shape(),other.shape()+this->num_dimensions(),
                            this->shape()));
    // iterator-based copy
    std::copy(other.begin(),other.end(),this->begin());
    return *this;
  }

  multi_array_ref_openfpm & operator=(const multi_array_ref_openfpm & other) {
    if (&other != this) {
      // make sure the dimensions agree

      BOOST_ASSERT(other.num_dimensions() == this->num_dimensions());
      BOOST_ASSERT(std::equal(other.shape(),
                              other.shape()+this->num_dimensions(),
                              this->shape()));
      // iterator-based copy
      std::copy(other.begin(),other.end(),this->begin());
    }
    return *this;
  }

  multi_array_ref_openfpm & operator=(multi_array_ref_openfpm && other) {

    this->base_ = other.base;
    this->storage_ = other.storage_;
    this->extent_list_ = other.extent_list_;
    this->stride_list_ = other.stride_list_;
    this->index_base_list_ = other.index_base_list_;
    this->origin_offset_ = other.origin_offset_;
    this->directional_offset_ = other.directional_offset_;
    this->num_elements_ = other.num_elements_;

    return *this;
  }

  void swap(multi_array_ref_openfpm & other)
  {
	T* base_tmp = this->base_;
	this->base_ = other.base_;
	other.base_ = base_tmp;

    storage_order_type storage_tmp = this->storage_;
    this->storage_ = other.storage_;
    other.storage_ = storage_tmp;

    size_list extent_list_tmp = this->extent_list_;
    this->extent_list_ = other.extent_list_;
    other.extent_list_ = extent_list_tmp;

    index_list stride_list_tmp = this->stride_list_;
    this->stride_list_ = other.stride_list_;
    other.stride_list_ = stride_list_tmp;

    index_list index_base_list_tmp = this->index_base_list_;
    this->index_base_list_ = other.index_base_list_;
    other.index_base_list_ = index_base_list_tmp;

    index origin_offset_tmp = this->origin_offset_;
    this->origin_offset_ = other.origin_offset_;
    other.origin_offset_ = origin_offset_tmp;

    index directional_offset_tmp = this->directional_offset_;
    this->directional_offset_ = other.directional_offset_;
    other.directional_offset_ = directional_offset_tmp;

    size_type num_elements_tmp = this->num_elements_;
    this->num_elements_ = other.num_elements_;
    other.num_elements_ = num_elements_tmp;
  }

  element* origin() { return super_type::base_+super_type::origin_offset_; }

  element* data() { return super_type::base_; }

  template <class IndexList>
  element& operator()(const IndexList& indices) {
    boost::function_requires<
      CollectionConcept<IndexList> >();
    return super_type::access_element(boost::type<element&>(),
                                      indices,origin(),
                                      this->shape(),this->strides(),
                                      this->index_bases());
  }


  reference operator[](index idx) {
    return super_type::access(boost::type<reference>(),
                              idx,origin(),
                              this->shape(),this->strides(),
                              this->index_bases());
  }


  // See note attached to generate_array_view in base.hpp
  template <int NDims>
  typename array_view<NDims>::type
  operator[](const detail::multi_array::
             index_gen<NumDims,NDims>& indices) {
    typedef typename array_view<NDims>::type return_type;
    return
      super_type::generate_array_view(boost::type<return_type>(),
                                      indices,
                                      this->shape(),
                                      this->strides(),
                                      this->index_bases(),
                                      origin());
  }


  iterator begin() {
    return iterator(*this->index_bases(),origin(),this->shape(),
                    this->strides(),this->index_bases());
  }

  iterator end() {
    return iterator(*this->index_bases()+(index)*this->shape(),origin(),
                    this->shape(),this->strides(),
                    this->index_bases());
  }

  // rbegin() and rend() written naively to thwart MSVC ICE.
  reverse_iterator rbegin() {
    reverse_iterator ri(end());
    return ri;
  }

  reverse_iterator rend() {
    reverse_iterator ri(begin());
    return ri;
  }

  // Using declarations don't seem to work for g++
  // These are the proxies to work around this.

  const element* origin() const { return super_type::origin(); }
  const element* data() const { return super_type::data(); }

  template <class IndexList>
  const element& operator()(const IndexList& indices) const {
    boost::function_requires<
      CollectionConcept<IndexList> >();
    return super_type::operator()(indices);
  }

  const_reference operator[](index idx) const {
    return super_type::access(boost::type<const_reference>(),
                              idx,origin(),
                              this->shape(),this->strides(),
                              this->index_bases());
  }

  // See note attached to generate_array_view in base.hpp
  template <int NDims>
  typename const_array_view<NDims>::type
  operator[](const detail::multi_array::
             index_gen<NumDims,NDims>& indices)
    const {
    return super_type::operator[](indices);
  }

  const_iterator begin() const {
    return super_type::begin();
  }

  const_iterator end() const {
    return super_type::end();
  }

  const_reverse_iterator rbegin() const {
    return super_type::rbegin();
  }

  const_reverse_iterator rend() const {
    return super_type::rend();
  }

protected:
  // This is only supplied to support multi_array's default constructor
  explicit multi_array_ref_openfpm(T* base,
                           const storage_order_type& so,
                           const index* index_bases,
                           const size_type* extents) :
    super_type(base,so,index_bases,extents) { }

};

} // namespace boost


#endif /* OPENFPM_DATA_SRC_GRID_BOOST_MULTI_ARRAY_OPENFPM_HPP_ */
