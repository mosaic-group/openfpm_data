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

#include "boost/multi_array/collection_concept.hpp"
#include "boost/multi_array/concept_checks.hpp"
#include "boost/multi_array/storage_order.hpp"
#include "boost_multi_array_view_openfpm.hpp"
#include "boost/multi_array/algorithm.hpp"
#include "boost/type_traits/is_integral.hpp"
#include "boost/utility/enable_if.hpp"
#include "boost_array_openfpm.hpp"
#include "boost/concept_check.hpp"
#include "boost/functional.hpp"
#include "boost/limits.hpp"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>

#include "util/boost/boost_multi_array_base_openfpm.hpp"
#include "util/boost/boost_multi_array_iterator_openfpm.hpp"
#include "util/boost/boost_multi_array_subarray_openfpm.hpp"

namespace boost {

  // RG - This is to make things work with VC++. So sad, so sad.
  class c_storage_order;
  class fortran_storage_order;
  class ofp_storage_order;

  template <std::size_t NumDims>
  class general_storage_order_ofp
  {
  public:
    typedef detail::multi_array::size_type size_type;
    template <typename OrderingIter, typename AscendingIter>
    general_storage_order_ofp(OrderingIter ordering,
                          AscendingIter ascending) {
      boost::detail::multi_array::copy_n(ordering,NumDims,ordering_.begin());
      boost::detail::multi_array::copy_n(ascending,NumDims,ascending_.begin());
    }

    // RG - ideally these would not be necessary, but some compilers
    // don't like template conversion operators.  I suspect that not
    // too many folk will feel the need to use customized
    // storage_order objects, I sacrifice that feature for compiler support.
    general_storage_order_ofp(const c_storage_order&) {
      for (size_type i=0; i != NumDims; ++i) {
        ordering_[i] = NumDims - 1 - i;
      }
      ascending_.assign(true);
    }

    general_storage_order_ofp(const fortran_storage_order&) {
      for (size_type i=0; i != NumDims; ++i) {
        ordering_[i] = i;
      }
      ascending_.assign(true);
    }

    general_storage_order_ofp(const ofp_storage_order&) {
      ordering_[NumDims - 1] = 0;

      for (size_type i=0; i != NumDims; ++i) {
        ordering_[i] = i + 1;
      }
      ascending_.assign(true);
    }

    size_type ordering(size_type dim) const { return ordering_[dim]; }
    bool ascending(size_type dim) const { return ascending_[dim]; }

    bool all_dims_ascending() const {
      return std::accumulate(ascending_.begin(),ascending_.end(),true,
                      std::logical_and<bool>());
    }

    bool operator==(general_storage_order_ofp const& rhs) const {
      return (ordering_ == rhs.ordering_) &&
        (ascending_ == rhs.ascending_);
    }

  protected:
    boost::array<size_type,NumDims> ordering_;
    boost::array<bool,NumDims> ascending_;
  };

  class ofp_storage_order
  {
    typedef detail::multi_array::size_type size_type;
  public:
    // This is the idiom for creating your own custom storage orders.
    // Not supported by all compilers though!
#ifndef __MWERKS__ // Metrowerks screams "ambiguity!"
    template <std::size_t NumDims>
    operator general_storage_order<NumDims>() const {
      boost::array<size_type,NumDims> ordering;
      boost::array<bool,NumDims> ascending;

      ordering[0] = 0;
      ascending[0] = true;
      for (size_type i=1; i != NumDims; ++i) {
        ordering[i] = NumDims - i;
        ascending[i] = true;
      }
      return general_storage_order<NumDims>(ordering.begin(),
                                            ascending.begin());
    }
#endif
  };

template <typename T, std::size_t NumDims,
          typename TPtr = const T*>
class const_multi_array_ref_openfpm : public detail::multi_array::multi_array_impl_base_openfpm<T,NumDims>
{
	typedef detail::multi_array::multi_array_impl_base_openfpm<T,NumDims> super_type;

	public:

    typedef typename super_type::value_type value_type;
    typedef typename super_type::const_reference const_reference;
    typedef typename super_type::const_iterator const_iterator;
    typedef typename super_type::const_reverse_iterator const_reverse_iterator;
    typedef typename super_type::element element;
    typedef typename super_type::size_type size_type;
    typedef typename super_type::difference_type difference_type;
    typedef typename super_type::index index;
    typedef typename super_type::extent_range extent_range;
    typedef general_storage_order<NumDims> storage_order_type;

    // template typedefs
    template <std::size_t NDims>
    struct const_array_view_openfpm
	{
    	typedef boost::detail::multi_array::const_multi_array_view_openfpm<T,NDims> type;
    };

    template <std::size_t NDims>
    struct array_view
	{
      typedef boost::detail::multi_array::multi_array_view<T,NDims> type;
    };

	#ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
    	// make const_multi_array_ref a friend of itself
    	template <typename,std::size_t,typename>
    	friend class const_multi_array_ref;
  	#endif

    // This ensures that const_multi_array_ref types with different TPtr
    // types can convert to each other
    template <typename OPtr>
    const_multi_array_ref_openfpm(const const_multi_array_ref_openfpm<T,NumDims,OPtr>& other)
      : base_(other.base_), storage_(other.storage_),
        extent_list_(other.extent_list_),
        stride_list_(other.stride_list_),
        index_base_list_(other.index_base_list_),
        origin_offset_(other.origin_offset_),
        directional_offset_(other.directional_offset_),
        num_elements_(other.num_elements_)  {  }

    template <typename ExtentList>
    explicit const_multi_array_ref_openfpm(TPtr base, const ExtentList& extents)
    :base_(base), storage_(c_storage_order())
    {
      boost::function_requires<
        CollectionConcept<ExtentList> >();

      index_base_list_.assign(0);
      init_multi_array_ref(extents.begin());
    }

    template <typename ExtentList>
    explicit const_multi_array_ref_openfpm(TPtr base, const ExtentList& extents,
                         const general_storage_order<NumDims>& so)
    :base_(base), storage_(so)
    {
      boost::function_requires<
        CollectionConcept<ExtentList> >();

      index_base_list_.assign(0);
      init_multi_array_ref(extents.begin());
    }

    explicit const_multi_array_ref_openfpm(TPtr base,
                           const detail::multi_array::
                           extent_gen<NumDims>& ranges)
    :base_(base), storage_(c_storage_order())
    {

      init_from_extent_gen(ranges);
    }

    explicit const_multi_array_ref_openfpm(TPtr base,
                             const detail::multi_array::
                             extent_gen<NumDims>& ranges,
                             const general_storage_order<NumDims>& so)
    :base_(base), storage_(so)
    {
      init_from_extent_gen(ranges);
    }

    template <class InputIterator>
    void assign(InputIterator begin, InputIterator end) {
      boost::function_requires<InputIteratorConcept<InputIterator> >();

      InputIterator in_iter = begin;
      T* out_iter = base_;
      std::size_t copy_count=0;
      while (in_iter != end && copy_count < num_elements_) {
        *out_iter++ = *in_iter++;
        copy_count++;
      }
    }

    template <class BaseList>
  #ifdef BOOST_NO_SFINAE
    void
  #else
    typename
    disable_if<typename boost::is_integral<BaseList>::type,void >::type
  #endif // BOOST_NO_SFINAE
    reindex(const BaseList& values) {
      boost::function_requires<
        CollectionConcept<BaseList> >();
      boost::detail::multi_array::
        copy_n(values.begin(),num_dimensions(),index_base_list_.begin());
      origin_offset_ =
        this->calculate_origin_offset(stride_list_,extent_list_,
                                storage_,index_base_list_);
    }

    void reindex(index value) {
      index_base_list_.assign(value);
      origin_offset_ =
        this->calculate_origin_offset(stride_list_,extent_list_,
                                storage_,index_base_list_);
    }

    template <typename SizeList>
    void reshape(const SizeList& extents) {
      boost::function_requires<
        CollectionConcept<SizeList> >();
      BOOST_ASSERT(num_elements_ ==
                   std::accumulate(extents.begin(),extents.end(),
                                   size_type(1),std::multiplies<size_type>()));

      std::copy(extents.begin(),extents.end(),extent_list_.begin());
      this->compute_strides(stride_list_,extent_list_,storage_);

      origin_offset_ = this->calculate_origin_offset(stride_list_,extent_list_, storage_,index_base_list_);
    }

    size_type num_dimensions() const { return NumDims; }

    size_type size() const { return extent_list_.front(); }

    // given reshaping functionality, this is the max possible size.
    size_type max_size() const { return num_elements(); }

    bool empty() const { return size() == 0; }

    __device__ __host__ const size_type* shape() const {
      return extent_list_.data();
    }

    __device__ __host__ const index* strides() const {
      return stride_list_.data();
    }

    __device__ __host__ const element* origin() const { return base_+origin_offset_; }
    __device__ __host__ const element* data() const { return base_; }

    size_type num_elements() const { return num_elements_; }

    __device__ __host__ const index* index_bases() const {
      return index_base_list_.data();
    }


    const storage_order_type& storage_order() const {
      return storage_;
    }

    template <typename IndexList>
    const element& operator()(IndexList indices) const {
      boost::function_requires<
        CollectionConcept<IndexList> >();
      return super_type::access_element(boost::type<const element&>(),
                                        indices,origin(),
                                        shape(),strides(),index_bases());
    }

    // Only allow const element access
    __device__ __host__ const_reference operator[](index idx) const {
      return super_type::access(boost::type<const_reference>(),
                                idx,origin(),
                                shape(),strides(),index_bases());
    }

    // see generate_array_view in base.hpp
    template <int NDims>
    __device__ __host__ typename const_array_view_openfpm<NDims>::type
    operator[](const detail::multi_array::
               index_gen<NumDims,NDims>& indices)
      const {
      typedef typename const_array_view_openfpm<NDims>::type return_type;
      return
        super_type::generate_array_view(boost::type<return_type>(),
                                        indices,
                                        shape(),
                                        strides(),
                                        index_bases(),
                                        origin());
    }

    const_iterator begin() const {
      return const_iterator(*index_bases(),origin(),
                            shape(),strides(),index_bases());
    }

    const_iterator end() const {
      return const_iterator(*index_bases()+(index)*shape(),origin(),
                            shape(),strides(),index_bases());
    }

    const_reverse_iterator rbegin() const {
      return const_reverse_iterator(end());
    }

    const_reverse_iterator rend() const {
      return const_reverse_iterator(begin());
    }


    template <typename OPtr>
    bool operator==(const
                    const_multi_array_ref_openfpm<T,NumDims,OPtr>& rhs)
      const {
      if(std::equal(extent_list_.begin(),
                    extent_list_.end(),
                    rhs.extent_list_.begin()))
        return std::equal(begin(),end(),rhs.begin());
      else return false;
    }

    template <typename OPtr>
    bool operator<(const
                   const_multi_array_ref_openfpm<T,NumDims,OPtr>& rhs)
      const {
      return std::lexicographical_compare(begin(),end(),rhs.begin(),rhs.end());
    }

    template <typename OPtr>
    bool operator!=(const
                    const_multi_array_ref_openfpm<T,NumDims,OPtr>& rhs)
      const {
      return !(*this == rhs);
    }

    template <typename OPtr>
    bool operator>(const
                   const_multi_array_ref_openfpm<T,NumDims,OPtr>& rhs)
      const {
      return rhs < *this;
    }

    template <typename OPtr>
    bool operator<=(const
                   const_multi_array_ref_openfpm<T,NumDims,OPtr>& rhs)
      const {
      return !(*this > rhs);
    }

    template <typename OPtr>
    bool operator>=(const
                   const_multi_array_ref_openfpm<T,NumDims,OPtr>& rhs)
      const {
      return !(*this < rhs);
    }


  #ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
  protected:
  #else
  public:
  #endif

    typedef boost::array_openfpm<size_type,NumDims> size_list;
    typedef boost::array_openfpm<index,NumDims> index_list;

    // This is used by multi_array, which is a subclass of this
    void set_base_ptr(TPtr new_base) { base_ = new_base; }


    // This constructor supports multi_array's default constructor
    // and constructors from multi_array_ref, subarray, and array_view
    explicit
    const_multi_array_ref_openfpm(TPtr base,
                          const storage_order_type& so,
                          const index * index_bases,
                          const size_type* extents) :
      base_(base), storage_(so), origin_offset_(0), directional_offset_(0)
   {
     // If index_bases or extents is null, then initialize the corresponding
     // private data to zeroed lists.
     if(index_bases) {
       boost::detail::multi_array::
         copy_n(index_bases,NumDims,index_base_list_.begin());
     } else {
       std::fill_n(index_base_list_.begin(),NumDims,0);
     }
     if(extents) {
       init_multi_array_ref(extents);
     } else {
       boost::array<index,NumDims> extent_list;
       extent_list.assign(0);
       init_multi_array_ref(extent_list.begin());
     }
   }


    TPtr base_;
    storage_order_type storage_;
    size_list extent_list_;
    index_list stride_list_;
    index_list index_base_list_;
    index origin_offset_;
    index directional_offset_;
    size_type num_elements_;

  private:
    // const_multi_array_ref cannot be assigned to (no deep copies!)
    const_multi_array_ref_openfpm& operator=(const const_multi_array_ref_openfpm & other);

    void init_from_extent_gen(const
                          detail::multi_array::
                          extent_gen<NumDims>& ranges) {

      typedef boost::array<index,NumDims> extent_list;

      // get the index_base values
      std::transform(ranges.ranges_.begin(),ranges.ranges_.end(),
                index_base_list_.begin(),
                boost::mem_fun_ref(&extent_range::start));

      // calculate the extents
      extent_list extents;
      std::transform(ranges.ranges_.begin(),ranges.ranges_.end(),
                extents.begin(),
                boost::mem_fun_ref(&extent_range::size));

      init_multi_array_ref(extents.begin());
    }


  #ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
  protected:
  #else
  public:
  #endif
    // RG - move me!
    template <class InputIterator>
    void init_multi_array_ref(InputIterator extents_iter) {
      boost::function_requires<InputIteratorConcept<InputIterator> >();

      boost::detail::multi_array::
        copy_n(extents_iter,num_dimensions(),extent_list_.begin());

      // Calculate the array size
      num_elements_ = std::accumulate(extent_list_.begin(),extent_list_.end(),
                              size_type(1),std::multiplies<size_type>());

      this->compute_strides(stride_list_,extent_list_,storage_);

      origin_offset_ =
        this->calculate_origin_offset(stride_list_,extent_list_,
                                storage_,index_base_list_);
      directional_offset_ =
        this->calculate_descending_dimension_offset(stride_list_,extent_list_,
                                              storage_);
    }
  };

template <typename T, std::size_t NumDims>
class multi_array_ref_openfpm :
  public const_multi_array_ref_openfpm<T,NumDims,T*>
{
  typedef const_multi_array_ref_openfpm<T,NumDims,T*> super_type;
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
  struct const_array_view_openfpm {
    typedef boost::detail::multi_array::const_multi_array_view_openfpm<T,NDims> type;
  };

  template <std::size_t NDims>
  struct array_view_openfpm {
    typedef boost::detail::multi_array::multi_array_view_openfpm<T,NDims> type;
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
                           const general_storage_order_ofp<NumDims>& so) :
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

  multi_array_ref_openfpm & bind_ref(const multi_array_ref_openfpm & other) {
    if (&other != this) {

        this->base_ = other.base_;
        this->storage_ = other.storage_;
        this->extent_list_ = other.extent_list_;
        this->stride_list_ = other.stride_list_;
        this->index_base_list_ = other.index_base_list_;
        this->origin_offset_ = other.origin_offset_;
        this->directional_offset_ = other.directional_offset_;
        this->num_elements_ = other.num_elements_;

    }
    return *this;
  }

  /* \brief Set the internal pointer
   *
   * \param base internal pointer
   *
   */
  void set_pointer(void * base)
  {
	  this->base_ = static_cast<T *>(base);
  }

  multi_array_ref_openfpm & operator=(multi_array_ref_openfpm && other) {

    this->base_ = other.base_;
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

  __device__ __host__ element* origin() { return super_type::base_+super_type::origin_offset_; }

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


  __device__ __host__ reference operator[](index idx) {
    return super_type::access(boost::type<reference>(),
                              idx,origin(),
                              this->shape(),this->strides(),
                              this->index_bases());
  }


  // See note attached to generate_array_view in base.hpp
  template <int NDims>
  typename array_view_openfpm<NDims>::type
  operator[](const detail::multi_array::
             index_gen<NumDims,NDims>& indices) {
    typedef typename array_view_openfpm<NDims>::type return_type;
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
  typename const_array_view_openfpm<NDims>::type
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
