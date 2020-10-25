/*
 * multi_array_ref_subarray_openfpm.hpp
 *
 *  Created on: Jul 1, 2018
 *      Author: i-bird
 */

#ifndef MULTI_ARRAY_REF_SUBARRAY_OPENFPM_HPP_
#define MULTI_ARRAY_REF_SUBARRAY_OPENFPM_HPP_

#include "multi_array_view_openfpm.hpp"
#include "multi_array_ref_base_openfpm.hpp"

/*! \brief return the dimension of the sub_array
 *
 *
 */
template<typename vmpl>
struct subar_dim
{
	typedef typename boost::mpl::at<vmpl,boost::mpl::int_<0>>::type type;
};

namespace openfpm {
namespace detail {
namespace multi_array {

//
// const_sub_array
//    multi_array's proxy class to allow multiple overloads of
//    operator[] in order to provide a clean multi-dimensional array
//    interface.
template <typename T, std::size_t NumDims, typename vector, typename TPtr>
class const_sub_array_openfpm : public openfpm::detail::multi_array::multi_array_impl_base_openfpm<T,NumDims,vector>
{
  typedef openfpm::detail::multi_array::multi_array_impl_base_openfpm<T,NumDims,vector> super_type;

  typedef typename boost::mpl::accumulate<vector,
								   typename boost::mpl::int_<1>,
								   typename boost::mpl::multiplies<typename boost::mpl::_2,typename boost::mpl::_1>  >::type size_ct;


public:

  typedef typename super_type::value_type value_type;
  typedef typename super_type::const_reference const_reference;
  typedef typename super_type::const_iterator const_iterator;
  typedef typename super_type::const_reverse_iterator const_reverse_iterator;
  typedef typename super_type::element element;
  typedef typename super_type::size_type size_type;
  typedef typename super_type::difference_type difference_type;
  typedef typename super_type::index index;

  // template typedefs
  template <std::size_t NDims>
  struct const_array_view {
    typedef openfpm::detail::multi_array::const_multi_array_view_openfpm<T,NDims> type;
  };

  template <std::size_t NDims>
  struct array_view {
    typedef openfpm::detail::multi_array::multi_array_view_openfpm<T,NDims> type;
  };

  // Allow default copy constructor as well.

  template <typename Ovector>
  const_sub_array_openfpm (const const_sub_array_openfpm<T,NumDims,Ovector>& rhs)
  :base_(rhs.base_), strides_(rhs.strides_)
  {}

  // const_sub_array always returns const types, regardless of its own
  // constness.
  inline __device__ __host__  const_reference operator[](index idx) const
  {
    return super_type::access(boost::type<const_reference>(),idx,strides(),base_);
  }

  template <typename IndexList>
  const element& operator()(const IndexList& indices) const
  {
    boost::function_requires<boost::CollectionConcept<IndexList> >();
    return super_type::access_element(boost::type<const element&>(),
                                      indices,origin(),strides());
  }

  // see generate_array_view in base.hpp
  template <int NDims>
  __device__ __host__   typename const_array_view<NDims>::type
  operator[](const boost::detail::multi_array::index_gen<NumDims,NDims>& indices) const
  {
    typedef typename const_array_view<NDims>::type return_type;
    return super_type::generate_array_view(boost::type<return_type>(),
                                           indices,
                                           base_);
  }

  template <typename OPtr>
  bool operator!=(const const_sub_array_openfpm<T,NumDims,OPtr>& rhs) const {
    return !(*this == rhs);
  }

  template <typename OPtr>
  bool operator>(const const_sub_array_openfpm<T,NumDims,OPtr>& rhs) const {
    return rhs < *this;
  }

  template <typename OPtr>
  bool operator<=(const const_sub_array_openfpm<T,NumDims,OPtr>& rhs) const {
    return !(*this > rhs);
  }

  template <typename OPtr>
  bool operator>=(const const_sub_array_openfpm<T,NumDims,OPtr>& rhs) const {
    return !(*this < rhs);
  }

  TPtr origin() const { return base_; }
  inline __host__ __device__ size_type size() const { return boost::mpl::at<vector,boost::mpl::int_<0>>::type::value; }
  size_type max_size() const { return num_elements(); }
  bool empty() const { return size() == 0; }
  size_type num_dimensions() const { return NumDims; }
  inline __host__ __device__ const index* strides() const { return strides_; }

  size_type num_elements() const
  {
    return size_ct::type::value;
  }

  __device__ __host__ const_sub_array_openfpm (TPtr base, const index* strides)
  :base_(base), strides_(strides)
  {}

#ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
protected:
  template <typename,std::size_t,typename> friend class value_accessor_n_openfpm;
  template <typename,std::size_t,typename,typename> friend class const_sub_array_openfpm;
#else
public:  // Should be protected
#endif


  TPtr base_;
  const index* strides_;
private:
  // const_sub_array cannot be assigned to (no deep copies!)
  const_sub_array_openfpm& operator=(const const_sub_array_openfpm&);
};

//
// sub_array
//    multi_array's proxy class to allow multiple overloads of
//    operator[] in order to provide a clean multi-dimensional array
//    interface.
template <typename T, std::size_t NumDims, typename vector>
class sub_array_openfpm : public const_sub_array_openfpm<T,NumDims,vector,T*>
{
  typedef const_sub_array_openfpm<T,NumDims,vector,T*> super_type;
public:
  typedef typename super_type::element element;
  typedef typename super_type::reference reference;
  typedef typename super_type::index index;
  typedef typename super_type::size_type size_type;
  typedef typename super_type::iterator iterator;
  typedef typename super_type::reverse_iterator reverse_iterator;
  typedef typename super_type::const_reference const_reference;
  typedef typename super_type::const_iterator const_iterator;
  typedef typename super_type::const_reverse_iterator const_reverse_iterator;
  typedef int yes_is_multi_array;

  // template typedefs
  template <std::size_t NDims>
  struct const_array_view {
    typedef openfpm::detail::multi_array::const_multi_array_view_openfpm<T,NDims> type;
  };

  template <std::size_t NDims>
  struct array_view {
    typedef openfpm::detail::multi_array::multi_array_view_openfpm<T,NDims> type;
  };

  // Assignment from other ConstMultiArray types.
  template <typename ConstMultiArray>
  __device__ __host__ inline sub_array_openfpm& operator=(const ConstMultiArray& other)
  {
#ifdef SE_CLASS1

    // make sure the dimensions agree
    BOOST_ASSERT(other.num_dimensions() == this->num_dimensions());

#endif
    // iterator-based copy
//    std::copy(other.begin(),other.end(),begin());

	for (int i = 0 ; i < (int)other.size() ; i++)
	{this->operator[](i) = other[i];}
    return *this;
  }

  __device__ __host__ sub_array_openfpm& operator=(const sub_array_openfpm& other)
  {
#ifdef SE_CLASS1
      // make sure the dimensions agree
      BOOST_ASSERT(other.num_dimensions() == this->num_dimensions());
//      BOOST_ASSERT(std::equal(other.shape(),
//                              other.shape()+this->num_dimensions(),
//                              this->shape()));
#endif
      // iterator-based copy
      //std::copy(other.begin(),other.end(),begin());

    	  for (int i = 0 ; i < (int)other.size() ; i++)
    	  {this->operator[](i) = other[i];}
    return *this;
  }

  __device__ __host__ T* origin() { return this->base_; }
  __device__ __host__ const T* origin() const { return this->base_; }

  __device__ __host__ reference operator[](index idx) {
    return super_type::access(boost::type<reference>(),
                              idx,
                              this->strides(),
                              this->base_);
  }

  __device__ __host__ iterator begin()
  {
    return iterator(*this->index_bases(),origin(),
                    this->shape(),this->strides(),this->index_bases());
  }

  __device__ __host__ iterator end()
  {
    return iterator(*this->index_bases()+(index)*this->shape(),origin(),
                    this->shape(),this->strides(),this->index_bases());
  }

  // RG - rbegin() and rend() written naively to thwart MSVC ICE.
  reverse_iterator rbegin() {
    reverse_iterator ri(end());
    return ri;
  }

  reverse_iterator rend() {
    reverse_iterator ri(begin());
    return ri;
  }

  //
  // proxies
  //

  template <class IndexList>
  const element& operator()(const IndexList& indices) const {
    boost::function_requires<
      boost::CollectionConcept<IndexList> >();
    return super_type::operator()(indices);
  }

  const_reference operator[](index idx) const {
    return super_type::operator[](idx);
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

#ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
private:
  template <typename,std::size_t,typename> friend class value_accessor_n_openfpm;
#endif
public: // should be private

  inline __device__ __host__ sub_array_openfpm (T* base,
            									const index* strides)
  :super_type(base,strides)
  {}

};

} // namespace multi_array
} // namespace detail
//
// traits classes to get sub_array types
//
template <typename Array, int N, typename vector>
class subarray_gen_openfpm {
  typedef typename Array::element element;
public:
  typedef openfpm::detail::multi_array::sub_array_openfpm<element,N,vector> type;
};

template <typename Array, int N, typename vector>
class const_subarray_gen_openfpm {
  typedef typename Array::element element;
public:
  typedef openfpm::detail::multi_array::const_sub_array_openfpm<element,N,vector> type;
};
} // namespace boost


#endif /* MULTI_ARRAY_REF_SUBARRAY_OPENFPM_HPP_ */
