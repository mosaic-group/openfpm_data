/*
 * multi_array_ref_openfpm.hpp
 *
 *  This is an heavily modified boost::multi_array version
 * 
 *  Created on: Jun 29, 2018
 *      Author: i-bird
 */

#ifndef MULTI_ARRAY_REF_OPENFPM_HPP_
#define MULTI_ARRAY_REF_OPENFPM_HPP_

#include "util/cuda_util.hpp"
#include "boost/multi_array/collection_concept.hpp"
#include "boost/multi_array/concept_checks.hpp"
#include "boost/multi_array/algorithm.hpp"
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/multiplies.hpp>
#include <boost/concept_check.hpp>
#include "array_openfpm.hpp"
#include "types.hpp"
#include "multi_array_ref_base_openfpm.hpp"
#include "multi_array_ref_subarray_openfpm.hpp"
#include "storage_order.hpp"
#include "multi_array_iterator_openfpm.hpp"
#include "util/common.hpp"

namespace openfpm {



template <typename T, std::size_t NumDims, typename vector, typename TPtr, typename StorageOrder>
class const_multi_array_ref_openfpm : public detail::multi_array::multi_array_impl_base_openfpm<T,NumDims,vector>
{
	typedef detail::multi_array::multi_array_impl_base_openfpm<T,NumDims,vector> super_type;

	typedef typename boost::mpl::accumulate<vector,
								   typename boost::mpl::int_<1>,
								   typename boost::mpl::multiplies<typename boost::mpl::_2,typename boost::mpl::_1>  >::type size_ct;

public:

/*    typedef typename super_type::value_type value_type;*/
    typedef typename super_type::const_reference const_reference;
    typedef typename super_type::const_iterator const_iterator;
/*    typedef typename super_type::const_reverse_iterator const_reverse_iterator;*/
    typedef T element;
    typedef size_t size_type;
/*    typedef typename super_type::difference_type difference_type;*/
    typedef typename super_type::index index;
/*    typedef typename super_type::extent_range extent_range;*/

    template <typename ExtentType>
    explicit const_multi_array_ref_openfpm(TPtr base, const ExtentType& extents)
    :base_(base)
    {
    	init_multi_array_ref(extents);
    }


    template <class InputIterator>
    void assign(InputIterator begin, InputIterator end)
    {
    	boost::function_requires<boost::InputIteratorConcept<InputIterator> >();

    	InputIterator in_iter = begin;
    	T* out_iter = base_;
    	std::size_t copy_count=0;
    	while (in_iter != end && copy_count < num_elements_)
    	{
    		*out_iter++ = *in_iter++;
    		copy_count++;
    	}
    }

    size_type num_dimensions() const { return NumDims; }

    size_type size() const { return extent_sz; }

    // given reshaping functionality, this is the max possible size.
    size_type max_size() const { return num_elements(); }

    bool empty() const { return size() == 0; }

    inline __device__ __host__ const index* strides() const {return stride_list_.data();}
    inline __device__ __host__ const element* origin() const { return base_; }
    inline __device__ __host__ const element* data() const { return base_; }

    size_type num_elements() const { return num_elements_; }

    const_iterator begin() const
    {
      return const_iterator(0,origin(),size(),strides());
    }

    const_iterator end() const
    {
      return const_iterator(size(),origin(),size(),strides());
    }

    typedef openfpm::array<index,NumDims> index_list;

    // This is used by multi_array, which is a subclass of this
    void set_base_ptr(TPtr new_base) { base_ = new_base; }


    TPtr base_;
    size_type extent_sz;
    size_type num_elements_;
    index_list stride_list_;

private:

    // const_multi_array_ref cannot be assigned to (no deep copies!)
    const_multi_array_ref_openfpm& operator=(const const_multi_array_ref_openfpm & other);

    void init_multi_array_ref(const index sz)
    {
      // calculate the extents
      extent_sz = sz;

      this->template compute_strides<StorageOrder>(stride_list_,extent_sz);

      num_elements_ = sz * size_ct::value;
    }
};


template <typename T, int NumDims, typename vector, typename StorageOrder>
class multi_array_ref_openfpm : public const_multi_array_ref_openfpm<T,NumDims,vector,T *,StorageOrder>
{
	typedef const_multi_array_ref_openfpm<T,NumDims,vector,T *,StorageOrder> super_type;
public:
/*  typedef typename super_type::value_type value_type;*/
  typedef typename super_type::reference reference;
  typedef typename super_type::iterator iterator;
/*  typedef typename super_type::reverse_iterator reverse_iterator;*/
  typedef typename super_type::const_reference const_reference;
  typedef typename super_type::const_iterator const_iterator;
/*  typedef typename super_type::const_reverse_iterator const_reverse_iterator;*/
  typedef typename super_type::element element;
  typedef typename super_type::size_type size_type;
//  typedef typename super_type::difference_type difference_type;
  typedef typename super_type::index index;
/*  typedef typename super_type::extent_range extent_range;*/

  typedef typename super_type::index_list index_list;

  //! indicate that this class is a multi dimensional array
  typedef int yes_is_multi_array;

/*  typedef typename super_type::size_list size_list;*/

  template <class ExtentType>
  explicit multi_array_ref_openfpm(T* base, const ExtentType r_sz)
  :super_type(base,r_sz)
  {
  }

  // Assignment from other ConstMultiArray types.
  template <typename ConstMultiArray>
  multi_array_ref_openfpm & operator=(const ConstMultiArray& other)
  {
    boost::function_requires<
      boost::multi_array_concepts::
      ConstMultiArrayConcept<ConstMultiArray,NumDims> >();

    // make sure the dimensions agree
    BOOST_ASSERT(other.num_dimensions() == this->num_dimensions());
    BOOST_ASSERT(std::equal(other.shape(),other.shape()+this->num_dimensions(),
                            this->shape()));
    // iterator-based copy
    std::copy(other.begin(),other.end(),this->begin());
    return *this;
  }

  multi_array_ref_openfpm & operator=(const multi_array_ref_openfpm & other)
  {
    if (&other != this)
    {
      // make sure the dimensions agree

      BOOST_ASSERT(other.num_dimensions() == this->num_dimensions());

      // iterator-based copy
      std::copy(other.begin(),other.end(),this->begin());
    }
    return *this;
  }

  multi_array_ref_openfpm & bind_ref(const multi_array_ref_openfpm & other)
  {
    if (&other != this) {

        this->base_ = other.base_;
        this->extent_sz = other.extent_sz;
        this->stride_list_ = other.stride_list_;
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

  /* \brief Get the internal pointer
   *
   * \return the internal pointer
   *
   */
  __device__ __host__ void * get_pointer()
  {
	  return this->base_;
  }

  /* \brief Get the internal pointer
   *
   * \return the internal pointer
   *
   */
  __device__ __host__ const void * get_pointer() const
  {
	  return this->base_;
  }

  multi_array_ref_openfpm & operator=(multi_array_ref_openfpm && other)
  {
	swap(other);

    return *this;
  }

  void swap(multi_array_ref_openfpm & other)
  {
    T* base_tmp = this->base_;
    this->base_ = other.base_;
    other.base_ = base_tmp;

    size_type extent_tmp = this->extent_sz;
    this->extent_sz = other.extent_sz;
    other.extent_sz = extent_tmp;

    index_list stride_list_tmp = this->stride_list_;
    this->stride_list_ = other.stride_list_;
    other.stride_list_ = stride_list_tmp;

    size_type num_elements_tmp = this->num_elements_;
    this->num_elements_ = other.num_elements_;
    other.num_elements_ = num_elements_tmp;
  }

  __device__ __host__ element* origin() { return super_type::base_; }

  __device__ __host__ const element* origin() const { return super_type::origin(); }

  __device__ __host__ reference operator[](index idx)
  {
    return super_type::access(boost::type<reference>(),
                              idx,
                              this->strides(),
                              this->origin());
  }



  iterator begin()
  {return iterator(0,origin(),this->size(),this->strides());}

  iterator end()
  {return iterator(this->size(),origin(),this->size(),this->strides());}

  __inline__ __device__ __host__ const_reference operator[](index idx) const
  {
	  return super_type::access(boost::type<const_reference>(),
                              idx,
                              this->strides(),
                              this->origin());
  }

  const_iterator begin() const
  {return super_type::begin();}

  const_iterator end() const
  {return super_type::end();}
};

template<typename T, typename Sfinae = void>
struct is_multi_array: std::false_type {};


/*! \brief has_noPointers check if a type has defined a
 * method called noPointers
 *
 * ### Example
 *
 * \snippet util_test.hpp Check no pointers
 *
 * return true if T::noPointers() is a valid expression (function pointers)
 * and produce a defined type
 *
 */
template<typename T>
struct is_multi_array<T, typename Void<typename T::yes_is_multi_array >::type> : std::true_type
{};

} // namespace openfpm


#endif /* MULTI_ARRAY_REF_OPENFPM_HPP_ */
