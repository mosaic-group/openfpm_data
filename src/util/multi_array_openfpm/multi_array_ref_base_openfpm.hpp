/*
 * multi_array_ref_openfpm_base.hpp
 *
 *  Created on: Jun 30, 2018
 *      Author: i-bird
 */

#ifndef MULTI_ARRAY_REF_OPENFPM_BASE_HPP_
#define MULTI_ARRAY_REF_OPENFPM_BASE_HPP_

#include "types.hpp"
#include <boost/mpl/size_t.hpp>
#include "boost/iterator/reverse_iterator.hpp"
#include "storage_order.hpp"
#include <boost/mpl/at.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/pop_front.hpp>
#include <boost/type.hpp>

namespace openfpm
{

/////////////////////////////////////////////////////////////////////////
// class declarations
/////////////////////////////////////////////////////////////////////////

//template<typename T, std::size_t NumDims, typename Allocator = std::allocator<T> >
//class multi_array_openfpm;

template <typename T, int NumDims, typename vector>
class multi_array_ref_openfpm;

// This is a public interface for use by end users!
namespace multi_array_types
{
	typedef openfpm::detail::multi_array::size_type size_type;
	typedef std::ptrdiff_t difference_type;
	typedef openfpm::detail::multi_array::index index;
}


namespace detail {
namespace multi_array {

template <typename T, std::size_t NumDims, typename vector>
class sub_array_openfpm;

template <typename T, std::size_t NumDims, typename vector, typename TPtr = const T*>
class const_sub_array_openfpm;

  template <typename T, typename TPtr, typename NumDims, typename vector, typename Reference,
            typename IteratorCategory>
class array_iterator_openfpm;

template <typename T, std::size_t NumDims, typename TPtr = const T*>
class const_multi_array_view_openfpm;

template <typename T, std::size_t NumDims>
class multi_array_view_openfpm;


/////////////////////////////////////////////////////////////////////////
// class interfaces
/////////////////////////////////////////////////////////////////////////

class multi_array_base_openfpm
{
public:
	typedef multi_array_types::size_type size_type;
	typedef multi_array_types::difference_type difference_type;
	typedef multi_array_types::index index;
};

//
// value_accessor_n
//  contains the routines for accessing elements from
//  N-dimensional views.
//
template<typename T, std::size_t NumDims, typename vector>
class value_accessor_n_openfpm : public multi_array_base_openfpm
{
	typedef multi_array_base_openfpm super_type;
public:
	typedef typename super_type::index index;

	//
	// public typedefs used by classes that inherit from this base
	//
	typedef T element;
	typedef openfpm::multi_array_ref_openfpm<T,NumDims-1,typename boost::mpl::pop_front<vector>::type> value_type;
	typedef sub_array_openfpm<T,NumDims-1,typename boost::mpl::pop_front<vector>::type> reference;
	typedef const_sub_array_openfpm<T,NumDims-1,typename boost::mpl::pop_front<vector>::type> const_reference;

protected:

	// used by array operator[] and iterators to get reference types.
	template <typename Reference, typename TPtr>
	__device__ __host__ inline Reference access(boost::type<Reference>,
										 index idx,
										 const index* strides,
										 TPtr base) const
	{
		TPtr newbase = base + idx * strides[0];
		return Reference(newbase,strides+1);
	}

	__device__ __host__ value_accessor_n_openfpm() { }
	__device__ __host__ ~value_accessor_n_openfpm() { }
};

template <class T> inline __device__ __host__ void ignore_unused_variable_warning_ofp(T const&) {}

//
// value_accessor_one
//  contains the routines for accessing reference elements from
//  1-dimensional views.
//
template<typename T, typename vector>
class value_accessor_one_openfpm : public multi_array_base_openfpm
{
	typedef multi_array_base_openfpm super_type;
public:
	typedef typename super_type::index index;
	//
	// public typedefs for use by classes that inherit it.
	//
	typedef T element;
	typedef T value_type;
	typedef T& reference;
	typedef T const& const_reference;

protected:

	// used by array operator[] and iterators to get reference types.
	template <typename Reference, typename TPtr>
	inline __device__ __host__ Reference access(boost::type<Reference>,index idx,
                   const index* strides,
                   TPtr base) const
	{
		return *(base + idx * strides[0]);
	}

	// used by array operator[] and iterators to get reference types.
	template <typename Reference, typename TPtr>
	inline __device__ __host__ Reference access(boost::type<Reference>,index idx,TPtr base, const index* strides) const
	{
		BOOST_ASSERT(size_type(idx < boost::mpl::at<vector,boost::mpl::int_<0>>::type::value));
		return *(base + idx * strides[0]);
	}

	__device__ __host__ value_accessor_one_openfpm() { }
	__device__ __host__ ~value_accessor_one_openfpm() { }
};


/////////////////////////////////////////////////////////////////////////
// choose value accessor begins
//

template <typename T, std::size_t NumDims,typename vector>
struct choose_value_accessor_n_openfpm
{
	typedef value_accessor_n_openfpm<T,NumDims,vector> type;
};

template <typename T,typename vector>
struct choose_value_accessor_one_openfpm
{
	typedef value_accessor_one_openfpm<T,vector> type;
};

template <typename T, typename NumDims, typename vector>
struct value_accessor_generator_openfpm
{
	BOOST_STATIC_CONSTANT(std::size_t, dimensionality = NumDims::value);

	typedef typename
			boost::mpl::eval_if_c<(dimensionality == 1),
                  choose_value_accessor_one_openfpm<T,vector>,
                  choose_value_accessor_n_openfpm<T,dimensionality,vector>
			>::type type;
};


template <class T, class NumDims, typename vector>
struct associated_types_openfpm: value_accessor_generator_openfpm<T,NumDims,vector>::type
{};

//
// choose value accessor ends
/////////////////////////////////////////////////////////////////////////

// Due to some imprecision in the C++ Standard,
// MSVC 2010 is broken in debug mode: it requires
// that an Output Iterator have output_iterator_tag in its iterator_category if
// that iterator is not bidirectional_iterator or random_access_iterator.
#if BOOST_WORKAROUND(BOOST_MSVC, >= 1600)
struct mutable_iterator_tag
 : boost::random_access_traversal_tag, std::input_iterator_tag
{
  operator std::output_iterator_tag() const {
    return std::output_iterator_tag();
  }
};
#endif


////////////////////////////////////////////////////////////////////////
// multi_array_base
////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t NumDims, typename vector>
class multi_array_impl_base_openfpm: public value_accessor_generator_openfpm<T,boost::mpl::size_t<NumDims>,vector>::type
{
	typedef associated_types_openfpm<T,boost::mpl::size_t<NumDims>,vector > types;
public:

	typedef typename types::index index;
	typedef typename types::size_type size_type;
	typedef typename types::element element;
	typedef typename types::value_type value_type;
	typedef typename types::reference reference;
	typedef typename types::const_reference const_reference;

	template <std::size_t NDims>
	struct subarray
	{
		typedef openfpm::detail::multi_array::sub_array_openfpm<T,NDims,vector> type;
	};

	template <std::size_t NDims>
	struct const_subarray
	{
		typedef openfpm::detail::multi_array::const_sub_array_openfpm<T,NDims,vector> type;
	};

	template <std::size_t NDims>
	struct array_view_openfpm
	{
		typedef openfpm::detail::multi_array::multi_array_view_openfpm<T,NDims> type;
	};

	template <std::size_t NDims>
	struct const_array_view_openfpm
	{
		public:
			typedef openfpm::detail::multi_array::const_multi_array_view_openfpm<T,NDims> type;
	};

#if BOOST_WORKAROUND(BOOST_MSVC, >= 1600)
	// Deal with VC 2010 output_iterator_tag requirement
	typedef array_iterator_openfpm<T,T*,mpl::size_t<NumDims>,vector,reference,
                         mutable_iterator_tag> iterator;
#else
	typedef array_iterator_openfpm<T,T*,boost::mpl::size_t<NumDims>, vector,reference,
                         boost::random_access_traversal_tag> iterator;
#endif
	typedef array_iterator_openfpm<T,T const*,boost::mpl::size_t<NumDims>,vector,const_reference,
                                   boost::random_access_traversal_tag> const_iterator;

	typedef ::boost::reverse_iterator<iterator> reverse_iterator;
	typedef ::boost::reverse_iterator<const_iterator> const_reverse_iterator;

	BOOST_STATIC_CONSTANT(std::size_t, dimensionality = NumDims);

protected:

	__device__ __host__ multi_array_impl_base_openfpm() { }
	__device__ __host__ ~multi_array_impl_base_openfpm() { }


	  template <typename Stride_list, typename Extent_type>
	  void compute_strides(Stride_list& stride_list, Extent_type& extent,
			               const general_storage_order<NumDims>& storage)
	  {
		  // invariant: stride = the stride for dimension n
		  index stride = 1;
		  for (size_type n = 0; n != NumDims; ++n)
		  {
			  // The stride for this dimension is the product of the
			  // lengths of the ranks minor to it.
			  stride_list[storage.ordering(n)] = stride;

			  if (storage.ordering(n) == 0)
			  {stride *= extent;}
			  else
			  {
				  switch(storage.ordering(n))
				  {
				  case 1:
					  stride *= at_impl<vector,1>::type::value;
					  break;
				  case 2:
					  stride *= at_impl<vector,2>::type::value;
					  break;
				  case 3:
					  stride *= at_impl<vector,3>::type::value;
					  break;
				  case 4:
					  stride *= at_impl<vector,4>::type::value;
					  break;
				  case 5:
					  stride *= at_impl<vector,5>::type::value;
					  break;
				  case 6:
					  stride *= at_impl<vector,6>::type::value;
					  break;
				  case 7:
					  stride *= at_impl<vector,7>::type::value;
					  break;
				  case 8:
					  stride *= at_impl<vector,8>::type::value;
					  break;
				  case 9:
					  stride *= at_impl<vector,9>::type::value;
					  break;
				  case 10:
					  stride *= at_impl<vector,10>::type::value;
					  break;
				  case 11:
					  stride *= at_impl<vector,11>::type::value;
					  break;
				  case 12:
					  stride *= at_impl<vector,12>::type::value;
					  break;
				  case 13:
					  stride *= at_impl<vector,13>::type::value;
					  break;
				  case 14:
					  stride *= at_impl<vector,14>::type::value;
					  break;
				  case 15:
					  stride *= at_impl<vector,15>::type::value;
					  break;
				  }
			  }
		  }
	  }

	// Used by operator() in our array classes
	template <typename Reference, typename IndexList, typename TPtr>
	Reference access_element(boost::type<Reference>,
                           const IndexList& indices,
                           TPtr base,
                           const size_type* extents,
                           const index* strides,
                           const index* index_bases) const
	{
		boost::function_requires<
			boost::CollectionConcept<IndexList> >();
		ignore_unused_variable_warning(index_bases);
		ignore_unused_variable_warning(extents);
#ifdef SE_CLASS1
		for (size_type i = 0; i != NumDims; ++i)
		{
			BOOST_ASSERT(indices[i] - index_bases[i] >= 0);
			BOOST_ASSERT(size_type(indices[i] - index_bases[i]) < extents[i]);
		}
#endif

		index offset = 0;
		{
			typename IndexList::const_iterator i = indices.begin();
			size_type n = 0;
			while (n != NumDims)
			{
				offset += (*i) * strides[n];
				++n;
				++i;
			}
		}
		return base[offset];
	}
};

} // namespace multi_array
} // namespace detail

} // namespace openfpm


#endif /* MULTI_ARRAY_REF_OPENFPM_BASE_HPP_ */
