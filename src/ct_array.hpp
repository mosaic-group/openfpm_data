/*
 * ct_array.hpp
 *
 *  Created on: Aug 21, 2014
 *      Author: Pietro Incardona
 */



/*! \brief These set of classes generate an array definition at compile-time
 *
 * These set of classes generate an array definition at compile-time
 *
 * \see generate_array
 *
 */

#ifndef CT_ARRAY_HPP_
#define CT_ARRAY_HPP_

#include <boost/fusion/mpl.hpp>

//! the array itself
template<class T, unsigned long... args> struct ArrayHolder {
    static const T data[sizeof...(args)];
};

//! initialize the array from variadic template
template<class T, unsigned long... args>
const T ArrayHolder<T,args...>::data[sizeof...(args)] = { args... };

//! Generate the array specializing ArrayHolder
template<class T,size_t N, size_t orig_N, template<size_t,size_t> class F, unsigned... args>
struct generate_array_impl {
    typedef typename generate_array_impl<T,N-1,orig_N, F, F<N,orig_N>::value, args...>::result result;
};

//! terminator of the variadic template
template<class T, size_t orig_N, template<size_t,size_t> class F, unsigned... args>
struct generate_array_impl<T,0,orig_N, F, args...> {
    typedef ArrayHolder<T,F<0,orig_N>::value, args...> result;
};

/*! \brief Main class to generate array
 *
 * Main class to generate array
 *
 * Usage:
 *
 * template<size_t index, size_t N> struct MetaFunc {
 *   enum { value = index + N };
 * };
 *
 * const size_t count = 5;
 * typedef generate_array<int,count, MetaFunc>::result A;
 *
 *  for (size_t i=0; i<count; ++i)
 *      std::cout << A::data[i] << "\n";
 * }
 *
 * A::data is an array of int that go from 1 to 5
 *
 * \param T is the type ot the output array
 * \param N size of the sequence
 * \param F Meta function it take two template arguments
 *
 */
template<class T, size_t N, template<size_t,size_t> class F>
struct generate_array {
    typedef typename generate_array_impl<T,N-1, N, F>::result result;
};

/////////////// Classes to generate at compile time arrays from a boost::mpl::vector

//! Generate the array specializing ArrayHolder
template<class T, size_t N ,class F, unsigned... args>
struct generate_array_vector_impl {
    typedef typename generate_array_vector_impl<T,N-1,F, boost::mpl::at<F,boost::mpl::int_<N> >::type::value , args...>::result result;
};

//! terminator of the variadic template
template<class T, class F, unsigned... args>
struct generate_array_vector_impl<T,1, F, args...> {
    typedef ArrayHolder<T,boost::mpl::at<F,boost::mpl::int_<1> >::type::value, args...> result;
};

/*! \brief Main class to generate an array from a boost::mpl::vector of numbers
 *
 * Main class to generate an array from a boost::mpl::vector of numbers
 *
 * Usage:
 *
 * boost::mpl::vector<float,int_<4>,int<4>,......> v
 *
 * typdef generate_array_vector<size_t,v>::result B;
 *
 * B is an array or size_t {4,4}
 *
 * \param T type of the output buffer
 * \param F boost::mpl::vector
 *
 */
template<class T, class F>
struct generate_array_vector {
    typedef typename generate_array_vector_impl<T,boost::mpl::size<F>::value-1 , F>::result result;
};

#endif
