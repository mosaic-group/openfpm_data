#ifndef BASE_TYPE_HPP
#define BASE_TYPE_HPP

#include <boost/multi_array.hpp>

template<class T>
struct base_type
{
	typedef T type;
};


template<class T>
struct base_type<boost::multi_array_ref<T,1>>
{
	typedef T type;
};

template<class T>
struct base_type<boost::multi_array_ref<T,2>>
{
	typedef T type;
};

template<class T>
struct base_type<boost::multi_array_ref<T,3>>
{
	typedef T type;
};

template<class T>
struct base_type<boost::multi_array_ref<T,4>>
{
	typedef T type;
};

template<class T>
struct base_type<boost::multi_array_ref<T,5>>
{
	typedef T type;
};

#endif
