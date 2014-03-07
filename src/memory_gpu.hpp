#ifndef MEMORY_GPU_HPP_
#define MEMORY_GPU_HPP_

#include "config.h"
#include <algorithm>
#include <cstdlib>
#include <boost/shared_array.hpp>
#include "base_type.hpp"

typedef long int mem_id;

#ifndef NVCC

template<typename T>
struct mem_reference<T>
{
	typedef T& type;
};

template<typename T>
class memory_gpu
{
	  boost::shared_array<T> mem;
	  // Thrust memory management

	  size_t sz;

	public:

	  typedef T type;

	  T & get(mem_id id)
	  {
	    return mem.get()[id];
	  }

	  inline void set(mem_id id, T m)
	  {
	    mem.get()[id] = m;
	  }

	  inline bool allocate(size_t sz)
	  {
	    mem = boost::shared_array<T>(new T[sz]);
	    return true;
	  }

	  inline void destroy()
	  {
	    mem.reset();
	  }

	  inline bool copy(memory_gpu m)
	  {
	    memcpy(mem.get(),m.mem.get(),m.size());
	    return true;
	  }

	  inline size_t size()
	  {
	    return sz;
	  }
};

#endif

#ifdef NVCC

#include <memory_gpu_thrust.hpp>
#include <boost/type_traits/remove_reference.hpp>

template<typename T>
struct mem_reference
{
	typedef typename boost::remove_reference<T>::type type;
};

template<> struct mem_reference<float&>
{
	typedef float& type;
};

template<> struct mem_reference<int&>
{
	typedef int& type;
};

template<> struct mem_reference<double&>
{
	typedef double& type;
};

template<> struct mem_reference<float>
{
	typedef float& type;
};

template<> struct mem_reference<int>
{
	typedef int& type;
};

template<> struct mem_reference<double>
{
	typedef double& type;
};

template<typename T> struct mem_reference<boost::multi_array_ref<T,1>>
{
	typedef typename boost::multi_array_ref<T,1>::reference type;
};

template<typename T> struct mem_reference<boost::multi_array_ref<T,2>>
{
	typedef typename boost::multi_array_ref<T,2>::reference type;
};

template<typename T> struct mem_reference<boost::multi_array_ref<T,3>>
{
	typedef typename boost::multi_array_ref<T,3>::reference type;
};

template<typename T> struct mem_reference<boost::multi_array_ref<T,1>&>
{
	typedef typename boost::multi_array_ref<T,1>::reference type;
};

template<typename T> struct mem_reference<boost::multi_array_ref<T,2>&>
{
	typedef typename boost::multi_array_ref<T,2>::reference type;
};

template<typename T> struct mem_reference<boost::multi_array_ref<T,3>&>
{
	typedef typename boost::multi_array_ref<T,3>::reference type;
};

template<typename T>
class memory_gpu
{
	  memory_gpu_thrust< typename base_type<T>::type > mem;
	  T * mem_ptr;
	  size_t sz;

	public:

	  typedef T type;

	  typename mem_reference<T>::type get(mem_id id)
	  {
	    return mem_ptr[id];
	  }

	  inline void set(mem_id id, T m)
	  {
	    mem_ptr[id] = m;
	  }

	  inline bool allocate(size_t sz)
	  {
	    mem.resize(sz);
	    mem_ptr = mem.getPointer();
	    return true;
	  }

	  inline void destroy()
	  {
	  }

	  inline bool copy(memory_gpu m)
	  {
		memory_gpu_thrust<T>::copy(mem.dv,m.mem.dv);
	    return true;
	  }

	  inline size_t size()
	  {
	    return mem->dv.size();
	  }
};

template<typename T, unsigned int i, unsigned int j>
class memory_gpu_array2D
{
	  memory_gpu_thrust<T> mem;
	  boost::multi_array_ref<T,3> * ma;
	  size_t sz;

	public:

	  typedef T type;

	  memory_gpu_array2D()
	  : ma(NULL), sz(0)
	  {}

	  typename boost::multi_array_ref<T,3>::reference get(mem_id id)
	  {
	    return (*ma)[id];
	  }

	  inline void set(mem_id id, T m)
	  {
	    mem.dv[id] = m;
	  }

	  inline bool allocate(size_t sz)
	  {
	    mem.resize(sz);

	    // shape on ma

	    bool ascending[] = {true,true,true};
	    typename boost::multi_array_ref<T,3>::size_type xDim(i), yDim(j), zDim(sz);

	    typename boost::multi_array<T,3>::size_type ordering[] = {2,1,0};
	    ma = new boost::multi_array_ref<T,3>(mem.getPointer(),boost::extents[xDim][yDim][zDim],boost::general_storage_order<3>(ordering,ascending));

	    return true;
	  }

	  inline void destroy()
	  {
	  }

	  inline bool copy(memory_gpu_array2D m)
	  {
		memory_gpu_thrust<T>::copy(mem.dv,m.mem.dv);
	    return true;
	  }

	  inline size_t size()
	  {
	    return mem.dv.size();
	  }
};

template<typename T, unsigned int i>
class memory_gpu_array1D
{
	  memory_gpu_thrust<T> mem;
	  boost::multi_array_ref<T,2> * ma;
	  size_t sz;

	public:

	  typedef T type;

	  memory_gpu_array1D()
	  : ma(NULL), sz(0)
	  {}

	  typename boost::multi_array_ref<T,2>::reference get(mem_id id)
	  {
	    return (*ma)[id];
	  }

	  inline void set(mem_id id, T m)
	  {
	    mem.dv[id] = m;
	  }

	  inline bool allocate(size_t sz)
	  {
	    mem.resize(sz);

	    // shape array

	    bool ascending[] = {true,true};

	    typename boost::multi_array_ref<T,2>::size_type xDim(i), zDim(sz);

	    typename boost::multi_array<T,2>::size_type ordering[] = {1,0};
	    ma = new boost::multi_array_ref<T,2>(mem.getPointer(),boost::extents[xDim][zDim],boost::general_storage_order<2>(ordering,ascending));

	    return true;
	  }

	  inline void destroy()
	  {
	  }

	  inline bool copy(memory_gpu_array1D m)
	  {
		memory_gpu_thrust<T>::copy(mem.dv,m.mem.dv);
	    return true;
	  }

	  inline size_t size()
	  {
	    return mem.dv.size();
	  }
};

#endif


////////////////////////////// compile based specialization

#include "Point.hpp"

template<typename T>
class memory_gpu_type
{
	typedef T type;
};


template<typename T>
struct memory_gpu_type< Point<T> >
{
	typedef typename Point<T>::type ptype;

	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<0> >::type ptype_0;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<1> >::type ptype_1;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<2> >::type ptype_2;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<3> >::type ptype_3;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<4> >::type ptype_4;
	typedef typename boost::fusion::result_of::at<ptype,boost::mpl::int_<5> >::type ptype_5;

	typedef typename boost::remove_reference<ptype_0>::type mt_0;
	typedef typename boost::remove_reference<ptype_1>::type mt_1;
	typedef typename boost::remove_reference<ptype_2>::type mt_2;
	typedef typename boost::remove_reference<ptype_3>::type mt_3;
	typedef typename boost::remove_reference<ptype_4>::type mt_4;
	typedef typename boost::remove_reference<ptype_5>::type mt_5;

	typedef typename mt_4::element mt_4ele;
	typedef typename mt_5::element mt_5ele;

	typedef boost::fusion::vector< memory_gpu<mt_0>,memory_gpu<mt_1>,memory_gpu<mt_2>,memory_gpu<mt_3>,memory_gpu_array1D<mt_4ele,3>,memory_gpu_array2D<mt_5ele,3,3> > type;
};

#endif
