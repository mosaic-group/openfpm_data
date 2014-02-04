#ifndef MEMORY_HPP_
#define MEMORY_HPP_

#include <boost/shared_array.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>

typedef long int mem_id;

template<typename T>
class memory_cpu
{
  boost::shared_array<T> mem;
  size_t sz;
  
public:
  
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
  
  inline bool copy(memory_cpu m)
  {
    memcpy(mem.get(),m.mem.get(),m.size());
    return true;
  }
  
  inline size_t size()
  {
    return sz;
  }
};


template<typename T>
class memory_gpu
{
	  boost::shared_array<T> mem;
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

////////////////////////////// compile based specialization

#include "Point.hpp"

template<typename T>
class memory_gpu_type
{
	typedef T type;
};


template<typename T>
struct memory_gpu_type<Point<T>>
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

	typedef boost::fusion::vector<memory_gpu<mt_0>,memory_gpu<mt_1>,memory_gpu<mt_2>,memory_gpu<mt_3>,memory_gpu<mt_4>,memory_gpu<mt_5>> type;
};

template<typename T>
struct memory_cpu_type
{
	typedef T type;
};

template<typename T>
struct memory_cpu_type<Point<T>>
{
	typedef typename Point<T>::type type;
};

////////////////// Compile based specialization

struct memory_null
{  
};


#endif
