#ifndef MEMORY_CPU_HPP_
#define MEMORY_CPU_HPP_

#include <boost/shared_array.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include <base_type.hpp>

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

////////////////////////////// compile based specialization

#include "Point.hpp"

template<typename T>
struct memory_cpu_type
{
	typedef typename T::type type;
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
