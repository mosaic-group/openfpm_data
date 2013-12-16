#include <boost/shared_array.hpp>

#ifndef MEMORY_HPP_
#define MEMORY_HPP_

typedef long int mem_id;

template<typename T>
class memory_cpu
{
  boost::shared_array<T> mem;
  size_t sz;
  
  T get(mem_id id)
  {
    return mem.get()[id];
  }
  
  inline void set(mem_id id, T m)
  {
    mem.get()[id] = m;
  }
  
  inline bool create(size_t sz)
  {
    mem = boost::shared_array<T>(new T[sz]);
  }
  
  inline void destroy()
  {
    mem.reset();
  }
  
  inline bool copy(memory_cpu m)
  {
    memcpy(mem.get(),m.mem.get(),m.size());
  }
  
  inline size_t size()
  {
    return sz;
  }
};


template<typename T>
class memory_gpu
{
};

class memory_null
{  
};


#endif
