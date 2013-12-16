#include <iostream>
#include <boost/shared_array.hpp>
#include <vector>
#include <memory.hpp>

#define HARDWARE 1

template<typename T, unsigned int s1, unsigned int s2, unsigned int s3>
class tensor
{
  T g_mem[s1][s2][s3];
 
public:
  tensor() {};
  ~tensor() {};
  
  size_t size()
  {
    return s1*s2*s3*g_mem[0][0][0].size();
  };
};

template<typename T, typename Mem>
class grid
{
  size_t size_tot;
  std::vector<size_t> sz;
  Mem g_mem;
 
public:
  
  // Static element to calculate total size

  size_t totalSize(std::vector<size_t> sz)
  {
    size_t tSz = 1;
    
    for (int i = 0 ;  i < sz.size() ; i++)
    {
      tSz *= sz[i];
    }
    
    return tSz;
  }
  
  grid(std::vector<size_t> sz) 
  : size_tot(totalSize(sz)), sz(sz)
  {
    g_mem.create(totalSize(sz));
  };
  
  ~grid() {};
  
  size_t size()
  {
    return size_tot;
  };
};

template<typename T>
class grid<T,memory_null>
{
  size_t size_tot;
  std::vector<size_t> sz;
 
public:
  
  // Static element to calculate total size

  size_t totalSize(std::vector<size_t> sz)
  {
    size_t tSz = 1;
    
    for (int i = 0 ;  i < sz.size() ; i++)
    {
      tSz *= sz[i];
    }
    
    return tSz;
  }
  
  grid(std::vector<size_t> sz) 
  : size_tot(totalSize(sz)), sz(sz)
  {
  };
  
  ~grid() {};
  
  size_t size()
  {
    return size_tot;
  };
};

template<unsigned int s1, unsigned int s2, unsigned int s3>
class tensor<int,s1,s2,s3>
{
  
public:
  size_t size()
  {
    return s1*s2*s3;
  }
};




class grid_key
{
  mem_id id[];
};

