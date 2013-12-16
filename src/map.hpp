
#include "grid.hpp"
#include "Point.hpp"
#include "memory.hpp"

#ifndef MAP_HPP_
#define MAP_HPP_

template<typename T, typename K> class layout_cpu
{
  mem_id get(K i)
  {
    return i.toId();
  }
};

template<typename T, typename K> class layout_gpu
{
  mem_id get(K i)
  {
    return i.toId();
  }
};

// Compiled based specialization

template<typename T>
class layout_cpu<grid<Point<T>,memory_cpu<T>>,grid_key>
{
  grid<Point<T>,memory_cpu<T>> & g1;
  
  layout_cpu(grid<Point<T>,memory_cpu<T>> & g1)
  : g1(g1)
  {
  }
  
  inline mem_id get(grid_key & v1)
  {
    return g1.idx(v1)*sizeof(Point<T>);
  }
  
  inline size_t size()
  {
    g1.size();
  }
};

/*template<>
map_cpu<grid<grid<Point>>>
{
  grid & g1;
  grid & g2;
  
  map(grid & g1, grid & g2)
  : g1(g1), g2(g2)
  {
  }
  
  inline map_id <unsigned int p>get(key & v1, key & v2)
  {
    return g1.idx(v1)*g2.size()+g2.idx(v2)*sizeof(Point);
  }
  
  inline map_id <>
  {
  }
}*/

template<typename T>
class layout_gpu<grid<Point<T>,memory_gpu<T>>,grid_key>
{
  grid<T,memory_gpu<T>> x;
  grid<T,memory_gpu<T>> y;
  grid<T,memory_gpu<T>> z;
  
  layout_gpu(grid<Point<T>,memory_gpu<T>> & g1)
  :x(g1),y(g1),z(g1)
  {
  }
  
  inline mem_id get_x(grid_key & v1)
  {
    return x.idx(v1)*sizeof(Point<T>);
  }
  
  inline mem_id get_y(grid_key & v1)
  {
    return x.idx(v1)*sizeof(Point<T>);
  }
  
  inline mem_id get_z(grid_key & v1)
  {
    return x.idx(v1)*sizeof(Point<T>);
  }
  
  inline size_t size()
  {
    return 3*x.size();
  }
};

#endif


