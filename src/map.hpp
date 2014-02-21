
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>

#include "grid.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "memory_cpu.hpp"
#include "memory_gpu.hpp"

#ifndef MAP_HPP_
#define MAP_HPP_

template<typename T, typename Mem> class layout_cpu
{
};

template<typename T, typename Mem> class layout_gpu
{
};

// Compiled based specialization

template<typename T, typename Mem>
class layout_cpu<grid<T>,Mem>
{
  grid<T> g1;
  Mem data;
  
public:
  
  layout_cpu(std::vector<size_t> & sz)
  :g1(sz)
  {
    data.allocate(g1.size());
  }
  
  template <unsigned int p>inline typename boost::fusion::result_of::at<Point<float>::type,boost::mpl::int_<p> >::type & get(grid_key<p> & v1)
  {
    return boost::fusion::at_c<p>(data.get(g1.LinId(v1.getId())));
  }
  
  template <unsigned int p>inline typename boost::fusion::result_of::at<Point<float>::type,boost::mpl::int_<p> >::type & get(grid_key_1<p> & v1)
  {
    return boost::fusion::at_c<p>(data.get(g1.LinId(v1.k[0])));
  }
  
  template <unsigned int p>inline typename boost::fusion::result_of::at<Point<float>::type,boost::mpl::int_<p> >::type & get(grid_key_2<p> & v1)
  {
    return boost::fusion::at_c<p>(data.get(g1.LinId(v1.k[1],v1.k[0])));
  }
  
  template <unsigned int p>inline typename boost::fusion::result_of::at<Point<float>::type,boost::mpl::int_<p> >::type & get(grid_key_3<p> & v1)
  {
    return boost::fusion::at_c<p>(data.get(g1.LinId(v1.k[2],v1.k[1],v1.k[0])));
  }
  
  template <unsigned int p>inline typename boost::fusion::result_of::at<Point<float>::type,boost::mpl::int_<p> >::type & get(grid_key_4<p> & v1)
  {
    return boost::fusion::at_c<p>(data.get(g1.LinId(v1.k[3],v1.k[2],v1.k[1],v1.k[0])));
  }
  
  template <unsigned int p, unsigned int dim>inline typename boost::fusion::result_of::at<Point<float>::type,boost::mpl::int_<p> >::type & get(grid_key_d<dim,p> & v1)
  {
    return boost::fusion::at_c<p>(data.get(g1.LinId(v1)));
  }
  
  
  template <unsigned int p, unsigned int dim>inline typename boost::fusion::result_of::at<Point<float>::type,boost::mpl::int_<p> >::type & get(grid_key_dx<dim> & v1)
  {
    return boost::fusion::at_c<p>(data.get(g1.LinId(v1)));
  }
  
  inline size_t size()
  {
    return g1.size();
  }
};

/*template<unsigned int dim, typename T>
class layout_cpu<Particles<dim,Point<T>,memory_cpu<T>>,part_key>
{
  grid<Point<T>,memory_cpu<Point<T>>> * g1;
  
public:
  
  layout_cpu(std::vector<size_t> sz)
  {
    g1 = new grid<Point<T>,memory_cpu<Point<T>>>(sz);
  }
  
  inline Point<T> & get(grid_key & v1)
  {
    return g1->get(v1);
  }
  
  inline size_t size()
  {
    return g1->size();
  }
};*/

/*template<typename T>
class layout_cpu<grid<tensor<T>,memory_cpu<T>>,grid_key>
{
  grid<Point<T>,memory_cpu<T>> & g1;
  
  layout_cpu(grid<Point<T>,memory_cpu<T>> & g1)
  : g1(g1)
  {
  }
  
  inline Point<T> get(grid_key & v1)
  {
    return g1.get(v1);
  }
  
  inline size_t size()
  {
    g1.size();
  }
};*/

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

template<unsigned int p>
struct Point_type_prop
{
	  typedef typename boost::fusion::result_of::at<Point<float>::type,boost::mpl::int_<p> >::type type;
};

template<typename T, typename Mem>
class layout_gpu<grid<Point<T>>,Mem>
{
  grid<T> g1;

  typename Mem::type data;

public:

  layout_gpu(std::vector<size_t> & sz)
  :g1(sz)
  {
	  boost::fusion::at_c<0>(data).allocate(g1.size());
	  boost::fusion::at_c<1>(data).allocate(g1.size());
	  boost::fusion::at_c<2>(data).allocate(g1.size());
	  boost::fusion::at_c<3>(data).allocate(g1.size());
	  boost::fusion::at_c<4>(data).allocate(g1.size());
	  boost::fusion::at_c<5>(data).allocate(g1.size());
  }


  template <unsigned int p>inline typename mem_reference<typename Point_type_prop<p>::type>::type get(grid_key<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.getId()));
  }

  template <unsigned int p>inline typename mem_reference<typename Point_type_prop<p>::type >::type get(grid_key_1<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.k[0]));
  }

  template <unsigned int p>inline typename mem_reference<typename Point_type_prop<p>::type >::type get(grid_key_2<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.k[1],v1.k[0]));
  }

  template <unsigned int p>inline typename mem_reference<typename Point_type_prop<p>::type >::type get(grid_key_3<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.k[2],v1.k[1],v1.k[0]));
  }

  template <unsigned int p>inline typename mem_reference<typename Point_type_prop<p>::type >::type get(grid_key_4<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.k[3],v1.k[2],v1.k[1],v1.k[0]));
  }

  template <unsigned int p, unsigned int dim>inline typename mem_reference<typename Point_type_prop<p>::type >::type get(grid_key_d<dim,p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1));
  }

  template <unsigned int p, unsigned int dim>inline typename mem_reference<typename Point_type_prop<p>::type >::type get(grid_key_dx<dim> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1));
  }

  inline size_t size()
  {
    return g1.size();
  }
};

#endif


