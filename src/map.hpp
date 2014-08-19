
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>

#include "grid.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "memory_cpu.hpp"
#include "memory_gpu.hpp"
#include "memory_array.hpp"

#ifndef MAP_HPP_
#define MAP_HPP_


// Compiled based specialization

template<unsigned int p>
struct Point_type_cpu_prop
{
	  typedef typename boost::fusion::result_of::at<memory_cpu_type<Point<float>>::rtype,boost::mpl::int_<p> >::type type;
};

/*!
 *
 * \brief This is an N-dimensional grid or an N-dimensional array working on CPU
 *
 * This is an N-Dimensional grid or an N-dimensional array working on CPU
 *
 *	\param dim Dimensionality of the grid
 *	\param T type of object the grid store
 *	\param Mem interface used to allocate memory
 *
 */

template<unsigned int dim, typename T, typename Mem>
class grid_cpu
{
		//! This is an header that store all information related to the grid
		grid<dim,T> g1;

		//! This is the interface to allocate an resize memory
		//! and give also a representation to the allocated memory
		Mem data;

	public:
  
		//! Constructor allocate memory and give them a representation
		grid_cpu(std::vector<size_t> & sz)
		:g1(sz)
		{
			//! allocate the memory
			data.mem.allocate(g1.size()*sizeof(T));

			//! set the pointer
			data.mem_r.set_pointer(data.getPointer());
		}
  
  template <unsigned int p>inline typename Point_type_cpu_prop<p>::type & get(grid_key<p> & v1)
  {
    return boost::fusion::at_c<p>(data.mem_r.get(g1.LinId(v1.getId())));
  }
  
  template <unsigned int p>inline typename Point_type_cpu_prop<p>::type & get(grid_key_1<p> & v1)
  {
    return boost::fusion::at_c<p>(data.mem_r.get(g1.LinId(v1.k[0])));
  }
  
  template <unsigned int p>inline typename Point_type_cpu_prop<p>::type & get(grid_key_2<p> & v1)
  {
    return boost::fusion::at_c<p>(data.mem_r.get(g1.LinId(v1.k[1],v1.k[0])));
  }
  
  template <unsigned int p>inline typename Point_type_cpu_prop<p>::type & get(grid_key_3<p> & v1)
  {
    return boost::fusion::at_c<p>(data.mem_a.get(g1.LinId(v1.k[2],v1.k[1],v1.k[0])));
  }
  
  template <unsigned int p>inline typename Point_type_cpu_prop<p>::type & get(grid_key_4<p> & v1)
  {
    return boost::fusion::at_c<p>(data.mem_a.get(g1.LinId(v1.k[3],v1.k[2],v1.k[1],v1.k[0])));
  }
  
  template <unsigned int p>inline typename Point_type_cpu_prop<p>::type & get(grid_key_d<dim,p> & v1)
  {
    return boost::fusion::at_c<p>(data.mem_a.get(g1.LinId(v1)));
  }
  
  
  template <unsigned int p>inline typename Point_type_cpu_prop<p>::type & get(grid_key_dx<dim> & v1)
  {
    return boost::fusion::at_c<p>(data.mem_a.get(g1.LinId(v1)));
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

// nested_struct

/*template<typename S, typename T>
struct nested_struct
{
	typedef typename boost::fusion::result_of::pop_back<T>::type popv;

	nested_struct< S, popv > ggn;
	S gn;
};

template<typename S>
struct nested_struct<S,boost::fusion::vector<>>
{
	nested_struct< S, boost::fusion::result_of::pop_back<T>::type > ggn;
	S<boost::fusion::result_of::end<T>> gn;
};*/



template<unsigned int p,typename T>
struct type_gpu_prop
{
	typedef typename memory_gpu_type<T>::rtype vtype;

	  typedef typename boost::fusion::result_of::at< vtype,boost::mpl::int_<p> >::type type;
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called
 *
 */

struct allocate
{
	size_t sz;

	allocate(size_t sz)
	:sz(sz){};

    template<typename T>
    void operator()(T& t) const
    {
        t.mem.allocate(sz*sizeof(t.mem_r::type));
        t.mem_r.set_pointer(t.mem.getPointer());
    }
};

template<unsigned int dim, typename T, typename Mem>
class grid_gpu
{
		//! It store all the information regarding the grid
		grid<dim,void> g1;

		//! This is the interface to allocate an resize memory
		//! and give also a representation to the allocated memory
		typename Mem::type data;

	public:

		//! Constructor it initialize the memory and give representation
		grid_gpu(std::vector<size_t> & sz)
		:g1(sz)
		{
			for_each(data,allocate(g1.size()));
		}

/*  template <unsigned int p>inline typename mem_reference<typename type_gpu_prop<p,T>::type>::type get(grid_key<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.getId()));
  }

  template <unsigned int p>inline typename mem_reference<typename type_gpu_prop<p,T>::type >::type get(grid_key_1<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.k[0]));
  }

  template <unsigned int p>inline typename mem_reference<typename type_gpu_prop<p,T>::type >::type get(grid_key_2<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.k[1],v1.k[0]));
  }

  template <unsigned int p>inline typename mem_reference<typename type_gpu_prop<p,T>::type >::type get(grid_key_3<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.k[2],v1.k[1],v1.k[0]));
  }

  template <unsigned int p>inline typename mem_reference<typename type_gpu_prop<p,T>::type >::type get(grid_key_4<p> & v1)
  {
    return boost::fusion::at_c<p>(data).get(g1.LinId(v1.k[3],v1.k[2],v1.k[1],v1.k[0]));
  }*/

		template <unsigned int p>inline typename mem_reference<typename type_gpu_prop<p,T>::type >::type get(grid_key_d<dim,p> & v1)
		{
			return boost::fusion::at_c<p>(data.mem_r).get(g1.LinId(v1));
		}

		template <unsigned int p>inline typename mem_reference<typename type_gpu_prop<p,T>::type >::type get(grid_key_dx<dim> & v1)
		{
			return boost::fusion::at_c<p>(data.mem_r).get(g1.LinId(v1));
		}

		inline size_t size()
		{
			return g1.size();
		}

		//! this function set the memory interface if required
		//! this operation is required when we define a void memory
		//! allocator
		void set_memory(memory & mem)
		{
			data.mem.set_memory(mem);
		}
};

#ifdef GPU
	use template<unsigned int n, typename T>grid_gpu<n,T> grid<n,T>
#else CPU

#endif

#endif


