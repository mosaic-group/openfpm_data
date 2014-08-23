

//! Warning: apparently you cannot used nested boost::mpl with boost::fusion
//! can create template circularity, this include avoid the problem
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include <boost/type_traits.hpp>

#include "grid.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "memory_array.hpp"
#include "memory_c.hpp"

#ifndef MAP_HPP_
#define MAP_HPP_


// Compiled based specialization

template<typename T>
struct mem_reference
{
	typedef T& type;
};

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

/*! \brief This class is an helper to get the return type of get for each property
 *
 * This class is an helper to get the return type of get for each property
 *
 * \param p id of the property
 * \param T original boost fusion vector, T is suppose to be a boost::fusion::vector<memory_c<...>,.....>
 *
 */

template<unsigned int p,typename T>
struct type_gpu_prop
{
	//! return a boost::fusion::vector<memory_c<....>....>
	typedef typename memory_gpu_type<T>::type vtype;
	//! return a memory_c<...>
	typedef typename boost::fusion::result_of::at< vtype,boost::mpl::int_<p> >::type rtype;
	//! remove the reference
	typedef typename boost::remove_reference<rtype>::type mtype;
	//! get the base type that the buffer is storing
	typedef typename mtype::type type;
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called
 *
 * \param T Type of memory allocator
 *
 */

template<typename S>
struct allocate
{
	//! size to allocate
	size_t sz;

	//! constructor it fix the size
	allocate(size_t sz)
	:sz(sz){};

	//! It call the allocate function for each member
    template<typename T>
    void operator()(T& t) const
    {
    	//! Create and set the memory allocator
    	t.setMemory(*new S());

    	//! Allocate the memory and create the reppresentation
        t.allocate(sz);
    }
};

template<unsigned int dim, typename T, typename Mem>
class grid_gpu
{
		//! It store all the information regarding the grid
		grid<dim,void> g1;

		//! This is the interface to allocate,resize ... memory
		//! and give also a representation to the allocated memory
		Mem data;

	public:

		//! Constructor it initialize the memory and give representation
		grid_gpu(std::vector<size_t> & sz)
		:g1(sz)
		{
		}

		/*! \brief Create the object that provide memory
		 *
		 * Create the object that provide memory
		 *
		 * \param T memory
		 *
		 */

		template<typename S> void setMemory()
		{
			//! Create an allocate object
			allocate<S> all(g1.size());

			for_each(data,all);
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

		template <unsigned int p>inline typename type_gpu_prop<p,T>::reference get(grid_key_d<dim,p> & v1)
		{
			return boost::fusion::at_c<p>(data).mem_r->operator[](g1.LinId(v1));
		}

		template <unsigned int p>inline typename type_gpu_prop<p,T>::type::reference get(grid_key_dx<dim> & v1)
		{
			return boost::fusion::at_c<p>(data).mem_r->operator[](g1.LinId(v1));
		}

		template <unsigned int p>inline typename type_gpu_prop<p,T>::type diocane(grid_key_dx<dim> & v1)
		{
			return boost::fusion::at_c<p>(data).mem_r->operator[](g1.LinId(v1));
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

#endif


