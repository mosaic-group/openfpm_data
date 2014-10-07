#ifndef MAP_HPP_
#define MAP_HPP_

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
#include <boost/fusion/include/for_each.hpp>
#include "memory_conf.hpp"

#include "grid.hpp"
#include "memory_array.hpp"
#include "memory_c.hpp"

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called to copy
 * all the properties from the source grid to the destination
 * grid
 *
 * \param S
 *
 */

template<unsigned int dim, typename S>
struct copy_cpu
{
	//! size to allocate
	grid_key_dx<dim> & key;

	//! grid where we have to store the data
	S & grid_src;

	//! type of the object we have to set
	typedef typename S::type obj_type;

	//! onject we have to store
	obj_type obj;

	//! constructor it fix the size
	copy_cpu(grid_key_dx<dim> & key, S grid_src, obj_type obj)
	:grid_src(grid_src),key(key),obj(obj){};

	//! It call the copy function for each member
    template<typename T>
    void operator()(T& t) const
    {
    	t = boost::fusion::at_c<T::value>(obj);
    }
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called to copy
 * all the properties from the source grid to the destination
 * grid
 *
 * \param S
 *
 */

template<unsigned int dim, typename S>
struct copy_cpu_sd
{
	//! size to allocate
	grid_key_dx<dim> & key;

	S & grid_src;

	//! constructor it fix the size
	copy_cpu_sd(grid_key_dx<dim> & key, S grid_src, S grid_dst)
	:grid_src(grid_src),key(key){};

	//! It call the copy function for each member
    template<typename T>
    void operator()(T& t) const
    {
    	t = grid_src.template get<T::value>(key);
    }
};


/*! \brief Metafunction take T and return a reference
 *
 * Metafunction take T and return a reference
 *
 * \param T type
 *
 */

template<typename T>
struct mem_reference
{
	typedef T& type;
};

/*template<unsigned int p>
struct Point_type_cpu_prop
{
	typedef typename boost::fusion::result_of::at<Point<float>::memory_lin::vtype,boost::mpl::int_<p> >::type type;
};*/

/*! \brief This class is an helper to get the return type for get method for each property
 *
 * This class is an helper to get the return type for get method for each property
 *
 * \param p id of the property
 * \param T original boost fusion vector, T is suppose to be a boost::fusion::vector<memory_c<...>,.....>
 *
 */

template<unsigned int p,typename T>
struct type_cpu_prop
{
	//! return a boost::fusion::vector<memory_c<....>....>
	typedef typename T::memory_lin::vtype vtype;
	//! return a memory_c<...>
	typedef typename boost::fusion::result_of::at< vtype,boost::mpl::int_<p> >::type type;
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

template<unsigned int dim, typename T, typename Mem = typename memory_traits_lin< typename T::type >::type >
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
		}

		/*! \brief Return the internal grid information
		 *
		 * Return the internal grid information
		 *
		 * \return the internal grid
		 *
		 */

		grid<dim,T> getGrid()
		{
			return g1;
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
	    	//! Create and set the memory allocator
	    	data.setMemory(*new S());

	    	//! Allocate the memory and create the reppresentation
	        data.allocate(g1.size());
		}

		/*! \brief Return a plain pointer to the internal data
		 *
		 * Return a plain pointer to the internal data
		 *
		 * \return plain data pointer
		 *
		 */

		void * getPointer()
		{
			return data.mem->getPointer();
		}

		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & get(grid_key<p> & v1)
		{
			return boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1.getId())));
		}
  
		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & get(grid_key_d<dim,p> & v1)
		{
			return boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1)));
		}
  
  
		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & get(grid_key_dx<dim> & v1)
		{
			return boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1)));
		}
  
		/*! \brief Resize the space
		 *
		 * Resize the space to a new grid, the element are retained on the new grid,
		 * if the new grid is bigger the new element are now initialized, if is smaller
		 * the data are cropped
		 *
		 */

		void resize(std::vector<size_t> & sz)
		{
			//! Create a completely new grid with sz

			grid_cpu<dim,T,Mem> grid_new;

	        // We know that, if it is 1D we can safely copy the memory
	        if (dim == 1)
	        {
	        	//! 1-D copy (This case is simple we use raw memory copy because the fastest option)
	        	grid_new.data.mem.copy(data.mem);
	        }
	        else
	        {
	        	//! N-D copy

	        	//! create a source grid iterator
	        	grid_key_dx_iterator<dim> it(g1);

	        	while(it.hasNext())
	        	{
	        		// get the grid key
	        		grid_key_dx<dim> key = it.next();

	        		// create a copy element

	        		copy_cpu_sd<dim,grid_cpu<dim,T,Mem>> cp(key,*this,grid_new);

	        		// copy each property

	        		for_each(T::type,cp);
	        	}
	        }
		}

		/*! \brief set/copy an element of the grid
		 *
		 * set an element of the grid
		 *
		 * \param dx is the grid key or the position to set
		 * \param obj value to set
		 *
		 */

		inline void set(grid_key_dx<dim> dx, T & obj)
		{
			// create the object to copy the properties
    		copy_cpu<dim,grid_cpu<dim,grid_cpu<dim,T,Mem> >> cp(dx,*this);

    		// copy each property
    		for_each(T::type,cp);
		}

		/*! \brief return the size of the grid
		 *
		 * Return the size of the grid
		 *
		 */

		inline size_t size()
		{
			return g1.size();
		}

		/*! \brief Return a grid iterator
		 *
		 * Return a grid iterator, to iterate through the grid
		 *
		 */

		inline grid_key_dx_iterator<dim> getIterator()
		{
			return grid_key_dx_iterator<dim>(g1);
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
	typedef typename T::memory_int vtype;
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

template<unsigned int dim, typename T, typename Mem = typename memory_traits_inte< typename T::type >::type >
class grid_gpu
{
		//! It store all the information regarding the grid
		grid<dim,void> g1;

		//! This is the interface to allocate,resize ... memory
		//! and give also a representation to the allocated memory
//		typename T::memory_int::type data;
		Mem data;

	public:

		//! Constructor it initialize the memory and give representation
		grid_gpu(std::vector<size_t> & sz)
		:g1(sz)
		{
		}

		/*! \brief Return the internal grid information
		 *
		 * Return the internal grid information
		 *
		 * \return the internal grid
		 *
		 */

		grid<dim,void> getGrid()
		{
			return g1;
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

			//! for each element in the vector allocate the buffer
			boost::fusion::for_each(data,all);
		}

		template <unsigned int p>inline typename type_gpu_prop<p,T>::reference get(grid_key_d<dim,p> & v1)
		{
			return boost::fusion::at_c<p>(data).mem_r->operator[](g1.LinId(v1));
		}

		template <unsigned int p>inline typename type_gpu_prop<p,T>::type::reference get(grid_key_dx<dim> & v1)
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

		/*! \brief Return a grid iterator
		 *
		 * Return a grid iterator, to iterate through the grid
		 *
		 */

		inline grid_key_dx_iterator<dim> getIterator()
		{
			return grid_key_dx_iterator<dim>(g1);
		}
};

/*! device selector struct
 *
 * device selector struct, it return the correct data type for each device
 *
 */

template<unsigned int dim, typename T>
struct device_g
{
	//! cpu
	typedef grid_cpu<dim,T> cpu;
	//! gpu
	typedef grid_gpu<dim,T> gpu;
};

#endif


