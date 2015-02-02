#ifndef MAP_HPP_
#define MAP_HPP_

#include "config.h"

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
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>
#include "memory_conf.hpp"
#include "meta_copy.hpp"
#include "Memleak_check.hpp"


#include "grid.hpp"
#include "Encap.hpp"
#include "memory_array.hpp"
#include "memory_c.hpp"
#include <vector>

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one object into one target
 * grid  element in a generic way for a
 * generic object T with variable number of property
 *
 * \param dim Dimensionality
 * \param S type of grid
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

	//! type of the object boost::sequence
	typedef typename S::type::type ov_seq;

	//! object we have to store
	obj_type & obj;

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param key which element we are modifying
	 * \param grid_src grid we are updating
	 * \param obj object we have to set in grid_src
	 *
	 */
	copy_cpu(grid_key_dx<dim> & key, S & grid_src, obj_type & obj)
	:key(key),grid_src(grid_src),obj(obj){};

	//! It call the copy function for each property
    template<typename T>
    void operator()(T& t) const
    {
    	// This is the type of the object we have to copy
    	typedef typename boost::fusion::result_of::at_c<ov_seq,T::value>::type copy_type;

    	// Remove the reference from the type to copy
    	typedef typename boost::remove_reference<copy_type>::type copy_rtype;

    	meta_copy<copy_rtype> cp(grid_src.template get<T::value>(key),boost::fusion::at_c<T::value>(obj.data));
    }
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one source grid element into one target
 * grid element in a generic way for an object T with variable
 * number of property
 *
 * \param dim Dimensionality
 * \param S grid type
 *
 */

template<unsigned int dim, typename S>
struct copy_cpu_sd
{
	//! size to allocate
	grid_key_dx<dim> & key;

	//! Source grid
	S & grid_src;

	//! Destination grid
	S & grid_dst;

	//! type of the object boost::sequence
	typedef typename S::type::type ov_seq;

	//! constructor it fix the size
	copy_cpu_sd(grid_key_dx<dim> & key, S grid_src, S grid_dst)
	:key(key),grid_src(grid_src),grid_dst(grid_dst){};

	//! It call the copy function for each member
    template<typename T>
    void operator()(T& t) const
    {
    	// This is the type of the object we have to copy
    	typedef typename boost::fusion::result_of::at_c<ov_seq,T::value>::type copy_type;

    	// Remove the reference from the type to copy
    	typedef typename boost::remove_reference<copy_type>::type copy_rtype;

    	meta_copy<copy_rtype> cp(grid_dst.template get<T::value>(key),grid_src.template get<T::value>(key));
    }
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one source grid element into one target
 * grid element in a generic way for an object T with variable
 * number of property
 *
 * \param dim Dimensionality
 * \param S grid type
 *
 */

template<unsigned int dim, typename S>
struct copy_cpu_sd_k
{
	//! source key
	grid_key_dx<dim> & key_s;

	//! destination key
	grid_key_dx<dim> & key_d;

	//! Source grid
	S & grid_src;

	//! Destination grid
	S & grid_dst;

	//! type of the object boost::sequence
	typedef typename S::type::type ov_seq;

	//! constructor it fix the size
	copy_cpu_sd_k(grid_key_dx<dim> & key_s, grid_key_dx<dim> & key_d, S & grid_src, S & grid_dst)
	:key_s(key_s),key_d(key_d),grid_src(grid_src),grid_dst(grid_dst){};

	//! It call the copy function for each member
    template<typename T>
    void operator()(T& t) const
    {
    	// This is the type of the object we have to copy
    	typedef typename boost::fusion::result_of::at_c<ov_seq,T::value>::type copy_type;

    	// Remove the reference from the type to copy
    	typedef typename boost::remove_reference<copy_type>::type copy_rtype;

    	meta_copy<copy_rtype> cp(grid_dst.template get<T::value>(key_s),grid_src.template get<T::value>(key_d));
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
		//! Access the key
		typedef grid_key_dx<dim> access_key;

		//! boost::vector that describe the data type
		typedef typename T::type T_type;

		//! This is an header that store all information related to the grid
		grid<dim,T> g1;

		//! This is the interface to allocate an resize memory
		//! and give also a representation to the allocated memory
		Mem data;

		/*! \brief Get 1D vector with the
		 *
		 * Get std::vector with element 0 to dim set to 0
		 *
		 */

		std::vector<size_t> getV()
		{
			std::vector<size_t> tmp;

			for (unsigned int i = 0 ; i < dim ; i++)
			{
				tmp.push_back(0);
			}

			return tmp;
		}

	public:

		//! Memory traits
		typedef Mem memory_t;

		//! Object container for T, it is the return type of get_o it return a object type trough
		// you can access all the properties of T
		typedef encapc<dim,T,Mem> container;

		// The object type the grid is storing
		typedef T type;

		//! Default constructor
		grid_cpu()
		:g1(getV())
	    {
	    }

		//! Set the grid dimensions
		void setDimensions(std::vector<size_t> & sz)
		{
			g1.setDimension(sz);
		}

		/*! \brief create a grid from another grid
		 *
		 * \param the grid to copy
		 * \param S memory type, used for template deduction
		 *
		 */
		template<typename S> grid_cpu(const grid_cpu & g, S & mem)
		{
			swap(g.duplicate<S>());
		}

		//! Constructor allocate memory and give them a representation
		grid_cpu(std::vector<size_t> & sz)
		:g1(sz)
		{
		}

		//! Constructor allocate memory and give them a representation
		grid_cpu(std::vector<size_t> && sz)
		:g1(sz)
		{
		}

		//! Constructor allocate memory and give them a representation
		grid_cpu(size_t (& sz)[dim])
		:g1(sz)
		{
		}

		/*! \brief create a duplicated version of the grid
		 *
		 */

		template<typename S> grid_cpu<dim,T,Mem> duplicate()
		{
			//! Create a completely new grid with sz

			grid_cpu<dim,T,Mem> grid_new(g1.getSize());

			//! Set the allocator and allocate the memory
			grid_new.template setMemory<S>();

	        // We know that, if it is 1D we can safely copy the memory
	        if (dim == 1)
	        {
	        	//! 1-D copy (This case is simple we use raw memory copy because is the fastest option)
	        	grid_new.data.mem->copy(*data.mem);
	        }
	        else
	        {
	        	//! N-D copy

	        	//! create a source grid iterator
	        	grid_key_dx_iterator<dim> it(g1);

	        	while(it.isNext())
	        	{
	        		// get the grid key
	        		grid_key_dx<dim> key = it.get();

	        		// create a copy element

	        		copy_cpu_sd<dim,grid_cpu<dim,T,Mem>> cp(key,*this,grid_new);

	        		// copy each property for each point of the grid

	        		boost::mpl::for_each< boost::mpl::range_c<int,0,T::max_prop> >(cp);

	        		++it;
	        	}
	        }

	        // copy grid_new to the base

	        return grid_new;
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

		/*! \brief Get the reference of the selected element
		 *
		 * Get the reference of the selected element
		 *
		 * \param p property to get (is an integer)
		 * \param v1 grid_key that identify the element in the grid
		 *
		 */
		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & get(grid_key<p> & v1)
		{
#ifdef MEMLEAK_CHECK
			check_valid(&boost::fusion::at_c<p>(&data.mem_r->operator[](g1.LinId(v1.getId()))));
#endif
			return boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1.getId())));
		}
  
		/*! \brief Get the reference of the selected element
		 *
		 * Get the reference of the selected element
		 *
		 * \param v1 grid_key that identify the element in the grid
		 *
		 */
		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & get(grid_key_d<dim,p> & v1)
		{
#ifdef MEMLEAK_CHECK
			check_valid(&boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1))));
#endif
			return boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1)));
		}
  
		/*! \brief Get the reference of the selected element
		 *
		 * Get the reference of the selected element
		 *
		 * \param v1 grid_key that identify the element in the grid
		 *
		 */
		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & get(grid_key_dx<dim> & v1)
		{
#ifdef MEMLEAK_CHECK
			check_valid(&boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1))),sizeof(typename type_cpu_prop<p,T>::type));
#endif
			return boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1)));
		}


		/*! \brief Get the of the selected element as a boost::fusion::vector
		 *
		 * Get the selected element as a boost::fusion::vector
		 *
		 * \param v1 grid_key that identify the element in the grid
		 *
		 */
		inline encapc<dim,T,Mem> get_o(const grid_key_dx<dim> & v1)
		{
#ifdef MEMLEAK_CHECK
			check_valid(&data.mem_r->operator[](g1.LinId(v1)),sizeof(typename type_cpu_prop<p,T>::type));
#endif
			return encapc<dim,T,Mem>(data.mem_r->operator[](g1.LinId(v1)));
		}


		/*! \brief Resize the space
		 *
		 * Resize the space to a new grid, the element are retained on the new grid,
		 * if the new grid is bigger the new element are now initialized, if is smaller
		 * the data are cropped
		 *
		 * \param sz reference to an array of dimension dim
		 *
		 */

		template<typename S> void resize(size_t (& sz)[dim])
		{
			//! Create a completely new grid with sz

			grid_cpu<dim,T,Mem> grid_new(sz);

			//! Set the allocator and allocate the memory
			grid_new.template setMemory<S>();

	        // We know that, if it is 1D we can safely copy the memory
	        if (dim == 1)
	        {
	        	//! 1-D copy (This case is simple we use raw memory copy because is the fastest option)
	        	grid_new.data.mem->copy(*data.mem);
	        }
	        else
	        {
	        	//! N-D copy

	        	//! create a source grid iterator
	        	grid_key_dx_iterator<dim> it(g1);

	        	while(it.isNext())
	        	{
	        		// get the grid key
	        		grid_key_dx<dim> key = it.get();

	        		// create a copy element

	        		copy_cpu_sd<dim,grid_cpu<dim,T,Mem>> cp(key,*this,grid_new);

	        		// copy each property for each point of the grid

	        		boost::mpl::for_each< boost::mpl::range_c<int,0,T::max_prop> >(cp);

	        		++it;
	        	}
	        }

	        // copy grid_new to the base

	        this->swap(grid_new);
		}

		/*! \brief Resize the space
		 *
		 * Resize the space to a new grid, the element are retained on the new grid,
		 * if the new grid is bigger the new element are now initialized, if is smaller
		 * the data are cropped
		 *
		 */

		template<typename S> void resize(std::vector<size_t> & sz)
		{
			// array containing the size of the grid
			size_t sz_a[dim];

			// fill the array
			for (int i = 0 ; i < dim ; i++)
			{sz_a[i] = sz[i];}

			// resize
			resize<S>(sz_a);
		}

		/*! \brief It move the allocated object from one grid to another
		 *
		 * It move the allocated object from one grid to another, after this
		 * call the argument grid is no longer valid
		 *
		 * \param grid to move/copy
		 *
		 */

		void swap(grid_cpu<dim,T,Mem> & grid)
		{
			// move the data
			data.move_copy(grid.data);

			// move the grid info
			g1 = grid.g1;
		}

		/*! \brief It move the allocated object from one grid to another
		 *
		 * It move the allocated object from one grid to another, after this
		 * call the argument grid is no longer valid
		 *
		 * \param grid to move/copy
		 *
		 */

		void swap(grid_cpu<dim,T,Mem> && grid)
		{
			swap(grid);
		}

		/*! \brief set an element of the grid
		 *
		 * set an element of the grid
		 *
		 * \param dx is the grid key or the position to set
		 * \param obj value to set
		 *
		 */

		inline void set(grid_key_dx<dim> dx, T & obj)
		{
#ifdef DEBUG
			// Check that the element exist

			for (int i = 0 ; i < dim ; i++)
			{
				if (dx.get(i) >= g1.size(i))
				{
					std::cerr << "Error map_grid.hpp: out of bound" << "\n";
				}
			}
#endif

			// create the object to copy the properties
    		copy_cpu<dim,grid_cpu<dim,T,Mem>> cp(dx,*this,obj);

    		// copy each property
    		boost::mpl::for_each< boost::mpl::range_c<int,0,T::max_prop> >(cp);
		}

		/*! \brief set an element of the grid
		 *
		 * set an element of the grid from another element of another grid
		 *
		 * \param key1 element of the grid to set
		 * \param g source grid
		 * \param element of the source grid to copy
		 *
		 */

		inline void set(grid_key_dx<dim> key1, grid_cpu<dim,T,Mem> & g, grid_key_dx<dim> key2)
		{
			//create the object to copy the properties
    		copy_cpu_sd_k<dim,grid_cpu<dim,T,Mem>> cp(key1,key2,*this,g);

    		// copy each property for each point of the grid

    		boost::mpl::for_each< boost::mpl::range_c<int,0,T::max_prop> >(cp);

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

		/*! \brief Return a sub-grid iterator
		 *
		 * Return a sub-grid iterator, to iterate through the grid
		 *
		 */

		inline grid_key_dx_iterator_sub<dim> getSubIterator(grid_key_dx<dim> & start, grid_key_dx<dim> & stop)
		{
			return grid_key_dx_iterator_sub<dim>(g1,start,stop);
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

		/*! \brief Return a grid iterator over all the point with the exception
		 *   of the ghost part
		 *
		 * Return a grid iterator over all the point with the exception of the
		 * ghost part
		 *
		 */

		inline grid_key_dx_iterator_sub<dim> getIteratorDomain()
		{
			// get the starting point and the end point of the real domain

			return grid_key_dx_iterator_sub<dim>(g1.getDomainStart(),g1.getDomainStop());
		}
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
	//! Access the key
	typedef grid_key_dx<dim> access_key;

	//! It store all the information regarding the grid
	grid<dim,void> g1;

	//! This is the interface to allocate,resize ... memory
	//! and give also a representation to the allocated memory
	Mem data;

public:

	//! Memory traits
	typedef Mem memory_t;

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef encapg<dim,T,Mem> container;

	// The object type the grid is storing
	typedef T type;

	//! Default constructor
	grid_gpu()
	{
	}

	//! Set the grid dimensions
	void setDimensions(std::vector<size_t> & sz)
	{
		g1.setDimension(sz);
	}

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

	template <unsigned int p>inline typename type_gpu_prop<p,T>::type::reference get(grid_key_d<dim,p> & v1)
	{
		return boost::fusion::at_c<p>(data).mem_r->operator[](g1.LinId(v1));
	}

	template <unsigned int p>inline typename type_gpu_prop<p,T>::type::reference get(grid_key_dx<dim> & v1)
	{
		return boost::fusion::at_c<p>(data).mem_r->operator[](g1.LinId(v1));
	}

	/*! \brief Get the of the selected element as a boost::fusion::vector
	 *
	 * Get the selected element as a boost::fusion::vector
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 */
	inline encapg<dim,T,Mem> get_o(grid_key_dx<dim> & v1)
	{
		return encapg<dim,T,Mem>(data,v1,g1);
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


	/*! \brief Return a sub-grid iterator
	 *
	 * Return a sub-grid iterator, to iterate through the grid
	 *
	 */

	inline grid_key_dx_iterator_sub<dim> getSubIterator(grid_key_dx<dim> & start, grid_key_dx<dim> & stop)
	{
		return grid_key_dx_iterator_sub<dim>(g1,start,stop);
	}

	/*! \brief Swap the memory of another grid
	 *
	 * Swap the memory of another grid
	 *
	 * \obj Memory to swap with
	 *
	 */
	void swap(grid_gpu<dim,T,Mem> & obj)
	{
		g1.swap(obj.g1);
		data.swap(obj.data);
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


