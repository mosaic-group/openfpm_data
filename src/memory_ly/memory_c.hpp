/*
 * memory_c.hpp
 *
 *  Created on: Aug 17, 2014
 *      Author: Pietro Incardona
 */

#include <boost/multi_array.hpp>
#include <boost/fusion/mpl.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/mpl/vector.hpp>
#include <array>
#include "util/ct_array.hpp"
#include "memory_array.hpp"
#include "memory/memory.hpp"

#ifndef MEMORY_C_HPP_
#define MEMORY_C_HPP_

/*!
 * \brief This class is a container for the memory interface like HeapMemory CudaMemory
 *
 * It store the object used to allocate memory and a representation of this memory as an array of objects T
 *
 * It is mainly used by memory_conf to create the correct memory layout
 *
 * \see memory_traits_inte memory_traits_lin
 *
 */

template<typename T, typename D = memory>
class memory_c
{
	public:

	//! define T
	typedef memory_c<T> type;

	//! define a reference to T
	typedef T& reference;

	//! define T
	typedef T vtype;

	//! object that allocate memory like HeapMemory or CudaMemory
	D * mem;

	//! object that represent the memory as an array of objects T
	memory_array<T> * mem_r;

	/*! \brief This function set the object that allocate memory
	 *
	 * \param mem the memory object
	 *
	 */
	void setMemory(memory & mem)
	{
		if (this->mem != NULL)
		{
			this->mem->decRef();

			if (this->mem->ref() == 0 && &mem != this->mem)
				delete(this->mem);
		}
//			this->mem->decRef();
		mem.incRef();
		this->mem = &mem;
	}

	/*! \brief This function get the object that allocate memory
	 *
	 * \return memory object to allocate memory
	 *
	 */

	memory& getMemory()
	{
		return *this->mem;
	}

	/*! \brief This function get the object that allocate memory
	 *
	 * \return memory object to allocate memory
	 *
	 */

	const memory& getMemory() const
	{
		return *this->mem;
	}

	/*! \brief This function allocate memory
	 *
	 */
	bool allocate(const size_t sz)
	{
		memory * mem = this->mem;

		//! We create a chunk of memory
	    mem->resize( sz*sizeof(T) );

	    //! we create the representation for this buffer
	    mem_r = new memory_array<T>(mem->getPointer(),sz,mem->isInitialized());

	    return true;
	}

	/*! \brief It absorb the allocated object from another memory_c
	 *
	 * \param mem_c Memory object
	 *
	 */
	void move_copy(memory_c & mem_c)
	{
		//if mem is already allocated, deallocate it
		if (mem != NULL)
			delete(mem);

		// if mem_r is allocated delete it
		if (mem_r != NULL)
			delete(mem_r);

		// move the pointer
		mem = mem_c.mem;

		// Null the pointer to avoid double deallocation
		mem_c.mem = NULL;

		// move the pointer
		mem_r = mem_c.mem_r;

		// Null the pointer to avoid double deallocation
		mem_c.mem_r = NULL;
	}

	//! constructor
	memory_c():mem(NULL),mem_r(NULL){}

	//! destructor
	~memory_c()
	{
		if (mem != NULL)
		{
			mem->decRef();

			if (mem->ref() == 0)
				delete(mem);
		}
		delete(mem_r);
	}

	/*! \brief swap the memory
	 *
	 * swap the memory between objects
	 *
	 */
	void swap(memory_c & mem_obj)
	{
		// Save on temporal

		void * mem_tmp = static_cast<void*>(mem);
		void * mem_r_tmp = static_cast<void*>(mem_r);

		// swap the memory between objects

		mem = mem_obj.mem;
		mem_r = mem_obj.mem_r;

		mem_obj.mem = static_cast<D*>(mem_tmp);
		mem_obj.mem_r = static_cast<memory_array<T>*>(mem_r_tmp);
	}

};

/*! \brief This class is a trick to indicate the compiler a specific
 *  specialization pattern
 *
 * In particular it say that a multidimensional array has been found and
 * need a special treatment, T is suppose to be a boost::mpl::vector of
 * unsigned int indicating each dimension. T has to be a type of known size at
 * compile time
 *
 */
template<typename T>
class multi_array
{
	typedef T type;
};

/*! \brief This class is a trick to indicate the compiler a specific
 *  specialization pattern
 *
 * This class is a trick to indicate the compiler a specific specialization
 * pattern, in particular it say that T is a key and
 * need a special treatment. T is suppose to be a boost::mpl::vector of
 * any type, (Actually does not have an application but is a generalization
 * of multi-array). T has to be a type of known size at compile time
 *
 */

template<typename T>
class key
{
	typedef T type;
};


/*! \brief this class multiply all the elements in a boost::mpl::vector excluding the first element
 *
 * \param T expecting a boost::mpl::vector
 *
 */
template<typename T, unsigned int N>
struct mult
{
	enum { value = mult<T,N-1>::value * boost::mpl::at<T,boost::mpl::int_<N>>::type::value };
};

template <typename T>
struct mult<T,1>
{
	enum { value = boost::mpl::at<T,boost::mpl::int_<1>>::type::value };
};

/*! \brief Specialization of memory_c for multi_array
 *
 * Specialization of memory_c for multi_array
 *
 * It is mainly used by memory_conf to create the correct layout
 *
 * \see memory_traits_inte memory_traits_lin
 *
 * \tparam T is suppose to be a boost::mpl::vector specifing at position 0 the type and at
 * position 1 to N the dimensions size of the multi_array
 *
 * \tparam D object that allocate memory
 *
 */

template<typename T, typename D>
class memory_c<multi_array<T>, D>
{
	//! define T
//	typedef T type;

	//! define boost::mpl::int_ without boost::mpl
	template<int S> using int_ = boost::mpl::int_<S>;

	//! define the template vector size it give a number at compile time
	typedef typename boost::mpl::size<T> size_p;

	//! Define "at" meta function without boost::mpl
	template< typename S, unsigned int n> using at = boost::mpl::at<S,boost::mpl::int_<n> >;

	typedef typename at<T,0>::type base;

	//! define size_type
	typedef typename boost::multi_array_ref<base,size_p::value>::size_type size_type;


	/*! \brief In combination with generate_array is used to produce array at compile-time
	 *
	 * In combination with generate_array is used to produce at compile-time
	 * arrays like {0,N-1,.........0} used in boost::multi_array for index ordering
	 *
	 */

	template<size_t index,size_t N> struct ordering
	{
	    enum { value = (N-index) % N };
	};

	/*! \brief In combination with generate_array is used to produce array at compile-time
	 *
	 * In combination with generate_array is used to produce at compile-time
	 * arrays like {true,true,.........true} used in boost::multi_array to
	 * define ascending order
	 *
	 */

	template<size_t index,size_t N> struct ascending
	{
	    enum { value = true };
	};

	public:

	/*! \brief This function set the object that allocate memory
	 *
	 * \param mem the memory object
	 *
	 */

	void setMemory(memory & mem)
	{
		if (this->mem != NULL)
		{
			this->mem->decRef();

			if (this->mem->ref() == 0 && &mem != this->mem)
				delete(this->mem);
		}
		mem.incRef();
		this->mem = &mem;
	}

	/*! \brief This function get the object that allocate memory
	 *
	 * \return memory object to allocate memory
	 *
	 */

	memory& getMemory()
	{
		return *this->mem;
	}

	/*! \brief This function allocate memory and associate the representation to mem_r
	 *
	 * This function allocate memory and associate the representation of that chunk of
	 * memory to mem_r
	 *
	 */
	bool allocate(const size_t sz)
	{
		memory * mem = this->mem;

		//! We create a chunk of memory
	    mem->resize( sz*mult<T,size_p::value-1>::value*sizeof(base) );

	    // We create an array dims from the boost::mpl::vector
	    typedef typename generate_array_vector<size_type,T>::result dims;

	    //! buffer to store the dimensions of the full multi_array buffer
	    std::array<size_type ,size_p::value> dimensions;

	    // fill runtime, and the other dimensions
	    dimensions[0] = sz;
	    for (int i = 0 ; i < size_p::value-1 ; i++)
	    {dimensions[i+1] = dims::data[i];}

	    //! we generate the ordering buffer ord::data = {0,N-1 ...... 1 }
	    typedef typename generate_array<typename boost::multi_array<T,size_p::value>::size_type,size_p::value, ordering>::result ord;

	    // we generate the ascending buffer
	    typedef typename generate_array<bool,size_p::value, ascending>::result asc;

	    //! we create the representation for the memory buffer
	    mem_r = new boost::multi_array_ref<base,size_p::value>(static_cast<base *>(mem->getPointer()),dimensions,boost::general_storage_order<size_p::value>(ord::data,asc::data));

	    return true;
	}

	//! define the type of the multi_array vector minus one of the original size
	//! basically we remove the index 0 of the multi_array
	typedef boost::multi_array<base,size_p::value> type;

	//! Reference to an object to allocate memory
	D * mem;

	//! object that represent the memory as an multi-dimensional array of objects T
	boost::multi_array_ref<base,boost::mpl::size<T>::value> * mem_r;

	//! constructor
	memory_c():mem(NULL),mem_r(NULL){}

	//! destructor
	~memory_c()
	{
		if (mem != NULL)
		{
			mem->decRef();
			if (mem->ref() == 0)
				delete(mem);
		}
		delete(mem_r);
	}

	//! set the device memory interface, the object that allocate memory
	void set_mem(memory & mem)
	{
		this->mem = &mem;
	}

	/*! \brief swap the memory
	 *
	 * swap the memory between objects
	 *
	 */
	void swap(memory_c & mem_obj)
	{
		// Save on temporal

		void * mem_tmp = static_cast<void*>(mem);
		void * mem_r_tmp = static_cast<void*>(mem_r);

		// swap the memory between objects

		mem = mem_obj.mem;
		mem_r = mem_obj.mem_r;

		mem_obj.mem = static_cast<D*>(mem_tmp);
		mem_obj.mem_r = static_cast<decltype(mem_obj.mem_r)>(mem_r_tmp);
	}
};

#endif

