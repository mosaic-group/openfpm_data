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
#include <boost/mpl/pop_front.hpp>
#include <boost/mpl/push_front.hpp>

//#include "util/boost/boost_multi_array_openfpm.hpp"
#include "util/multi_array_openfpm/multi_array_ref_openfpm.hpp"
#include "util/ct_array.hpp"
#include "memory_array.hpp"
#include "memory/memory.hpp"

#ifndef MEMORY_C_HPP_
#define MEMORY_C_HPP_

#define MEMORY_C_STANDARD 1
#define MEMORY_C_REDUCED  2

template<typename T, unsigned int impl = MEMORY_C_STANDARD, typename D = memory>
class memory_c
{

};

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
template<typename T, typename D>
class memory_c<T,MEMORY_C_STANDARD,D>
{
	//! indicate this object manage the memory.
	//! Call to functions like toKernels create views of the main data-structure. Views use memory_c create different ways to access the same sata. 
	//! We have mainly have two options. Either views manage memory, either they do not. When a view manage the memory if the main-data structure
	//! get destroyed the view is still valid, the memory retained and the last view destroy the memory.
	//! When a view does not manage the memory, if the main-data structure get destroyed the view become invalid.
	bool manage_memory = true;

	public:

	//! define T
	typedef memory_c<T> type;

	//! define a reference to T
	typedef T& reference;

	//! define T
	typedef T vtype;

	//! object that allocate memory like HeapMemory or CudaMemory
	memory * mem;

	//! object that represent the memory as an array of objects T
	memory_array<T> mem_r;

	/*! \brief This function set the object that allocate memory
	 *
	 * \param mem the memory object
	 *
	 */
	void setMemory(memory & mem)
	{
		if (manage_memory)
		{
			if (this->mem != NULL)
			{
				this->mem->decRef();

				if (this->mem->ref() == 0 && &mem != this->mem)
					delete(this->mem);
			}
			mem.incRef();
		}
		this->mem = &mem;
	}

	/*! \brief This function bind the memory_c to this memory_c as reference
	 *
	 * Bind ad reference it mean this this object does not create new memory but use the one from ref
	 * as a reference.
	 *
	 */
	bool bind_ref(const memory_c<T,MEMORY_C_STANDARD,D> & ref)
	{
	    mem = ref.mem;

		if (manage_memory)
	    {mem->incRef();}

	    //! we create the representation for the memory buffer
	    mem_r = ref.mem_r;

	    return true;
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
	bool allocate(const size_t sz, bool skip_initialization = false)
	{
		memory * mem = this->mem;

		//! We create a chunk of memory
	    mem->resize( sz*sizeof(T) );

	    //! we create the representation for this buffer
	    mem_r.initialize(mem->getPointer(),sz,mem->isInitialized() | skip_initialization);

	    return true;
	}

	//! constructor
	memory_c(bool manage_memory = true)
	:manage_memory(manage_memory),mem(NULL){}

	//! destructor
	~memory_c()
	{
		// deinitialixe mem_r
		if (manage_memory)
		{
			mem_r.deinit();
			if (mem != NULL)
			{
				mem->decRef();

				if (mem->ref() == 0)
					delete(mem);
			}
		}
	}

	/*! \brief Disable the management of memory (it is used for toKernel views)
	 * 
	 */
	void disable_manage_memory()
	{
		manage_memory = false;
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
		mem = mem_obj.mem;
		mem_obj.mem = static_cast<memory*>(mem_tmp);

		mem_obj.mem_r.swap(mem_r);
	}

	/*! \brief swap the memory
	 *
	 * swap the memory between objects
	 *
	 */
	template<typename Mem_type>
	__host__ void swap_nomode(memory_c & mem_obj)
	{
		// It would be better a dynamic_cast, unfortunately nvcc
		// does not accept it. It seems that for some reason nvcc want to
		// produce device code out of this method. While it is true
		// that this method is called inside a generic __device__ __host__
		// function, tagging this method only __host__ does not stop
		// nvcc from the intention to produce device code.
		// The workaround is to use static_cast. Another workaround (to test)
		// could be create a duplicate for_each function tagged only __host__ ,
		// but I have to intention to duplicate code

		Mem_type * mem_tmp = static_cast<Mem_type*>(mem);
		mem_tmp->swap(*static_cast<Mem_type*>(mem_obj.mem));

		mem_obj.mem_r.swap(mem_r);
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
class memory_c<multi_array<T>, MEMORY_C_STANDARD, D>
{
	//! indicate this object manage the memory.
	//! Call to functions like toKernels create views of the main data-structure. Views use memory_c create different ways to access the same sata. 
	//! We have mainly have two options. Either views manage memory, either they do not. When a view manage the memory if the main-data structure
	//! get destroyed the view is still valid, the memory retained and the last view destroy the memory.
	//! When a view does not manage the memory, if the main-data structure get destroyed the view become invalid.
	bool manage_memory = true;

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

	//! Remove the first element
	typedef typename boost::mpl::push_front<typename boost::mpl::pop_front<T>::type,boost::mpl::int_<-1>>::type Tv;

	//! define boost::mpl::int_ without boost::mpl
	template<int S> using int_ = boost::mpl::int_<S>;

	//! define the template vector size it give a number at compile time
	typedef typename boost::mpl::size<T> size_p;

	//! Define "at" meta function without boost::mpl
	template< typename S, unsigned int n> using at = boost::mpl::at<S,boost::mpl::int_<n> >;

	typedef typename at<T,0>::type base;

	//! define size_type
	typedef typename openfpm::multi_array_ref_openfpm<base,size_p::value,Tv>::size_type size_type;

	public:

	/*! \brief This function set the object that allocate memory
	 *
	 * \param mem the memory object
	 *
	 */

	void setMemory(memory & mem)
	{
		if (manage_memory)
		{
			if (this->mem != NULL)
			{
				this->mem->decRef();

				if (this->mem->ref() == 0 && &mem != this->mem)
					delete(this->mem);
			}
			mem.incRef();
		}
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


	/*! \brief This function bind the memory_c to this memory_c as reference
	 *
	 * Bind ad reference it mean this this object does not create new memory but use the one from ref
	 * as a reference.
	 *
	 */
	bool bind_ref(const memory_c<multi_array<T>, MEMORY_C_STANDARD, D> & ref)
	{
	    mem = ref.mem;
		if (manage_memory)
	    {mem->incRef();}

	    //! we create the representation for the memory buffer
	    mem_r = ref.mem_r;

	    return true;
	}

	/*! \brief This function allocate memory and associate the representation to mem_r
	 *
	 * This function allocate memory and associate the representation of that chunk of
	 * memory to mem_r
	 *
	 */
	bool allocate(const size_t sz, bool skip_initialization = false)
	{
		memory * mem = this->mem;

		//! We create a chunk of memory
		mem->resize( sz*mult<T,size_p::value-1>::value*sizeof(base) );

		openfpm::multi_array_ref_openfpm<base,size_p::value,Tv,storage_order> tmp(static_cast<base *>(mem->getPointer()), sz);

	    //! we create the representation for the memory buffer
	    mem_r.swap(tmp);

	    return true;
	}

	//! define the type of the multi_array vector minus one of the original size
	//! basically we remove the index 0 of the multi_array
	typedef boost::multi_array<base,size_p::value> type;

	//! define the storage order used to reference multidimensional arrays
	typedef typename openfpm::ofp_storage_order<size_p::value>::value storage_order;

	//! Reference to an object to allocate memory
	D * mem;

	//! object that represent the memory as a multi-dimensional array of objects T
	openfpm::multi_array_ref_openfpm<base,boost::mpl::size<T>::value,Tv,storage_order> mem_r;

	//! constructor
	memory_c(bool manage_memory = true)
	:manage_memory(manage_memory),mem(NULL),
	 mem_r(static_cast<base *>(NULL),0)
	{}

	//! destructor
	~memory_c()
	{
		// Indicate that this memory_c does not manage the memory
		if (manage_memory)
		{
			if (mem != NULL)
			{
				mem->decRef();
				if (mem->ref() == 0)
					delete(mem);
			}
		}
	}

	/*! \brief Disable the management of memory (it is used for toKernel views)
	 * 
	 */
	void disable_manage_memory()
	{
		manage_memory = false;
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
		mem = mem_obj.mem;
		mem_obj.mem = static_cast<D*>(mem_tmp);

		mem_r.swap(mem_obj.mem_r);
	}


	/*! \brief swap the memory
	 *
	 * While the previous one swap the mem the swap_nomode call the swap member of the mem object
	 *
	 * \note calling the swap member require knowledge of the object type, we cannot work on abtsract objects
	 *
	 * swap the memory between objects
	 *
	 */
	template<typename Mem_type>
	__host__ void swap_nomode(memory_c & mem_obj)
	{
		// It would be better a dynamic_cast, unfortunately nvcc
		// does not accept it. It seems that for some reason nvcc want to
		// produce device code out of this method. While it is true
		// that this function is called inside a generic __device__ __host__
		// function, tagging this method only __host__ does not stop
		// nvcc from the intention to produce device code.
		// The workaround is to use static_cast. Another workaround (to test)
		// could be create a duplicate for_each function tagged only __host__ ,
		// but I have to intention to duplicate code

		Mem_type * mem_tmp = static_cast<Mem_type*>(mem);
		mem_tmp->swap(*static_cast<Mem_type*>(mem_obj.mem));

		mem_obj.mem_r.swap(mem_r);
	}
};

//! Partial specialization for scalar N=0
template<typename T>
struct array_to_vmpl
{
};


//! Partial specialization for N=1
template<typename T,size_t N1>
struct array_to_vmpl<T[N1]>
{
	//! the internal array primitive information represented into a boost mpl vector
	typedef boost::mpl::vector<T,boost::mpl::int_<N1>> prim_vmpl;
};


//! Partial specialization for N=2
template<typename T,size_t N1,size_t N2>
struct array_to_vmpl<T[N1][N2]>
{
	//! the internal array primitive information represented into a boost mpl vector
	typedef boost::mpl::vector<T,boost::mpl::int_<N1>,
			                     boost::mpl::int_<N2>> prim_vmpl;
};

/*! \brief OpenFPM use memory_c<multi_array<T> ..... > to implement the structure of array layout
 *
 * This mean that the object returned by mem_r are complex objects that represent the memory view, these view has
 * the purpose to hook at compile time the operator[] to give the feeling of using an array
 *
 * This view depend from the template parameter Tv in the member mem_r, that as you can see is difficult to reconstruct
 * In some case deduction does not work because too complex. So we have to compute this type. this function does this
 * given a type like float[3], it produce the Tv parameter
 *
 *
 */
template<typename T>
struct to_memory_multi_array_ref_view
{
	// first we convert the type into a boost vector containing the primitive, followed by array dimension
	typedef typename array_to_vmpl<T>::prim_vmpl prim_vmpl;

	// Than we operate at compile-time the same operation memory_c<multi_array does
	//! Remove the first element (this is the Tv parameter of )
	typedef typename boost::mpl::push_front<typename boost::mpl::pop_front<prim_vmpl>::type,boost::mpl::int_<-1>>::type vmpl;


};

#endif

