/*
 * map_vector.hpp
 *
 *  Created on: Aug 30, 2014
 *      Author: Pietro Incardona
 */

#ifndef MAP_VECTOR_HPP
#define MAP_VECTOR_HPP

#include "Grid/map_grid.hpp"
#include "memory/HeapMemory.hpp"
#include "vect_isel.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#define PAGE_ALLOC 1

namespace openfpm
{
	/*! \brief Vector iterator
	 *
	 */

	class vector_key_iterator
	{
		//! Linearized end element
		size_t end;

	protected:

		//! Actual key
		size_t gk;

	public:

		/*! \brief Constructor require the size of the vector
		 *
		 * \param end size of the vector
		 */
		vector_key_iterator(size_t end)
		: end(end),gk(0)
		{}

		/*! \brief Get the next element
		 *
		 * Get the next element
		 *
		 * \return the next grid_key
		 *
		 */
		vector_key_iterator operator++()
		{
			//! increment the first index

			gk++;

			return *this;
		}

		/*! \brief Set the dimension
		 *
		 * \param d is the dimension (IGNORED is by default 0)
		 * \param sz set the counter to sz
		 *
		 */
		void set(int d, size_t sz)
		{
			// set the counter dim to sz

			gk = sz;
		}

		/*! \brief Check if there is the next element
		 *
		 * Check if there is the next element
		 *
		 * \return true if there is the next, false otherwise
		 *
		 */
		bool isNext()
		{
			if (gk < end)
			{
				//! we did not reach the end of the grid

				return true;
			}

			//! we reach the end of the grid
			return false;
		}

		/*! \brief Get the actual key
		 *
		 * Get the actual key
		 *
		 * \return the actual key
		 *
		 */
		size_t get()
		{
			return gk;
		}
	};

	/*! \brief Grow policy define how the vector should grow every time we exceed the size
	 *
	 * In this case it return the requested size
	 *
	 */

	class grow_policy_identity
	{
	public:

		/*! \brief It say how much the vector must grow
		 *
		 * \param original size
		 * \param requested size
		 *
		 * \return how much to grow
		 *
		 */
		static size_t grow(size_t original, size_t requested)
		{
			return requested;
		}
	};

	/*! \brief Grow policy define how the vector should grow every time we exceed the size
	 *
	 * In this case it double up the size
	 *
	 */

	class grow_policy_double
	{
	public:

		/*! \brief It say how much the vector must grow
		 *
		 * \param original size
		 * \param requested size
		 *
		 * \return how much to grow
		 *
		 */
		static size_t grow(size_t original, size_t requested)
		{
			size_t grow = (original == 0)?1:original;
			while (grow < requested)	{grow *= 2;}
			return grow;
		}
	};

	//! default grow policy
	typedef grow_policy_double vector_grow_policy_default;

	/*! \brief Grow policy define how the vector should grow every time we exceed the size
	 *
	 * In this case it increase of 4096 elements
	 *
	 */

	class grow_policy_page
	{
	public:

		/*! \brief It say how much the vector must grow
		 *
		 * \param original size
		 * \param requested size
		 *
		 * \return how much to grow
		 *
		 */
		static size_t grow(size_t original, size_t requested)
		{
			return (requested / PAGE_ALLOC) * PAGE_ALLOC + PAGE_ALLOC;
		}
	};

	/*! \brief Select the correct vector implementation based on the type checking
	 *
	 * \see device_cpu
	 *
	 * case OPENFPM native
	 *
	 */
	template<typename T, unsigned int impl>
	struct device_cpu_impl
	{
		//! grid
		typedef grid_cpu<1,T> type;
	};

	/*! \brief Select the correct vector implementation based on the type checking
	 *
	 * \see device_cpu
	 *
	 * case STD_VECTOR native
	 *
	 */
	template<typename T>
	struct device_cpu_impl<T,STD_VECTOR>
	{
		//! void
		typedef void type;
	};

	/*! \brief device selector struct
	 *
	 * device selector struct, it return the correct data type for each device
	 *
	 */
	template<typename T>
	struct device_cpu
	{
		//! cpu
		typedef typename device_cpu_impl<T,vect_isel<T>::value>::type type;
	};

	/*! device selector struct
	 *
	 * device selector struct, it return the correct data type for each device
	 *
	 */
	template<typename T>
	struct device_gpu
	{
		//! gpu
		typedef grid_gpu<1,T> type;
	};


	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * Stub object look at the various implementations
	 *
	 * \snippet vector_unit_tests.hpp Create add and access
	 *
	 * \param T type of structure the vector has to store
	 * \param device type of layout to use
	 * \param Memory allocator to use
	 * \param grow_p grow policy, how this vector should grow
	 *
	 * \see vector<T,device_cpu<T>,HeapMemory,grow_policy_double,STD_VECTOR>
	 * \see vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE>
	 * \see vector<T,device_gpu<T>,Memory,grow_p,OPENFPM_NATIVE>
	 *
	 */
	template<typename T, typename device=device_cpu<T>, typename Memory=HeapMemory, typename grow_p=grow_policy_double, unsigned int impl=vect_isel<T>::value>
	class vector
	{
	};

	#include "map_vector_std.hpp"


	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * The layout is memory_traits_lin
	 *
	 * ### Add and access elements
	 * \snippet vector_unit_tests.hpp Create add and access
	 *
	 * \tparam T type of object the vector store
	 * \tparam base memory layout to use
	 * \tparam Memory allocator to use
	 * \tparam grow_p grow policy for vector in case of reallocation
	 *
	 * OPENFPM_NATIVE implementation
	 *
	 */
	template<typename T,typename Memory, typename grow_p>
	class vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE>
	{
		//! This structure use this layout
		typedef typename grid_cpu<1,T>::memory_lin memory_lin;

		//! Actual size of the vector, warning: it is not the space allocated in grid
		//! grid size increase by a fixed amount every time we need a vector bigger than
		//! the actually allocated space
		size_t v_size;

		//! 1-D static grid
		grid_cpu<1,T> base;

		/*! \brief Create a 1D vector that contain the vector size
		 *
		 * Used to construct the underline 1D-grid
		 *
		 * \param sz size of the vector
		 *
		 * \return a vector with one element
		 *
		 */
		std::vector<size_t> getV(size_t sz)
		{
			std::vector<size_t> tmp;

			tmp.push_back(sz);

			return tmp;
		}

		/*! \brief If the argument is zero return 1 otherwise return the argument
		 *
		 * \param sz output
		 * \param arg argument
		 *
		 */
		void non_zero_one(size_t sz[1], size_t arg)
		{
			if (arg == 0)
			{sz[0] = 1;}
			else
			{sz[0] = arg;}
		}

	public:

		//! Type of the encapsulation memory parameter
		typedef typename grid_cpu<1,T>::memory_t memory_t;

		//! iterator for the vector
		typedef vector_key_iterator iterator_key;

		//! Object container for T, it is the return type of get_o it return a object type trough
		// you can access all the properties of T
		typedef typename grid_cpu<1,T>::container container;

		//! Type of the value the vector is storing
		typedef T value_type;

		/*! \brief Return the size of the vector
		 *
		 * \return the size
		 *
		 */
		size_t size() const
		{
			return v_size;
		}

		/*! \brief Reserve slots in the vector to avoid reallocation
		 *
		 * Reserve slots in the vector to avoid reallocation
		 *
		 * \param sp number of slot to reserve
		 *
		 */

		void reserve(size_t sp)
		{
			if (sp > base.size())
			{
				//! Resize the memory
				size_t sz[1] = {sp};
				base.template resize<Memory>(sz);
			}
		}

		/*! \brief Resize the vector
		 *
		 * Resize the vector and allocate n elements
		 *
		 * \param slot number of elements
		 *
		 */
		void resize(size_t slot)
		{
			// If we need more space than what we allocated, allocate new memory

			if (slot > base.size())
			{
				size_t gr = grow_p::grow(base.size(),slot);

				//! Resize the memory
				size_t sz[1] = {gr};
				base.template resize<Memory>(sz);
			}

			// update the vector size
			v_size = slot;
		}

		//! Access key for the vector
		typedef size_t access_key;

		/*! \brief It insert a new emtpy object on the vector, eventually it reallocate the grid
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 */
		void add()
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.template resize<Memory>(sz);
			}

			//! increase the vector size
			v_size++;
		}

		/*! \brief It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * \param v element to add
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 */
		void add(const T & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.template resize<Memory>(sz);
			}

			//! copy the element
			base.set(v_size,v);

			//! increase the vector size
			v_size++;
		}

		/*! \brief It insert a new object on the vector, eventually it reallocate the vector
		 *
		 * \param v object (encapsulated)
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 *
		 */
		void add(const typename grid_cpu<1,T>::container & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.template resize<Memory>(sz);
			}

			//! copy the added element
			base.set(v_size,v);

			//! increase the vector size
			v_size++;
		}

		/*! \brief Remove one entry from the vector
		 *
		 * \param key element to remove
		 *
		 */
		void remove(size_t key)
		{
			size_t d_k = key;
			size_t s_k = key + 1;

			// keys
			while (s_k < size())
			{
				set(d_k,get(s_k));
				d_k++;
				s_k++;
			}

			// re-calculate the vector size

			v_size--;
		}

		/*! \brief Remove several entries from the vector
		 *
		 * \param keys objects id to remove
		 * \param start key starting point
		 *
		 */
		void remove(openfpm::vector<size_t> keys, size_t start = 0)
		{
			// Nothing to remove return
			if (keys.size() <= start )
				return;

			size_t a_key = start;
			size_t d_k = keys.get(a_key);
			size_t s_k = keys.get(a_key) + 1;

			// keys
			while (s_k < size())
			{
				// s_k should always point to a key that is not going to be deleted
				while (a_key+1 < keys.size() && s_k == keys.get(a_key+1))
				{
					a_key++;
					s_k = keys.get(a_key) + 1;
				}

				// In case of overflow
				if (s_k >= size())
					break;

				set(d_k,get(s_k));
				d_k++;
				s_k++;
			}

			// re-calculate the vector size

			v_size -= keys.size() - start;
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \tparam p Property to get
		 * \param id Element to get
		 *
		 * \return the element value requested
		 *
		 */

		template <unsigned int p>inline typename type_cpu_prop<p,memory_lin>::type & get(size_t id)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " overflow" <<  "\n";}
#endif
			grid_key_dx<1> key(id);

			return base.template get<p>(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */

		inline auto get(size_t id) -> decltype(base.get_o(grid_key_dx<1>()))
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error "  << __FILE__ << "  " << __LINE__ << " vector overflow" << "\n";}
#endif
			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Constructor from vector
		 *
		 * \param mv vector
		 *
		 */
		vector(vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE> && mv)
		:v_size(mv.v_size),base(mv.base)
		{
		}

		/*! \brief It duplicate the vector
		 *
		 * \return a duplicated vector
		 *
		 */
		vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE> duplicate()
		{
			vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE> dup;

			dup.v_size = v_size;
			dup.base.swap(base.template duplicate<Memory>());

			return dup;
		}

		/*! \brief Constructor from another vector
		 *
		 * \param v the vector
		 *
		 */

		vector(const vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE> & v)
		:v_size(v.v_size),base(v.base)
		{}

		//! Constructor, vector of size 0
		vector():v_size(0),base(getV(0))
		{
			base.template setMemory<Memory>();
		}

		//! Constructor, vector of size sz
		vector(size_t sz):v_size(sz),base(getV(sz))
		{
			base.template setMemory<Memory>();
		}

		/*! \brief Set the object id to obj
		 *
		 * \param id element
		 * \param obj object (encapsulated)
		 *
		 */
		void set(size_t id, const typename grid_cpu<1,T>::container & obj)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif
			//! copy the element
			base.set(id,obj);
		}

		/*! \brief Set the object id to obj
		 *
		 * \param id
		 * \param obj
		 *
		 */
		void set(size_t id, T & obj)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif
			//! copy the element
			base.set(id,obj);
		}

		/*! \brief Set the element of the vector v from another element of another vector
		 *
		 * \param id element id
		 * \param v vector source
		 * \param src source element
		 *
		 */
		void set(size_t id, vector<T,device_cpu<T>,Memory,grow_p,OPENFPM_NATIVE> & v, size_t src)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif
			base.set(id,v.base,src);
		}

		/*! \brief Assignment operator
		 *
		 * move semantic movement operator=
		 *
		 * \param mv vector
		 *
		 * \return itself
		 *
		 */
		vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE> operator=(vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE> && mv)
		{
			v_size = mv.v_size;
			base.swap(mv.base);

			return *this;
		}

		/*! \brief Assignment operator
		 *
		 * it copy
		 *
		 * \param mv vector
		 *
		 * \return itself
		 *
		 */
		vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE> operator=(const vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE> & mv)
		{
			v_size = mv.v_size;
			size_t rsz[1] = {v_size};
			base.template resize<Memory>(rsz);

			// copy the object
			for (size_t i = 0 ; i < v_size ; i++ )
			{
				grid_key_dx<1> key(i);
				base.set(key,mv.base,key);
			}

			return *this;
		}

		/*! \brief Swap the memory with another vector
		 *
		 * \param v vector
		 *
		 */
		void swap(openfpm::vector<T,device_cpu<T>,Memory,grow_p,OPENFPM_NATIVE> & v)
		{
			size_t sz_sp = v_size;

			// swap the v_size
			v_size = v.v_size;

			base.swap(v.base);
			v.v_size = sz_sp;
		}

		/*! \brief Swap the memory with another vector
		 *
		 * \param v vector
		 *
		 */
		void swap(openfpm::vector<T,device_cpu<T>,Memory,grow_p,OPENFPM_NATIVE> && v)
		{
			size_t sz_sp = v_size;

			// swap the v_size
			v_size = v.v_size;

			base.swap(v.base);
			v.v_size = sz_sp;
		}

		/*! \brief Get the vector elements iterator
		 *
		 * Get the vector elements iterator
		 *
		 * \return an iterator to iterate through all the elements of the vector
		 *
		 */

		vector_key_iterator getIterator()
		{
			return vector_key_iterator(v_size);
		}

		/*! \brief Return the size of the message needed to pack this object
		 *
		 * \return The size
		 *
		 */

		size_t packObjectSize()
		{
			return base.packObjectSize();
		}

		/*! \brief Pack the object into the given pointer
		 *
		 * \param mem pointer
		 *
		 * \return the size of the packed message
		 *
		 */
		size_t packObject(void * mem)
		{
			return base.packObject(mem);
		}

		/*! \brief Calculate the memory size required to allocate n elements
		 *
		 * Calculate the total size required to store n-elements in a vector
		 *
		 * \param n number of elements
		 * \param e unused
		 *
		 * \return the size of the allocation number e
		 *
		 */
		inline static size_t calculateMem(size_t n, size_t e)
		{
			return grow_p::grow(0,n) * sizeof(T);
		}

		/*! \brief How many allocation are required to create n-elements
		 *
		 * \param n number of elements
		 *
		 * \return the number of allocations
		 *
		 */
		inline static size_t calculateNMem(size_t n)
		{
			return 1;
		}

		/*! \brief Set the memory of the base structure using an object
		 *
		 * \param mem Memory object to use for allocation
		 *
		 */
		void setMemory(Memory & mem)
		{
			base.template setMemory<Memory>(mem);
		}

		/*! \brief Return the pointer that store the data
		 *
		 * \return the pointer that store the data
		 *
		 */
		void * getPointer()
		{
			return base.getPointer();
		}
	};


	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * The layout is memory_traits_inte
	 *
	 * ### Add and access elements
	 * \snippet vector_unit_tests.hpp Create add and access
	 *
	 * \tparam T type of object the vector store
	 * \tparam base memory layout to use
	 * \tparam Memory allocator to use
	 * \tparam grow_p grow policy for vector in case of reallocation
	 *
	 * OPENFPM_NATIVE implementation
	 *
	 */
	template<typename T,typename Memory,typename grow_p>
	class vector<T,device_gpu<T>,Memory,grow_p,OPENFPM_NATIVE>
	{
		//! Actual size of the vector, warning: it is not the space allocated in grid
		//! grid size increase by a fixed amount every time we need a vector bigger than
		//! the actually allocated space
		size_t v_size;

		//! 1-D static grid
		grid_gpu<1,T> base;

	public:

		//! Type of the encapsulation memory parameter
		typedef typename grid_gpu<1,T>::memory_t memory_t;

		//! iterator for the vector
		typedef vector_key_iterator iterator_key;

		//! return the size of the vector
		size_t size()
		{
			return v_size;
		}

		/*! \brief Reserve slots in the vector to avoid reallocation
		 *
		 * Reserve slots in the vector to avoid reallocation
		 *
		 * \param sp number of slot to reserve
		 *
		 */

		void reserve(size_t sp)
		{
			// If we need more space than what we allocated allocate new memory

			if (sp > base.size())
			{
				//! Resize the memory
				std::vector<size_t> sz;
				sz.push_back(sp);
				base.resize(sz);
			}
		}

		/*! \brief Resize the vector
		 *
		 * Resize the vector to the specified number of elements
		 *
		 * \param slot number of elements
		 *
		 */
		void resize(size_t slot)
		{
			// If we need more space than what we allocated, allocate new memory

			if (slot > base.size())
			{
				//! Resize the memory
				std::vector<size_t> sz;
				sz.push_back(slot);
				base.resize(sz);
			}

			// update the vector size
			v_size = slot;
		}

		/*! \brief It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * \param v element to add
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 *
		 */
		void add(T & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory
				size_t sz[1] = {2*base.size()};
				base.resize(sz);
			}

			// copy the element
			base.set(v_size,v);
		}

		/*! \brief Set the object id to obj
		 *
		 * \param id element id
		 * \param obj object
		 *
		 */
		void set(size_t id, T & obj)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error "  << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif

			//! copy the element
			base.set(id,obj);
		}

		/*! \brief Set the element of the vector v from another element of another vector
		 *
		 * \param id element
		 * \param v source vector
		 * \param s_id source element id
		 *
		 */
		void set(size_t id, vector<T,device_gpu<T>,Memory,grow_p,OPENFPM_NATIVE> & v, size_t s_id)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error "  << __FILE__ << "  " << __LINE__ <<  " id overflow your vector" << "\n";}
#endif
			base.set(id,v.base,s_id);
		}

		/*! \brief Get the vector elements iterator
		 *
		 * \return an iterator to iterate through all the elements of the vector
		 *
		 */
		auto getIterator() -> decltype(base.getIterator())
		{
			return base.getIterator();
		}

		/*! \brief Calculate the memory size required to allocate n elements
		 *
		 * Calculate the total size required to store n-elements in a vector
		 *
		 * \param n number of elements
		 * \param e unused
		 *
		 * \return the size of the allocation number e
		 *
		 */
		inline static size_t calculateMem(size_t n, size_t e)
		{
			return n*sizeof(T);
		}

		/*! \brief How many allocation are required to create n-elements
		 *
		 * \param n number of elements
		 *
		 * \return the number of allocations
		 *
		 */
		inline static size_t calculateNMem(size_t n)
		{
			return 1;
		}

		/*! \brief Set the memory of the base structure using an object
		 *
		 * \param mem Memory object to use for allocation
		 *
		 */
		void setMemory(Memory & mem)
		{
			base.template setMemory<Memory>(mem);
		}
	};
}

#endif
