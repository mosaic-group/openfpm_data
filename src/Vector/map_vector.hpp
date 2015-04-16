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
		size_t end;

	protected:

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
		 * Set the dimension
		 *
		 * \param dim is the dimension
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

	// default grow policy
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

	/*! \brief Select the correct implementation based on the type checking
	 *
	 * case OPENFPM native
	 *
	 */

	template<typename T, unsigned int impl>
	struct device_cpu_impl
	{
		typedef grid_cpu<1,T> type;
	};

	/*! \brief Select the correct implementation based on the type checking
	 *
	 * case STD_VECTOR native
	 *
	 */

	template<typename T>
	struct device_cpu_impl<T,STD_VECTOR>
	{
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
	 * Implementation of 1-D std::vector like structure, empty structure
	 * when I do not know how to specialize, should be never selected by the
	 * compiler
	 *
	 * \param T type of structure the vector has to store
	 * \param device type of layout to use
	 * \param Memory allocator to use
	 * \param grow_p grow policy, how this vector should grow
	 *
	 */

	template<typename T, typename device=device_cpu<T>, typename Memory=HeapMemory, typename grow_p=grow_policy_double, unsigned int impl=vect_isel<T>::value>
	class vector
	{
	};

	#include "map_vector_std.hpp"


	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * Implementation of 1-D std::vector like structure, the memory allocated
	 * increase by PAGE_ALLOC constant every time we need more space
	 *
	 * \param T type of structure the vector has to store
	 * \param device type of layout to use
	 * \param Memory allocator to use
	 * \param grow_p grow policy for vector
	 *
	 */

	template<typename T,typename Memory, typename grow_p>
	class vector<T,device_cpu<T>, Memory,grow_p,OPENFPM_NATIVE>
	{
		//! Actual size of the vector, warning: it is not the space allocated in grid
		//! grid size increase by a fixed amount every time we need a vector bigger than
		//! the actually allocated space
		size_t v_size;

		//! 1-D static grid
		grid_cpu<1,T> base;


		/*! \brief Get 1D vector to create the grid in the constructor
		 *
		 * Get 1D vector to create the grid in the constructor
		 *
		 * \param size_t sizeof the vector
		 *
		 */

		std::vector<size_t> getV(size_t sz)
		{
			std::vector<size_t> tmp;

			tmp.push_back(sz);

			return tmp;
		}

	public:

		//! iterator for the vector
		typedef vector_key_iterator iterator_key;

		//! Object container for T, it is the return type of get_o it return a object type trough
		// you can access all the properties of T
		typedef typename grid_cpu<1,T>::container container;

		//! Type of the value the vector is storing
		typedef T value_type;

		/*! \brief Return the size of the vector
		 *
		 * Return the size of the vector
		 *
		 */
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
			if (sp > base.size())
			{
				//! Resize the memory
				size_t sz[1] = {sp};
				base.template resize<Memory>(sz);
			}
		}

		/*! \brief Resize the vector
		 *
		 * Resize the vector and allocate n slot
		 *
		 * \param resize the vector of n-slot
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

		// Access key for the vector
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
				size_t sz[1] = {2*base.size()};
				base.template resize<Memory>(sz);
			}

			//! increase the vector size
			v_size++;
		}

		/*! \brief It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 */
		void add(T & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				if (base.size() == 0)
				{sz[0] = 1;}
				else
				{sz[0] = 2*base.size();}
				base.template resize<Memory>(sz);
			}

			//! copy the element
			base.set(v_size,v);

			//! increase the vector size
			v_size++;
		}

		/*! \brief It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 *
		 */
		void add(typename grid_cpu<1,T>::container & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1] = {2*base.size()};
				base.template resize<Memory>(sz);
			}

			//! copy the element
			base.set(v_size,v);

			//! increase the vector size
			v_size++;
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \param id Element to get
		 * \param p Property to get
		 *
		 */

		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & get(size_t id)
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
		 * \param p Property to get
		 *
		 */

		inline auto get(size_t id) -> decltype(base.template get_o(grid_key_dx<1>()))
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error "  << __FILE__ << "  " << __LINE__ << " vector overflow" << "\n";}
#endif
			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/* \brief Constructor from vector
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

		/*! \brief Constructor for a temporal object
		 *
		 * \param the temporal object
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

		/*! \brief Set the object
		 *
		 * Set the object id to obj
		 *
		 * \param id
		 * \param object (encapsulated)
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

		/*! \brief Set the object
		 *
		 * Set the object id to obj
		 *
		 * \param id
		 * \param object
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

		/* \brief Set the element of the vector v from another element of another vector
		 *
		 * Set the element of the vector v from another element of another vector
		 */

		void set(size_t id, vector<T,device_cpu<T>,Memory,grow_p,OPENFPM_NATIVE> & v, size_t src)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif
			base.set(id,v.base,src);
		}

		/* \brief operator=
		 *
		 * move semantic movement operator=
		 *
		 * \param mv Vector operator
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

		/* \brief operator=
		 *
		 * Real copy
		 *
		 * \param mv Vector operator
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

		/*! \brief Swap the memory of another vector
		 *
		 * Swap the memory of another vector
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

		/*! \brief Swap the memory of another vector
		 *
		 * Swap the memory of another vector
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
		 * \return The size of the object to pack this object
		 *
		 */

		size_t packObjectSize()
		{
			return base.packObjectSize();
		}

		/*! \brief Return the size of the message needed to pack this object
		 *
		 * \return The size of the object to pack this object
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
	 * Implementation of 1-D std::vector like structure, the memory allocated
	 * increase by PAGE_ALLOC constant every time we need more space
	 *
	 * \param T base type
	 *
	 * \warning this implementation work on cpu and gpu, but on gpu its functionality
	 * is limited:
	 *
	 * 1) add is limited because reallocation cannot be handled on gpu,
	 *    a potential need for reallocation is signaled with an overflow
	 * 2) add is not thread safe so each thread on gpu should operate on a different
	 *    vector, or a thread at time should operate on the vector
	 *
	 * \param T type of structure the vector has to store
	 * \param device type of layout to use
	 * \param Memory allocator to use
	 * \param grow_p Grow policy of the vector
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

		// iterator for the vector
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
		 * Resize the vector and allocate n slot
		 *
		 * \param resize the vector of n-slot
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
		 * It insert a new object on the vector, eventually it reallocate the grid
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

		/*! \brief Set the object
		 *
		 * Set the object id to obj
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

		/* \brief Set the element of the vector v from another element of another vector
		 *
		 * Set the element of the vector v from another element of another vector
		 */

		void set(size_t id, vector<T,device_gpu<T>,Memory,grow_p,OPENFPM_NATIVE> & v, size_t src)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error "  << __FILE__ << "  " << __LINE__ <<  " id overflow your vector" << "\n";}
#endif
			base.set(id,v.base,src);
		}

		/*! \brief Get the vector elements iterator
		 *
		 * Get the vector elements iterator
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
