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

#define PAGE_ALLOC 1024

namespace openfpm
{
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
			size_t grow = original;
			while (grow < requested)	{grow *= 2;}
			return grow;
		}
	};

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

	/*! \brief device selector struct
	 *
	 * device selector struct, it return the correct data type for each device
	 *
	 */

	template<typename T>
	struct device_cpu
	{
		//! cpu
		typedef grid_cpu<1,T> type;
	};

	/*! \brief specialization for size_t
	 *
	 */

	template<>
	struct device_cpu<size_t>
	{};

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

	template<typename T, typename device=device_cpu<T>, typename Memory=HeapMemory, typename grow_p=grow_policy_double>
	class vector
	{
	};

	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * this implementation is just a wrapper for the std::vector in the case
	 * of the primitive size_t
	 *
	 * \param T base type
	 *
	 */

	template<typename Memory>
	class vector<size_t,device_cpu<size_t>,Memory>
	{
		//! Actual size of the vector, warning: it is not the space allocated in grid
		//! grid size increase by a fixed amount every time we need a vector bigger than
		//! the actually allocated space
		size_t v_size;

		//! 1-D static grid
		std::vector<size_t> base;

	public:

		// iterator for the vector
		typedef vector_key_iterator iterator_key;

		//! return the size of the vector
		inline size_t size()
		{
			return base.size();
		}


		/*! \ brief Resize the vector
		 *
		 * Resize the vector
		 *
		 * \param how many slot to reserve
		 *
		 */

		inline void resize(size_t slot)
		{
			base.resize(slot);
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
		inline void add(const size_t & v)
		{
			base.push_back(v);
		}

		/*! \brief Duplicate the vector
		 *
		 * \return the duplicated vector
		 *
		 */

		std::vector<size_t> duplicate()
		{
			return base;
		}

		/*! \brief swap the memory between the two vector
		 *
		 * \param vector to swap
		 *
		 */

		void swap(std::vector<size_t> && v)
		{
			base.swap(v);
		}


		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \tparam must be 0
		 *
		 * \param id Element to get
		 * \param p Property to get
		 *
		 */
		template <unsigned int p>inline size_t & get(size_t id)
		{
#ifdef DEBUG
			if (p != 0)
			{std::cerr << "Error the property does not exist" << "\n";}
#endif

			return base[id];
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \param id Element to get
		 * \param p Property to get
		 *
		 */
		inline size_t & get(size_t id)
		{
			return base[id];
		}

		/*! \brief it fill all the memory of fl patterns
		 *
		 * WARNING does not assign a value to each element but it fill the memory
		 * Useful to fast set the memory to zero
		 *
		 * \param fl byte to fill
		 *
		 */

		inline void fill(unsigned char fl)
		{
			memset(&base[0],0,base.size());
		}

		/*! \brief reserve a memory space in advance to avoid reallocation
		 *
		 * \param ns number of element the memory has to store
		 *
		 */

		inline void reserve(size_t ns)
		{
			base.reserve(ns);
		}

		//! Constructor, vector of size 0
		vector() {}

		//! Constructor, vector of size sz
		vector(size_t sz):base(sz) {}

		/*! swap the content of the vector
		 *
		 * \param vector to be swapped with
		 *
		 */

		void swap(openfpm::vector<size_t,device_cpu<size_t>,Memory> & v)
		{
			base.swap(v.base);
		}
	};

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
	class vector<T,device_cpu<T>, Memory,grow_p>
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

		// iterator for the vector
		typedef vector_key_iterator iterator_key;

		//! Object container for T, it is the return type of get_o it return a object type trough
		// you can access all the properties of T
		typedef typename grid_cpu<1,T>::container container;

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
				//! Resize the memory, double up the actual memory allocated for the vector
				std::vector<size_t> sz;
				sz.push_back(2*base.size());
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
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector file" <<  "\n";}
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

		/*! \brief It duplicate the vector
		 *
		 * \return a duplicated vector
		 *
		 */
		vector<T,device_cpu<T>, Memory,grow_p> duplicate()
		{
			vector<T,device_cpu<T>, Memory,grow_p> dup;

			dup.v_size = v_size;
			dup.base.swap(base.template duplicate<Memory>());

			return dup;
		}

		/*! \brief Constructor for a temporal object
		 *
		 * \param the temporal object
		 *
		 */

		vector(vector<T,device_cpu<T>, Memory,grow_p> && v)
		:v_size(v.v_size),base(v.base)
		{}

		//! Constructor, vector of size 0
		vector():v_size(0),base(getV(PAGE_ALLOC))
		{
			base.template setMemory<Memory>();
		}

		//! Constructor, vector of size sz
		vector(size_t sz):v_size(sz),base(getV(PAGE_ALLOC))
		{
			base.template setMemory<Memory>();
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
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif
			//! copy the element
			base.set(id,obj);
		}

		/* \brief Set the element of the vector v from another element of another vector
		 *
		 * Set the element of the vector v from another element of another vector
		 */

		void set(size_t id, vector<T,device_cpu<T>,Memory> & v, size_t src)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif
			base.set(id,v.base,src);
		}


		/*! \brief Swap the memory of another vector
		 *
		 * Swap the memory of another vector
		 *
		 */
		void swap(openfpm::vector<T,device_cpu<T>,Memory> & v)
		{
			base.swap(v.base);
		}

		/*! \brief Swap the memory of another vector
		 *
		 * Swap the memory of another vector
		 *
		 */
		void swap(openfpm::vector<T,device_cpu<T>,Memory> && v)
		{
			base.swap(v.base);
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
	 *
	 */

	template<typename T,typename Memory>
	class vector<T,device_gpu<T>,Memory>
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
				std::vector<size_t> sz;
				sz.push_back(base.size()*2);
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

		void set(size_t id, vector<T,device_gpu<T>,Memory> & v, size_t src)
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

	};
}

#endif
