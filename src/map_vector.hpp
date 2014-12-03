/*
 * map_vector.hpp
 *
 *  Created on: Aug 30, 2014
 *      Author: Pietro Incardona
 */

#ifndef MAP_VECTOR_HPP
#define MAP_VECTOR_HPP

#include "map_grid.hpp"
#include "memory/HeapMemory.hpp"

#define PAGE_ALLOC 1024

namespace openfpm
{

	/*! device selector struct
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
	 *
	 */

	template<typename T, typename device=device_cpu<T>, typename Memory=HeapMemory>
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
		//! Indicate if reallocation is needed on cpu is always false;
		bool need_reallocation;

		//! Actual size of the vector, warning: it is not the space allocated in grid
		//! grid size increase by a fixed amount every time we need a vector bigger than
		//! the actually allocated space
		size_t v_size;

		//! 1-D static grid
		std::vector<size_t> base;

		//! return the size of the vector
		size_t size()
		{
			return base.size();
		}

	public:

		/*! \ brief Resize the vector
		 *
		 * Resize the vector
		 *
		 * \param how many slot to reserve
		 *
		 */

		void resize(size_t slot)
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
		void add(size_t & v)
		{
			base.push_back(v);
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
		void add(size_t v)
		{
			base.push_back(v);
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
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

		//! Constructor, vector of size 0
		vector() {}

		//! Constructor, vector of size sz
		vector(size_t sz):base(sz) {}

	};

	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * Implementation of 1-D std::vector like structure, the memory allocated
	 * increase by PAGE_ALLOC constant every time we need more space
	 *
	 * \param T type of structure the vector has to store
	 * \param device type of layout to use
	 * \param Memory allocator to use
	 *
	 */

	template<typename T,typename Memory>
	class vector<T,device_cpu<T>, Memory>
	{
		//! Indicate if reallocation is needed on cpu is always false;
		bool need_reallocation;

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
				std::vector<size_t> sz;
				sz.push_back(sp);
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
				//! Resize the memory
				std::vector<size_t> sz;
				sz.push_back(slot);
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
			grid_key_dx<1> key(id);

			return base.template get<p>(key);
		}

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
			//! copy the element
			base.set(id,obj);
		}

		/* \brief Set the element of the vector v from another element of another vector
		 *
		 * Set the element of the vector v from another element of another vector
		 */

		void set(size_t id, vector<T,device_cpu<T>,Memory> & v, size_t src)
		{
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
			//! copy the element
			base.set(id,obj);
		}

		/* \brief Set the element of the vector v from another element of another vector
		 *
		 * Set the element of the vector v from another element of another vector
		 */

		void set(size_t id, vector<T,device_gpu<T>,Memory> & v, size_t src)
		{
			base.set(id,v.base,src);
		}
	};
}

#endif
