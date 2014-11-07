/*
 * map_vector.hpp
 *
 *  Created on: Aug 30, 2014
 *      Author: Pietro Incardona
 */

#ifndef MAP_VECTOR_HPP
#define MAP_VECTOR_HPP

#include "map_grid.hpp"

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
		typedef grid_cpu<1,T> cpu;
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
		typedef grid_gpu<1,T> gpu;
	};


	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * Implementation of 1-D std::vector like structure, empty structure
	 * when I do not know how to specialize, should be never selected by the
	 * compiler
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

		/*! \brief It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 *
		 */
		void push_back(size_t & v)
		{
			base.push_back(v);
		}

	public:

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
	 * \param T base type
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

		//! return the size of the vector
		size_t size()
		{
			return v_size;
		}


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

		/*! \brief It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 *
		 */
		void push_back(T & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory of PAGE_ALLOC elements
				std::vector<size_t> sz;
				sz.push_back(PAGE_ALLOC);
				base.template resize<Memory>(sz);
			}

			//! copy the element
			base.set(v_size,v);
		}

		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & get(size_t id)
		{
			grid_key_dx<1> key(id);

			return base.template get<p>(key);
		}

		//! Constructor, vector of size 0
		vector():v_size(0),base(getV(0)) {}

		//! Constructor, vector of size sz
		vector(size_t sz):v_size(sz),base(getV(sz)) {}

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
	 * 1) push_back is limited because reallocation cannot be handled on gpu,
	 *    a potential need for reallocation is signaled with an overflow
	 * 2) push_back is not thread safe so each thread on gpu should operate on a different
	 *    vector, or a thread at time should operate on the vector
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

		/*! \brief It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 *
		 */
		void push_back(T & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory
				std::vector<size_t> sz;
				sz.push_back(PAGE_ALLOC);
				base.resize(sz);
			}

			// copy the element
			base.set(v_size,v);
		}
	};
}

#endif
