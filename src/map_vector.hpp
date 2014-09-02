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
	struct device_v
	{
		//! cpu
		typedef grid_cpu<1,T> cpu;
		//! gpu
		typedef grid_gpu<1,T> gpu;
	};


	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * Implementation of 1-D std::vector like structure, empty structure
	 * when I do not know how to specialize
	 *
	 */

	template<typename T, typename device=device_v<T>::cpu>
	class vector
	{
	};

	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * Implementation of 1-D std::vector like structure, the memory allocated
	 * increase by PAGE_ALLOC constant every time we need more space
	 *
	 * \param T base type
	 *
	 */

	template<typename T>
	class vector<T,device_v<T>::cpu>
	{
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
			//! Here we are doing computation, reallocation is not permitted

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

	template<typename T>
	class vector<T,device_v<T>::gpu>
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
			//! Here we are doing computation, reallocation is not permitted

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
