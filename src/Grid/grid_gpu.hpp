/*
 * grid_gpu.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_GPU_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_GPU_HPP_

/*! \brief This is an N-dimensional grid or an N-dimensional array with memory_traits_inte layout
 *
 * it is basically an N-dimensional Cartesian grid
 *
 *	\tparam dim Dimensionality of the grid
 *	\tparam T type of object the grid store
 *	\tparam Mem memory layout
 *
 * ### Definition and allocation of a 3D grid on GPU memory
 * \snippet grid_unit_tests.hpp Definition and allocation of a 3D grid on GPU memory
 * ### Access a grid c3 of size sz on each direction
 * \snippet grid_unit_tests.hpp Access a grid c3 of size sz on each direction
 * ### Access to an N-dimensional grid with an iterator
 * \snippet grid_unit_tests.hpp Access to an N-dimensional grid with an iterator
 *
 */
//template<unsigned int dim, typename T, typename S>
//class grid_cpu<dim,T,S,typename memory_traits_inte<T>::type> : public grid_base_impl<dim,T,S,typename memory_traits_inte<T>::type, memory_traits_inte>
//{
//	//! grid layout
//	typedef typename memory_traits_inte<T>::type layout;
//
//public:
//
//	//! Object container for T, it is the return type of get_o it return a object type trough
//	// you can access all the properties of T
//	typedef typename grid_base_impl<dim,T,S,typename memory_traits_inte<T>::type, memory_traits_inte>::container container;
//
//	//! Default constructor
//	inline grid_cpu() THROW
//	:grid_base_impl<dim,T,S,layout,memory_traits_inte>()
//	{
//	}
//
//	/*! \brief create a grid from another grid
//	 *
//	 * \param g the grid to copy
//	 *
//	 */
//	inline grid_cpu(const grid_cpu & g) THROW
//	:grid_base_impl<dim,T,S,layout,memory_traits_inte>(g)
//	{
//	}
//
//	/*! \brief create a grid of size sz on each direction
//	 *
//	 * \param sz grid size in each direction
//	 *
//	 */
//	inline grid_cpu(const size_t & sz) THROW
//	:grid_base_impl<dim,T,S,layout,memory_traits_inte>(sz)
//	{
//	}
//
//	//! Constructor allocate memory and give them a representation
//	inline grid_cpu(const size_t (& sz)[dim]) THROW
//	:grid_base_impl<dim,T,S,layout,memory_traits_inte>(sz)
//	{
//	}
//};
//
////! short formula for a grid on gpu
//template <unsigned int dim, typename T> using grid_gpu = grid_cpu<dim,T,CudaMemory,typename memory_traits_inte<T>::type>;


#endif /* OPENFPM_DATA_SRC_GRID_GRID_GPU_HPP_ */
