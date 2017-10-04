/*
 * grid_gpu.hpp
 *
 *  Created on: Oct 31, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_GPU_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_GPU_HPP_

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */

template<typename T_type>
struct copy_memory_c
{
	//! encapsulated source object
	const typename memory_traits_inte<T_type>::type & src;
	//! encapsulated destination object
	typename memory_traits_inte_red<T_type>::type & dst;


	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline copy_memory_c(const typename memory_traits_inte<T_type>::type & src,
			                   typename memory_traits_inte_red<T_type>::type & dst)
	:src(src),dst(dst)
	{
	};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		boost::fusion::at_c<T::value>(dst).mem_r = boost::fusion::at_c<T::value>(src).mem_r;
	}
};

struct dim3_
{
	//! size in x dimension
	unsigned int x;

	//! size in y dimension
	unsigned int y;

	//! size in z dimension
	unsigned int z;
};

template<unsigned int dim>
struct device_grid
{
	//! number of treads in each block
	dim3_ threads;

	//! number of grid for the kernel execution
	dim3_ grids;
};


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
template<unsigned int dim, typename T, typename S>
class grid_cpu<dim,T,S,typename memory_traits_inte<T>::type> : public grid_base_impl<dim,T,S,typename memory_traits_inte<T>::type, memory_traits_inte>
{
	//! grid layout
	typedef typename memory_traits_inte<T>::type layout;

public:

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef typename grid_base_impl<dim,T,S,typename memory_traits_inte<T>::type, memory_traits_inte>::container container;

	//! Default constructor
	inline grid_cpu() THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_inte>()
	{
	}

	/*! \brief create a grid from another grid
	 *
	 * \param g the grid to copy
	 *
	 */
	inline grid_cpu(const grid_cpu & g) THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_inte>(g)
	{
	}

	/*! \brief create a grid of size sz on each direction
	 *
	 * \param sz grid size in each direction
	 *
	 */
	inline grid_cpu(const size_t & sz) THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_inte>(sz)
	{
	}

	//! Constructor allocate memory and give them a representation
	inline grid_cpu(const size_t (& sz)[dim]) THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_inte>(sz)
	{
	}

	/*! \brief It return the properties arrays.
	 *
	 * In case of Cuda memory it return the device pointers to pass to the kernels
	 *
	 *
	 */
	template<unsigned int id> void * getDeviceBuffer()
	{
		return boost::fusion::at_c<id>(this->data_).mem->getDevicePointer();
	}

	/*! \brief Synchronize the memory buffer in the device with the memory in the host
	 *
	 *
	 */
	template<unsigned int id> void deviceToHost()
	{
		return boost::fusion::at_c<id>(this->data_).mem->deviceToHost();
	}

	/*! \brief Convert the grid into a data-structure compatible for computing into GPU
	 *
	 *  The object created can be considered like a reference of the original
	 *
	 */
	grid_cpu<dim,T,S,typename memory_traits_inte<T>::type> toGPU()
	{
		grid_cpu<dim,T,S,typename memory_traits_inte<T>::type> g;
		copy_memory_c<T> cp_mc(this->data_,g.data_);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp_mc);

		return g;
	}

	/*! \brief When we switch to GPU mode the data-structure can be safely given to
	 *         a kernel for computation
	 *
	 *
	 */
	void switchToGPU()
	{

	}

	void switchToCPU()
	{

	}
};

//! short formula for a grid on gpu
template <unsigned int dim, typename T> using grid_gpu = grid_cpu<dim,T,CudaMemory,typename memory_traits_inte<T>::type>;


#endif /* OPENFPM_DATA_SRC_GRID_GRID_GPU_HPP_ */
