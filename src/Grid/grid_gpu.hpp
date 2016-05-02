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
 * element of the boost::vector the operator() is called
 *
 * \param T Type of memory allocator
 *
 */

template<typename S>
struct allocate
{
	//! size to allocate
	size_t sz;

	//! constructor it fix the size
	allocate(size_t sz)
	:sz(sz){};

	//! It call the allocate function for each member
	template<typename T>
	void operator()(T& t) const
	{
		//! Create and set the memory allocator
		t.setMemory(*new S());

		//! Allocate the memory and create the reppresentation
		t.allocate(sz);
	}
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
class grid_cpu<dim,T,S,typename memory_traits_inte<T>::type>
{
	typedef typename memory_traits_inte<T>::type layout;

	// Indicate if set memory has been called
	bool is_mem_init = false;

	//! Access the key
	typedef grid_key_dx<dim> access_key;

	//! It store all the information regarding the grid
	grid_sm<dim,void> g1;

	//! This is the interface to allocate,resize ... memory
	//! and give also a representation to the allocated memory
	layout data_;

public:

	//! it define that it is a grid
	typedef int yes_i_am_grid;

	//! Definition of the layout
	typedef typename memory_traits_inte<T>::type memory_int;

	//! Memory traits
	typedef layout memory_conf;

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef encapc<dim,T,layout> container;

	// The object type the grid is storing
	typedef T type;

	//! Default constructor
	grid_cpu() THROW
	{
	}

	//! Set the grid dimensions
	void setDimensions(size_t (& sz)[dim])
	{
		g1.setDimension(sz);
	}

	//! Constructor it initialize the memory and give representation
	grid_cpu(size_t sz) THROW
	:g1(sz)
	{
	}

	//! Constructor it initialize the memory and give representation
	grid_cpu(size_t (& sz)[dim]) THROW
	:g1(sz)
	{
	}

	/*! \brief Return the internal grid information
	 *
	 * \return the internal grid information
	 *
	 */

	const grid_sm<dim,void> & getGrid() const
	{
		return g1;
	}

	/*! \brief Create the object that provide memory
	 *
	 * \tparam S memory object type
	 *
	 */

	void setMemory()
	{
		//! Create an allocate object
		allocate<S> all(g1.size());

		//! for each element in the vector allocate the buffer
		boost::fusion::for_each(data_,all);

		is_mem_init = true;
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference to the element
	 *
	 */
	template <unsigned int p>inline typename type_gpu_prop<p,memory_int>::type::reference get(const grid_key_d<dim,p> & v1)
	{
		return boost::fusion::at_c<p>(data_).mem_r->operator[](g1.LinId(v1));
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference to the element
	 *
	 */
	template <unsigned int p>inline typename type_gpu_prop<p,memory_int>::type::reference get(const grid_key_dx<dim> & v1)
	{
		return boost::fusion::at_c<p>(data_).mem_r->operator[](g1.LinId(v1));
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference to the element
	 *
	 */
	template <unsigned int p>inline const typename type_gpu_prop<p,memory_int>::type::reference get(const grid_key_d<dim,p> & v1) const
	{
		return boost::fusion::at_c<p>(data_).mem_r->operator[](g1.LinId(v1));
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference to the element
	 *
	 */
	template <unsigned int p>inline const typename type_gpu_prop<p,memory_int>::type::reference get(const grid_key_dx<dim> & v1) const
	{
		return boost::fusion::at_c<p>(data_).mem_r->operator[](g1.LinId(v1));
	}

	/*! \brief Get the of the selected element as a boost::fusion::vector
	 *
	 * Get the selected element as a boost::fusion::vector
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 */
	inline encapc<dim,T,layout> get_o(const grid_key_dx<dim> & v1)
	{
		return encapc<dim,T,layout>(data_,g1.LinId(v1));
	}

	/*! \brief Get the of the selected element as a boost::fusion::vector
	 *
	 * Get the selected element as a boost::fusion::vector
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 */
	inline const encapc<dim,T,layout> get_o(const grid_key_dx<dim> & v1) const
	{
		return encapc<dim,T,layout>(data_,g1.LinId(v1));
	}

	/*! \brief Return the size of the grid
	 *
	 * \return the size of the grid
	 *
	 */
	inline size_t size()
	{
		return g1.size();
	}

	//! this function set the memory interface if required
	//! this operation is required when we define a void memory
	//! allocator
	void set_memory(memory & mem)
	{
		data_.mem.set_memory(mem);
	}

	/*! \brief Return a grid iterator
	 *
	 * Return a grid iterator, to iterate through the grid
	 *
	 */

	inline grid_key_dx_iterator<dim> getIterator()
	{
		return grid_key_dx_iterator<dim>(g1);
	}


	/*! \brief Return a sub-grid iterator
	 *
	 * Return a sub-grid iterator, to iterate through the grid
	 *
	 */

	inline grid_key_dx_iterator_sub<dim> getSubIterator(grid_key_dx<dim> & start, grid_key_dx<dim> & stop)
	{
		return grid_key_dx_iterator_sub<dim>(g1,start,stop);
	}

	/*! \brief Swap the memory of another grid
	 *
	 * Swap the memory of another grid
	 *
	 * \param obj Memory to swap with
	 *
	 */
	void swap(grid_cpu<dim,T,S,layout> & obj)
	{
		g1.swap(obj.g1);
		data_.swap(obj.data_);
	}
};

//! short formula for a grid on gpu
template <unsigned int dim, typename T> using grid_gpu = grid_cpu<dim,T,CudaMemory,typename memory_traits_inte<T>::type>;


#endif /* OPENFPM_DATA_SRC_GRID_GRID_GPU_HPP_ */
