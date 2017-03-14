/*
 * grid_base_implementation.hpp
 *
 *  Created on: May 2, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPLEMENTATION_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPLEMENTATION_HPP_

#include "grid_base_impl_layout.hpp"

/*! \brief
 *
 * Implementation of a N-dimensional grid
 *
 * \tparam dim dimansionality of the grid
 * \tparam T type store by the grid
 * \tparam S Memory pool from where to take the memory
 * \tparam layout_ memory layout
 * \tpaeam layout_base layout memory meta-function (the meta-function used to construct layout_)
 *
 */
template<unsigned int dim, typename T, typename S, typename layout_, template<typename> class layout_base >
class grid_base_impl
{
	//! memory layout
	typedef layout_ layout;

public:

	//! memory layout
	typedef layout_ layout_type;

	//! expose the dimansionality as a static const
	static constexpr unsigned int dims = dim;

	//! Access key
	typedef grid_key_dx<dim> access_key;

	//! boost::vector that describe the data type
	typedef typename T::type T_type;

protected:

	//! Memory layout specification + memory chunk pointer
	layout data_;

	//! This is a structure that store all information related to the grid and how indexes are linearized
	grid_sm<dim,T> g1;

private:

	//! Is the memory initialized
	bool is_mem_init = false;

	//! The memory allocator is not internally created
	bool isExternal;

#ifdef SE_CLASS1

	/*! \brief Check that the key is inside the grid
	 *
	 *
	 */
	inline void check_init() const
	{
		if (is_mem_init == false)
		{
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " you must call SetMemory before access the grid\n";
			ACTION_ON_ERROR(GRID_ERROR_OBJECT);
		}
	}

	/*! \brief Check that the key is inside the grid
	 *
	 * \param key
	 *
	 */
	inline void check_bound(const grid_key_dx<dim> & v1) const
	{
		for (long int i = 0 ; i < dim ; i++)
		{
			if (v1.get(i) >= (long int)getGrid().size(i))
			{
				std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << v1.get(i) << " >= " << getGrid().size(i) << "\n";
				ACTION_ON_ERROR(GRID_ERROR_OBJECT);
			}
			else if (v1.get(i) < 0)
			{
				std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << v1.get(i) << " is negative " << "\n";
				ACTION_ON_ERROR(GRID_ERROR_OBJECT);
			}
		}
	}

	/*! \brief Check that the key is inside the grid
	 *
	 * check if key2 is inside the g grid boundary
	 *
	 * \param g grid
	 * \param key2
	 *
	 */
	template<typename Mem> inline void check_bound(const grid_base_impl<dim,T,Mem,layout,layout_base> & g,const grid_key_dx<dim> & key2) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (key2.get(i) >= (long int)g.getGrid().size(i))
			{
				std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << key2.get(i) << " >= " << g.getGrid().size(i) << "\n";
				ACTION_ON_ERROR(GRID_ERROR_OBJECT);
			}
			else if (key2.get(i) < 0)
			{
				std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << key2.get(i) << " is negative " << "\n";
				ACTION_ON_ERROR(GRID_ERROR_OBJECT);
			}
		}
	}

#endif

public:

	// Implementation of packer and unpacker for grid
	#include "grid_pack_unpack.ipp"

	//! it define that it is a grid
	typedef int yes_i_am_grid;

	//! Definition of the layout
	typedef typename memory_traits_lin<typename T::type>::type memory_lin;

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef encapc<dim,T,layout> container;

	//! The object type the grid is storing
	typedef T value_type;

	//! Default constructor
	grid_base_impl() THROW
	:g1(0),isExternal(false)
	{
		// Add this pointer
#ifdef SE_CLASS2
		check_new(this,8,GRID_EVENT,1);
#endif
	}

	/*! \brief create a grid from another grid
	 *
	 *
	 * \param g the grid to copy
	 *
	 */
	grid_base_impl(const grid_base_impl & g) THROW
	:isExternal(false)
	{
		this->operator=(g);
	}

	/*! \brief create a grid of size sz on each direction
	 *
	 * \param sz size of the grid on each dimensions
	 *
	 */
	grid_base_impl(const size_t & sz) THROW
	:g1(sz),isExternal(false)
	{
		// Add this pointer
#ifdef SE_CLASS2
		check_new(this,8,GRID_EVENT,1);
#endif
	}

	//! Constructor allocate memory and give them a representation
	grid_base_impl(const size_t (& sz)[dim]) THROW
	:g1(sz),isExternal(false)
	{
		// Add this pointer
#ifdef SE_CLASS2
		check_new(this,8,GRID_EVENT,1);
#endif
	}

	//! Destructor
	~grid_base_impl() THROW
	{
		// delete this pointer
#ifdef SE_CLASS2
		check_delete(this);
#endif
	}

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 */
	grid_base_impl<dim,T,S,layout,layout_base> & operator=(const grid_base_impl<dim,T,S,layout,layout_base> & g)
	{
		// Add this pointer
#ifdef SE_CLASS2
		check_new(this,8,GRID_EVENT,1);
#endif
		swap(g.duplicate());

		return *this;
	}

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 * \return itself
	 *
	 */
	grid_base_impl<dim,T,S,layout,layout_base> & operator=(grid_base_impl<dim,T,S,layout,layout_base> && g)
	{
		// Add this pointer
#ifdef SE_CLASS2
		check_new(this,8,GRID_EVENT,1);
#endif

		swap(g);

		return *this;
	}

	/*! \brief Compare two grids
	 *
	 * \param g grid to check
	 *
	 * \return true if they match
	 *
	 */
	bool operator==(const grid_base_impl<dim,T,S,layout,layout_base> & g)
	{
		// check if the have the same size
		if (g1 != g.g1)
			return false;

		auto it = getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			if (this->get_o(key) != this->get_o(key))
				return false;

			++it;
		}

		return true;
	}

	/*! \brief create a duplicated version of the grid
	 *
	 * \return a duplicated version of the grid
	 *
	 */
	grid_base_impl<dim,T,S,layout,layout_base> duplicate() const THROW
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		//! Create a completely new grid with sz

		grid_base_impl<dim,T,S,layout,layout_base> grid_new(g1.getSize());

		//! Set the allocator and allocate the memory
		grid_new.setMemory();

		// We know that, if it is 1D we can safely copy the memory
//		if (dim == 1)
//		{
			//! 1-D copy (This case is simple we use raw memory copy because is the fastest option)
//			grid_new.data_.mem->copy(*data_.mem);
//		}
//		else
//		{
			//! N-D copy

			//! create a source grid iterator
			grid_key_dx_iterator<dim> it(g1);

			while(it.isNext())
			{
				grid_new.set(it.get(),*this,it.get());

				++it;
			}
//		}

		// copy grid_new to the base

		return grid_new;
	}

	/*! \brief Return the internal grid information
	 *
	 * Return the internal grid information
	 *
	 * \return the internal grid
	 *
	 */

	const grid_sm<dim,T> & getGrid() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return g1;
	}

	/*! \brief Create the object that provide memory
	 *
	 * Create the object that provide memory
	 *
	 * \tparam S memory type to allocate
	 *
	 */

	void setMemory()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif

		mem_setm<S,layout_base<T>,decltype(this->data_),decltype(this->g1)>::setMemory(data_,g1,is_mem_init);
	}

	/*! \brief Get the object that provide memory
	 *
	 * An external allocator is useful with allocator like PreAllocHeapMem
	 * to have contiguous in memory vectors. Or to force the system to retain
	 * memory
	 *
	 * \tparam S memory type
	 *
	 * \param m external memory allocator
	 *
	 */

	void setMemory(S & m)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		//! Is external
		isExternal = true;

		//! Create and set the memory allocator
		data_.setMemory(m);

		//! Allocate the memory and create the reppresentation
		if (g1.size() != 0) data_.allocate(g1.size());

		is_mem_init = true;
	}

	/*! \brief Return a plain pointer to the internal data
	 *
	 * Return a plain pointer to the internal data
	 *
	 * \return plain data pointer
	 *
	 */

	void * getPointer()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		if (data_.mem_r == NULL)
			return NULL;

		return data_.mem_r->get_pointer();
	}

	/*! \brief Return a plain pointer to the internal data
	 *
	 * Return a plain pointer to the internal data
	 *
	 * \return plain data pointer
	 *
	 */

	const void * getPointer() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		if (data_.mem_r == NULL)
			return NULL;

		return data_.mem_r->get_pointer();
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(mem_get<p,layout_base<T>,layout,grid_sm<dim,T>,grid_key_dx<dim>>::get(data_,g1,grid_key_dx<dim>()))> inline r_type get(const grid_key_dx<dim> & v1)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return mem_get<p,layout_base<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(data_,g1,v1);
	}

	/*! \brief Get the const reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the const reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(mem_get<p,layout_base<T>,layout,grid_sm<dim,T>,grid_key_dx<dim>>::get(data_,g1,grid_key_dx<dim>()))> inline const r_type get(const grid_key_dx<dim> & v1) const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return mem_get<p,layout_base<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(data_,g1,v1);
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
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return mem_geto<dim,T,layout_base<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(data_,g1,v1);
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
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return mem_geto<dim,T,layout_base<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(const_cast<decltype(this->data_) &>(data_),g1,v1);
	}

	/*! \brief Fill the memory with the selected byte
	 *
	 * \warning It is a low level memory operation it ignore any type and semantic safety
	 *
	 * \param fl byte pattern to fill
	 *
	 */
	void fill(unsigned char fl)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		memset(getPointer(),fl,size() * sizeof(T));
	}

	/*! \brief Resize the space
	 *
	 * Resize the space to a new grid, the element are retained on the new grid,
	 * if the new grid is bigger the new element are now initialized, if is smaller
	 * the data are cropped
	 *
	 * \param sz reference to an array of dimension dim
	 *
	 */
	void resize(const size_t (& sz)[dim])
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		//! Create a completely new grid with sz

		grid_base_impl<dim,T,S,layout,layout_base> grid_new(sz);

		//! Set the allocator and allocate the memory
		if (isExternal == true)
		{
/*			grid_new.setMemory(static_cast<S&>(data_.getMemory()));

			// Create an empty memory allocator for the actual structure

			setMemory();*/
			mem_setext<decltype(grid_new),S,layout_base<T>,decltype(data_)>::set(grid_new,*this,this->data_);
		}
		else
			grid_new.setMemory();


		// We know that, if it is 1D we can safely copy the memory
//		if (dim == 1)
//		{
//			//! 1-D copy (This case is simple we use raw memory copy because is the fastest option)
//			grid_new.data_.mem->copy(*data_.mem);
//		}
//		else
//		{
		// It should be better to separate between fast and slow cases

			//! N-D copy

			size_t sz_c[dim];
			for (size_t i = 0 ; i < dim ; i++)
				sz_c[i] = (g1.size(i) < sz[i])?g1.size(i):sz[i];

			grid_sm<dim,void> g1_c(sz_c);

			//! create a source grid iterator
			grid_key_dx_iterator<dim> it(g1_c);

			while(it.isNext())
			{
				// get the grid key
				grid_key_dx<dim> key = it.get();

				// create a copy element

				grid_new.get_o(key) = this->get_o(key);

				++it;
			}
//		}

		// copy grid_new to the base

		this->swap(grid_new);
	}

	/*! \brief Remove one element valid only on 1D
	 *
	 *
	 */
	void remove(size_t key)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		if (dim != 1)
		{
#ifdef SE_CLASS1
			std::cerr << "Error: " << __FILE__ << " " << __LINE__ << " remove work only on dimension == 1 " << "\n";
#endif
			return;
		}

		// It is safe to do a memory copy

		data_.move(&this->template get<0>());
	}


	/*! \brief It move the allocated object from one grid to another
	 *
	 * It move the allocated object from one grid to another, after this
	 * call the argument grid is no longer valid
	 *
	 * \param grid to move/copy
	 *
	 */

	void swap(grid_base_impl<dim,T,S,layout,layout_base> & grid)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		// move the data
//		data_.swap(grid.data_);

		mem_swap<T,layout_base<T>,decltype(data_),decltype(grid)>::swap(data_,grid.data_);

		// exchange the grid info
		g1.swap(grid.g1);

		// exchange the init status
		bool exg = is_mem_init;
		is_mem_init = grid.is_mem_init;
		grid.is_mem_init = exg;

		// exchange the is external status
		exg = isExternal;
		isExternal = grid.isExternal;
		grid.isExternal = exg;
	}

	/*! \brief It move the allocated object from one grid to another
	 *
	 * It move the allocated object from one grid to another, after this
	 * call the argument grid is no longer valid
	 *
	 * \param grid to move/copy
	 *
	 */

	void swap(grid_base_impl<dim,T,S,layout,layout_base> && grid)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		swap(grid);
	}

	/*! \brief set an element of the grid
	 *
	 * set an element of the grid
	 *
	 * \param dx is the grid key or the position to set
	 * \param obj value to set
	 *
	 */

	template<typename Memory> inline void set(grid_key_dx<dim> dx, const encapc<1,T,Memory> & obj)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		check_init();
		check_bound(dx);
#endif

		// create the object to copy the properties
		copy_cpu_encap<dim,grid_base_impl<dim,T,S,layout,layout_base>,layout> cp(dx,*this,obj);

		// copy each property
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(cp);
	}

	/*! \brief set an element of the grid
	 *
	 * set an element of the grid
	 *
	 * \param dx is the grid key or the position to set
	 * \param obj value to set
	 *
	 */

	inline void set(grid_key_dx<dim> dx, const T & obj)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		check_init();
		check_bound(dx);
#endif

		this->get_o(dx) = obj;
	}


	/*! \brief Set an element of the grid from another element of another grid
	 *
	 * \param key1 element of the grid to set
	 * \param g source grid
	 * \param key2 element of the source grid to copy
	 *
	 */

	inline void set(const grid_key_dx<dim> & key1,const grid_base_impl<dim,T,S,layout,layout_base> & g, const grid_key_dx<dim> & key2)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		check_init();
		check_bound(key1);
		check_bound(g,key2);
#endif

		this->get_o(key1) = g.get_o(key2);
	}

	/*! \brief Set an element of the grid from another element of another grid
	 *
	 * \param key1 element of the grid to set
	 * \param g source grid
	 * \param key2 element of the source grid to copy
	 *
	 */

	template<typename Mem> inline void set(const grid_key_dx<dim> & key1,const grid_base_impl<dim,T,Mem,layout,layout_base> & g, const grid_key_dx<dim> & key2)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
#ifdef SE_CLASS1
		check_init();
		check_bound(key1);
		check_bound(g,key2);
#endif

		this->get_o(key1) = g.get_o(key2);
	}

	/*! \brief return the size of the grid
	 *
	 * Return the size of the grid
	 *
	 */
	inline size_t size() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return g1.size();
	}

	/*! \brief Return a sub-grid iterator
	 *
	 * Return a sub-grid iterator, to iterate through the grid
	 *
	 * \param start start point
	 * \param stop stop point
	 *
	 * \return a sub-grid iterator
	 *
	 */
	inline grid_key_dx_iterator_sub<dim> getSubIterator(grid_key_dx<dim> & start, grid_key_dx<dim> & stop) const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return g1.getSubIterator(start,stop);
	}

	/*! \brief Return a sub-grid iterator
	 *
	 * Return a sub-grid iterator, to iterate through the grid
	 *
	 * \param m Margin
	 *
	 * \return a sub-grid iterator
	 *
	 */
	inline grid_key_dx_iterator_sub<dim> getSubIterator(size_t m)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return grid_key_dx_iterator_sub<dim>(g1,m);
	}

	/*! \brief Return a grid iterator
	 *
	 * Return a grid iterator, to iterate through the grid
	 *
	 * \return a grid iterator
	 *
	 */
	inline grid_key_dx_iterator<dim> getIterator() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return grid_key_dx_iterator<dim>(g1);
	}

	/*! \brief Return a grid iterator over all the point with the exception
	 *   of the ghost part
	 *
	 * Return a grid iterator over all the point with the exception of the
	 * ghost part
	 *
	 * \return a sub-grid iterator
	 *
	 */
	inline grid_key_dx_iterator_sub<dim> getIterator(const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop) const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		// get the starting point and the end point of the real domain

		return grid_key_dx_iterator_sub<dim>(g1,start,stop);
	}

	/*! \brief Return the size of the message needed to pack this object
	 *
	 * TODO They just return 0 for now
	 *
	 * \return The size of the object to pack this object
	 *
	 *
	 */

	size_t packObjectSize()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return 0;
	}

	/*! \brief It fill the message packet
	 *
	 * TODO They just return 0 doing nothing
	 *
	 * \return The packet size
	 *
	 *
	 */
	size_t packObject(void * mem)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return 0;
	}

	/* \brief It return the id of structure in the allocation list
	 *
	 * \see print_alloc and SE_CLASS2
	 *
	 */
	long int who()
	{
#ifdef SE_CLASS2
		return check_whoami(this,8);
#else
		return -1;
#endif
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPLEMENTATION_HPP_ */
