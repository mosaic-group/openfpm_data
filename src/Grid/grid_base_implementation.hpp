/*
 * grid_base_implementation.hpp
 *
 *  Created on: May 2, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPLEMENTATION_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPLEMENTATION_HPP_

#include "grid_base_impl_layout.hpp"

template<unsigned int dim, typename T, typename S, typename layout_, template<typename> class layout_base >
class grid_base_impl
{
	typedef layout_ layout;

public:

	typedef layout_ layout_type;

	// expose the dimansionality as a static const
	static constexpr unsigned int dims = dim;

	//! Access key
	typedef grid_key_dx<dim> access_key;

	//! boost::vector that describe the data type
	typedef typename T::type T_type;

	/*! \brief Pack the object into the memory given an iterator
	 *
	 * \tparam dim Dimensionality of the grid
	 * \tparam prp properties to pack
	 *
	 * \param mem preallocated memory where to pack the objects
	 * \param obj object to pack
	 * \param sts pack statistic
	 *
	 */

	static bool pack()
	{
		return false;
	}

	static bool packRequest()
	{
		return false;
	}

	static bool calculateMem()
	{
		return false;
	}

	template<int ... prp> void pack(ExtPreAlloc<S> & mem, Pack_stat & sts)
	{
#ifdef DEBUG
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif

		// Sending property object and vector
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,ExtPreAlloc<S>, typename layout_base<prp_object>::type,layout_base, openfpm::grow_policy_identity> dtype;
		dtype dvect;

		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(size());

		auto it = getIterator();

		pack_with_iterator<decltype(it),dtype,prp...>(it,dest);

		// Update statistic
		sts.incReq();
	}


	/*! \brief Pack the object into the memory given an iterator
	 *
	 * \tparam prp properties to pack
	 *
	 * \param mem preallocated memory where to pack the objects
	 * \param sub_it sub grid iterator ( or the elements in the grid to pack )
	 * \param sts pack statistic
	 *
	 */
	template<int ... prp> void pack(ExtPreAlloc<S> & mem, grid_key_dx_iterator_sub<dims> & sub_it, Pack_stat & sts)
	{
#ifdef DEBUG
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif

		// Sending property object
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,ExtPreAlloc<S>,typename layout_base<prp_object>::type, layout_base,openfpm::grow_policy_identity> dtype;

		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(sub_it.getVolume());

		pack_with_iterator<grid_key_dx_iterator_sub<dims>,dtype,prp...>(sub_it,dest);

		// Update statistic
		sts.incReq();
	}

	/*! \brief Insert an allocation request
	 *
	 * \param vector of requests
	 *
	 */
	template<int ... prp> void packRequest(std::vector<size_t> & v)
	{
		// Sending property object
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,ExtPreAlloc<S>,openfpm::grow_policy_identity> dtype;
		dtype dvect;

		// Calculate the required memory for packing
		size_t alloc_ele = dvect.calculateMem(size(),0);

		v.push_back(alloc_ele);
	}

	/*! \brief Insert an allocation request
	 *
	 * \tparam prp set of properties to pack
	 *

	 * \param sub sub-grid iterator
	 * \param vector of requests
	 *
	 */
	template<int ... prp> void packRequest(grid_key_dx_iterator_sub<dims> & sub, std::vector<size_t> & v)
	{
		typedef openfpm::vector<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type,ExtPreAlloc<S>,layout,layout_base,openfpm::grow_policy_identity> dtype;
		dtype dvect;

		// Calculate the required memory for packing
		size_t alloc_ele = dvect.template calculateMem<prp...>(sub.getVolume(),0);

		v.push_back(alloc_ele);
	}

	/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
	 *
	 * \tparam it type of iterator of the grid-structure
	 * \tparam dtype type of the structure B
	 * \tparam dim Dimensionality of the grid
	 * \tparam properties to pack
	 *
	 * \param it Grid iterator
	 * \param obj object to pack
	 * \param dest where to pack
	 *
	 */
	template <typename it, typename dtype, int ... prp> void pack_with_iterator(it & sub_it, dtype & dest)
	{
		// Sending property object
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;

		size_t id = 0;

		// Packing the information
		while (sub_it.isNext())
		{
			// copy all the object in the send buffer
			typedef encapc<dims,value_type,layout > encap_src;
			// destination object type
			typedef encapc<1,prp_object,typename dtype::layout_type > encap_dst;

			// Copy only the selected properties
			object_si_d<encap_src,encap_dst,OBJ_ENCAP,prp...>(get_o(sub_it.get()),dest.get(id));

			++id;
			++sub_it;
		}
	}

	/*! \brief unpack the grid given an iterator
	 *
	 * \tparam it type of iterator
	 * \tparam prp of the grid object to unpack
	 *
	 */
	template <typename it, typename stype, unsigned int ... prp> void unpack_with_iterator(ExtPreAlloc<S> & mem, it & sub_it, stype & src, Unpack_stat & ps)
	{
		size_t id = 0;

		// Sending property object
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;

		// unpacking the information
		while (sub_it.isNext())
		{
			// copy all the object in the send buffer
			typedef encapc<dims,grid_base_impl<dim,T,S,layout,layout_base>::value_type,layout > encap_dst;
			// destination object type
			typedef encapc<1,prp_object,typename memory_traits_lin<prp_object>::type > encap_src;

			// Copy only the selected properties
			object_s_di<encap_src,encap_dst,OBJ_ENCAP,prp...>(src.get(id),this->get_o(sub_it.get()));

			++id;
			++sub_it;
		}
	}

	/*! \brief unpack the grid object
	 *
	 * \tparam prp properties to unpack
	 *
	 * \param ext preallocated memory from where to unpack the grid
	 * \param obj object where to unpack
	 *
	 */
	template<unsigned int ... prp> void unpack(ExtPreAlloc<S> & mem, Unpack_stat & ps)
	{
		// object that store the information in mem
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,PtrMemory,openfpm::grow_policy_identity> stype;
		stype svect;


		// Calculate the size to pack the object
		size_t size = svect.calculateMem(this->size(),0);

		// Create an object over a pointer (No allocation is produced)
		stype src;
		src.setMemory(mem);
		src.resize(this->size());

		auto it = this->getIterator();

		unpack_with_iterator<decltype(it),stype,prp...>(mem,it,src,ps);

		ps.addOffset(size);
	}

	/*! \brief unpack the sub-grid object
	 *
	 * \tparam prp properties to unpack
	 *
	 * \param mem preallocated memory from where to unpack the object
	 * \param sub sub-grid iterator
	 * \param obj object where to unpack
	 *
	 */
	template<unsigned int ... prp> void unpack(ExtPreAlloc<S> & mem, grid_key_dx_iterator_sub<dims> & sub_it, Unpack_stat & ps)
	{
		// object that store the information in mem
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,PtrMemory, typename memory_traits_lin<prp_object>::type, memory_traits_lin ,openfpm::grow_policy_identity> stype;

		size_t size = stype::template calculateMem(sub_it.getVolume(),0);

		// Create an object over the preallocated memory (No allocation is produced)
		PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));

		// Create an object of the packed information over a pointer (No allocation is produced)
		stype src;
		src.setMemory(ptr);
		src.resize(sub_it.getVolume());

		unpack_with_iterator<grid_key_dx_iterator_sub<dims>,stype,prp...>(mem,sub_it,src,ps);

		ps.addOffset(size);
	}

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

	//! Error code
	size_t err_code;

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
			size_t * err_code_pointer = (size_t *)&this->err_code;
			*err_code_pointer = 1001;
			ACTION_ON_ERROR(GRID_ERROR);
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
				size_t * err_code_pointer = (size_t *)&this->err_code;
				*err_code_pointer = 1002;
				ACTION_ON_ERROR(GRID_ERROR);
			}
			else if (v1.get(i) < 0)
			{
				std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << v1.get(i) << " is negative " << "\n";
				size_t * err_code_pointer = (size_t *)&this->err_code;
				*err_code_pointer = 1003;
				ACTION_ON_ERROR(GRID_ERROR);
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
	inline void check_bound(const grid_base_impl<dim,T,S,layout,layout_base> & g,const grid_key_dx<dim> & key2) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (key2.get(i) >= (long int)g.g1.size(i))
			{
				std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << key2.get(i) << " >= " << g.g1.size(i) << "\n";
				size_t * err_code_pointer = (size_t *)&this->err_code;
				*err_code_pointer = 1004;
				ACTION_ON_ERROR(GRID_ERROR);
			}
			else if (key2.get(i) < 0)
			{
				std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << key2.get(i) << " is negative " << "\n";
				size_t * err_code_pointer = (size_t *)&this->err_code;
				*err_code_pointer = 1005;
				ACTION_ON_ERROR(GRID_ERROR);
			}
		}
	}

#endif

public:

	//! it define that it is a grid
	typedef int yes_i_am_grid;

	//! Definition of the layout
	typedef typename memory_traits_lin<typename T::type>::type memory_lin;

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef encapc<dim,T,layout> container;

	// The object type the grid is storing
	typedef T value_type;

	//! Default constructor
	grid_base_impl() THROW
	:g1(0),isExternal(false),err_code(0)
	{
		// Add this pointer
#ifdef SE_CLASS2
		check_new(this,8,GRID_EVENT,1);
#endif
	}

	/*! \brief create a grid from another grid
	 *
	 * \tparam S memory type for allocation
	 *
	 * \param g the grid to copy
	 * \param mem memory object (only used for template deduction)
	 *
	 */
	grid_base_impl(const grid_base_impl & g) THROW
	:isExternal(false),err_code(0)
	{
		this->operator=(g);
	}

	/*! \brief create a grid of size sz on each direction
	 *
	 * \tparam S memory type for allocation
	 *
	 * \param g the grid to copy
	 * \param mem memory object (only used for template deduction)
	 *
	 */
	grid_base_impl(const size_t & sz) THROW
	:g1(sz),isExternal(false),err_code(0)
	{
		// Add this pointer
#ifdef SE_CLASS2
		check_new(this,8,GRID_EVENT,1);
#endif
	}

	//! Constructor allocate memory and give them a representation
	grid_base_impl(const size_t (& sz)[dim]) THROW
	:g1(sz),isExternal(false),err_code(0)
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

	/*! \brief It copy an operator
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

	/*! \brief It copy an operator
	 *
	 * \param g grid to copy
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
	 * \param g1 grid to check
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
	 * to have contiguous in memory vectors.
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
#ifdef SE_CLASS2
		if (check_valid(&boost::fusion::at_c<p>(data_.mem_r->operator[](g1.LinId(v1))),sizeof(typename type_cpu_prop<p,memory_lin>::type)) == false) {ACTION_ON_ERROR()};
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
#ifdef SE_CLASS2
		if (check_valid(&boost::fusion::at_c<p>(data_.mem_r->operator[](g1.LinId(v1))),sizeof(typename type_cpu_prop<p,memory_lin>::type)) == false) {ACTION_ON_ERROR()};
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
#ifdef SE_CLASS2
		check_valid(&data_.mem_r->operator[](g1.LinId(v1)),sizeof(T));
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
#ifdef SE_CLASS2
		if (check_valid(&data_.mem_r->operator[](g1.LinId(v1)),sizeof(T)) == false)	{ACTION_ON_ERROR()}
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

			//! create a source grid iterator
			grid_key_dx_iterator<dim> it(g1);

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
	 */

	inline grid_key_dx_iterator_sub<dim> getIterator(const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop) const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		// get the starting point and the end point of the real domain

		return grid_key_dx_iterator_sub<dim>(g1,start,stop);
	}

	/*! \brief Return the last error
	 *
	 */
	size_t getLastError()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return err_code;
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
