/*
 * grid_base_impl_layout.hpp
 *
 *  Created on: May 2, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPL_LAYOUT_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPL_LAYOUT_HPP_

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
		if (sz != 0)	t.allocate(sz);
	}
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called
 *
 * \param T Type of memory allocator
 *
 */

template<typename s_m>
struct frswap
{
	s_m & swap_src;
	s_m & swap_dst;

	//! constructor
	frswap(s_m & swap_dst, s_m & swap_src)
	:swap_src(swap_src),swap_dst(swap_dst)
	{};

	//! It call the allocate function for each member
	template<typename T>
	void operator()(T& t) const
	{
		boost::fusion::at_c<T::value>(swap_dst).swap(boost::fusion::at_c<T::value>(swap_src));
	}
};


//! Case memory_traits_lin
template<unsigned int p, typename layout, typename data_type, typename g1_type, typename key_type, unsigned int sel = 2*is_layout_mlin<layout>::value + is_layout_inte<layout>::value >
struct mem_get
{
        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	__host__ __device__ static inline auto get(data_type & data_, const g1_type & g1, const key_type & v1) -> decltype(boost::fusion::at_c<p>(data_.mem_r.operator[](g1.LinId(v1)))) &
	{
		return boost::fusion::at_c<p>(data_.mem_r.operator[](g1.LinId(v1)));
	}

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	__host__ __device__ static inline auto get_lin(data_type & data_, const g1_type & g1, const size_t lin_id) -> decltype(boost::fusion::at_c<p>(data_.mem_r.operator[](lin_id))) &
	{
		return boost::fusion::at_c<p>(data_.mem_r.operator[](lin_id));
	}

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	__host__ __device__ static inline auto get_c(const data_type & data_, const g1_type & g1, const key_type & v1) -> decltype(boost::fusion::at_c<p>(data_.mem_r.operator[](g1.LinId(v1)))) &
	{
		return boost::fusion::at_c<p>(data_.mem_r.operator[](g1.LinId(v1)));
	}

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	__host__ __device__ static inline auto get_lin_c(const data_type & data_, const g1_type & g1, const size_t lin_id) -> decltype(boost::fusion::at_c<p>(data_.mem_r.operator[](lin_id))) &
	{
		return boost::fusion::at_c<p>(data_.mem_r.operator[](lin_id));
	}
};

//! Case memory_traits_inte
template<unsigned int p, typename layout, typename data_type, typename g1_type, typename key_type>
struct mem_get<p,layout,data_type,g1_type,key_type,1>
{
        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	__host__ __device__ static inline auto get(data_type & data_, const g1_type & g1, const key_type & v1) -> decltype(boost::fusion::at_c<p>(data_).mem_r.operator[](g1.LinId(v1)))
	{
		return boost::fusion::at_c<p>(data_).mem_r.operator[](g1.LinId(v1));
	}

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	__host__ __device__ static inline auto get_lin(data_type & data_, const g1_type & g1, size_t lin_id) -> decltype(boost::fusion::at_c<p>(data_).mem_r.operator[](lin_id))
	{
		return boost::fusion::at_c<p>(data_).mem_r.operator[](lin_id);
	}

	 /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a const reference to the object selected
         *
         */
	__host__ __device__ static inline auto get_c(const data_type & data_, const g1_type & g1, const key_type & v1) -> decltype(boost::fusion::at_c<p>(data_).mem_r.operator[](g1.LinId(v1)))
	{
		return boost::fusion::at_c<p>(data_).mem_r.operator[](g1.LinId(v1));
	}

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a const reference to the object selected
         *
         */
	__host__ __device__ static inline auto get_lin_c(const data_type & data_, const g1_type & g1, size_t lin_id) -> decltype(boost::fusion::at_c<p>(data_).mem_r.operator[](lin_id))
	{
		return boost::fusion::at_c<p>(data_).mem_r.operator[](lin_id);
	}
};


//! Case memory_traits_lin
template<typename S, typename layout, typename data_type, typename g1_type, unsigned int sel = 2*is_layout_mlin<layout>::value + is_layout_inte<layout>::value >
struct mem_setm
{
	static inline void setMemory(data_type & data_, const g1_type & g1, bool & is_mem_init)
	{
		S * mem = new S();

		//! Create and set the memory allocator
		data_.setMemory(*mem);

		//! Allocate the memory and create the representation
		if (g1.size() != 0) data_.allocate(g1.size());

		is_mem_init = true;
	}
};

//! Case memory_traits_inte
template<typename S, typename layout, typename data_type, typename g1_type>
struct mem_setm<S,layout,data_type,g1_type,1>
{
	static inline void setMemory(data_type & data_, const g1_type & g1, bool & is_mem_init)
	{
		//! Create an allocate object
		allocate<S> all(g1.size());

		//! for each element in the vector allocate the buffer
		boost::fusion::for_each(data_,all);

		is_mem_init = true;
	}
};


//! Case memory_traits_lin
template<unsigned int dim , typename T, typename layout, typename data_type, typename g1_type, typename key_type, unsigned int sel = 2*is_layout_mlin<layout>::value + is_layout_inte<layout>::value >
struct mem_geto
{
	__device__ __host__ static inline encapc<dim,T,typename layout::type> get(data_type & data_, const g1_type & g1, const key_type & v1)
	{
		return encapc<dim,T,typename layout::type>(data_.mem_r.operator[](g1.LinId(v1)));
	}

	static inline encapc<dim,T,typename layout::type> get_lin(data_type & data_, const size_t & v1)
	{
		return encapc<dim,T,typename layout::type>(data_.mem_r.operator[](v1));
	}
};

//! Case memory_traits_inte
template<unsigned int dim, typename T,typename layout, typename data_type, typename g1_type, typename key_type>
struct mem_geto<dim,T,layout,data_type,g1_type,key_type,1>
{
	__device__ __host__ static inline encapc<dim,T,typename layout::type> get(data_type & data_, const g1_type & g1, const key_type & v1)
	{
		return encapc<dim,T,typename layout::type>(data_,g1.LinId(v1));
	}

	static inline encapc<dim,T,typename layout::type> get_lin(data_type & data_, const size_t & v1)
	{
		return encapc<dim,T,typename layout::type>(data_,v1);
	}
};


//! Case memory_traits_lin
template<typename grid_type, typename S , typename layout, typename data_type, unsigned int sel = 2*is_layout_mlin<layout>::value + is_layout_inte<layout>::value >
struct mem_setext
{
	static inline void set(grid_type & grid_new, grid_type & old, data_type & data_)
	{
		grid_new.setMemory(static_cast<S&>(data_.getMemory()));

		// Create an empty memory allocator for the actual structure

		old.setMemory();
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to set each property with external memory
 *
 *
 */
template<typename grid_type, typename data_type, typename S>
struct mem_setext_prp
{
	grid_type & grid_new;

	grid_type & old;

	data_type & data_;

	/*! \brief constructor
	 *
	 *
	 *
	 */
	inline mem_setext_prp(grid_type & g_new, grid_type & old, data_type & data_)
	:grid_new(g_new),old(old),data_(data_)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		grid_new.template setMemory<T::value>(static_cast<S&>(boost::fusion::at_c<T::value>(data_).getMemory()));


	}
};

//! Case memory_traits_inte
template<typename grid_type, typename S , typename layout, typename data_type>
struct mem_setext<grid_type,S,layout,data_type,1>
{
	static inline void set(grid_type & grid_new, grid_type & old, data_type & data_)
	{
//		std::cerr << __FILE__ << ":" << "__LINE__" << " Error, " << demangle(typeid(grid_type).name()) << " this structure does not support setting from external memory" << std::endl;
		mem_setext_prp<grid_type,data_type,S> mst(grid_new,old,data_);

		// Create an empty memory allocator for the actual structure
		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,grid_type::value_type::max_prop>>(mst);

		old.setMemory();
	}
};


//! Case memory_traits_lin
template<typename T , typename layout, typename data_type, typename grid_type, unsigned int sel = 2*is_layout_mlin<layout>::value + is_layout_inte<layout>::value >
struct mem_swap
{
	static inline void swap(data_type & data_dst, data_type & data_src)
	{
		// move the data
		data_dst.swap(data_src);
	}
};

//! Case memory_traits_inte
template<typename T , typename layout, typename data_type, typename grid_type>
struct mem_swap<T,layout,data_type,grid_type,1>
{
	static inline void swap(data_type & data_dst, data_type & data_src)
	{
		// swap the data for each property
		frswap<decltype(data_dst)> sw(data_dst,data_src);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(sw);
	}
};

template<typename data_type, typename layout, unsigned int sel = 2*is_layout_mlin<layout>::value + is_layout_inte<layout>::value>
struct mem_getpointer
{
	template<unsigned int d> __device__ __host__ static void * getPointer(data_type & data_)
	{
		return data_.mem_r.get_pointer();
	}

	template<unsigned int d> __device__ __host__ static void * getPointer(const data_type & data_)
	{
		return data_.mem_r.get_pointer();
	}
};

template<typename data_type, typename layout>
struct mem_getpointer<data_type,layout,1>
{
	template<unsigned int p> __device__ __host__ static void * getPointer(data_type & data_)
	{
		return boost::fusion::at_c<p>(data_).mem_r.get_pointer();
	}

	template<unsigned int p> __device__ __host__ static const void * getPointer(const data_type & data_)
	{
		return boost::fusion::at_c<p>(data_).mem_r.get_pointer();
	}
};

template<typename data_type, typename Mem_type, typename layout, unsigned int sel = 2*is_layout_mlin<layout>::value + is_layout_inte<layout>::value>
struct mem_setmemory
{
	template<unsigned int d> static void setMemory(data_type & data_, Mem_type & m, size_t sz, bool np)
	{
		//! Create and set the memory allocator
		data_.setMemory(m);

		//! Allocate the memory and create the reppresentation
		if (sz != 0) data_.allocate(sz,np);
	}
};

template<typename data_type, typename Mem_type, typename layout>
struct mem_setmemory<data_type,Mem_type,layout,1>
{
	template<unsigned int p> static void setMemory(data_type & data_, Mem_type & m, size_t sz,bool np)
	{
		//! Create and set the memory allocator
		boost::fusion::at_c<p>(data_).setMemory(m);

		//! Allocate the memory and create the reppresentation
		if (sz != 0) boost::fusion::at_c<p>(data_).allocate(sz,np);
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPL_LAYOUT_HPP_ */
