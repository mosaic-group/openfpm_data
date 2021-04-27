/*
 * grid_base_implementation.hpp
 *
 *  Created on: May 2, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPLEMENTATION_HPP_
#define OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPLEMENTATION_HPP_

#include "grid_base_impl_layout.hpp"
#include "util/cuda_util.hpp"
#include "cuda/cuda_grid_gpu_funcs.cuh"
#include "util/create_vmpl_sequence.hpp"
#include "util/cuda_launch.hpp"
#include "util/object_si_di.hpp"

constexpr int DATA_ON_HOST = 32;
constexpr int DATA_ON_DEVICE = 64;
constexpr int EXACT_RESIZE = 128;

template<bool np,typename T>
struct skip_init
{
	static bool skip_()
	{
		return true;
	}
};

template<typename T>
struct skip_init<true,T>
{
	static bool skip_()
	{
		return T::noPointers();
	}
};

#ifdef CUDA_GPU

#define GRID_ID_3_RAW(start,stop) int x[3] = {threadIdx.x + blockIdx.x * blockDim.x + start.get(0),\
    									 threadIdx.y + blockIdx.y * blockDim.y + start.get(1),\
										 threadIdx.z + blockIdx.z * blockDim.z + start.get(2)};\
										 \
										 if (x[0] > stop.get(0) || x[1] > stop.get(1) || x[2] > stop.get(2))\
    									 {return;}

#define GRID_ID_3_TRAW(start,stop) int tx = threadIdx.x + blockIdx.x * blockDim.x + start.get(0);\
    							   int ty = threadIdx.y + blockIdx.y * blockDim.y + start.get(1);\
								   int tz = threadIdx.z + blockIdx.z * blockDim.z + start.get(2);\
										 \
										 if (tx > stop.get(0) || ty > stop.get(1) || tz > stop.get(2))\
    									 {return;}

#define GRID_ID_3(ite_gpu) grid_key_dx<3,int> key;\
							  key.set_d(0,threadIdx.x + blockIdx.x * blockDim.x + ite_gpu.start.get(0));\
    						  key.set_d(1,threadIdx.y + blockIdx.y * blockDim.y + ite_gpu.start.get(1));\
							  key.set_d(2,threadIdx.z + blockIdx.z * blockDim.z + ite_gpu.start.get(2));\
										 \
										 if (key.get(0) > ite_gpu.stop.get(0) || key.get(1) > ite_gpu.stop.get(1) || key.get(2) > ite_gpu.stop.get(2))\
    									 {return;}

#define GRID_ID_2(ite_gpu) grid_key_dx<2,int> key;\
							  key.set_d(0,threadIdx.x + blockIdx.x * blockDim.x + ite_gpu.start.get(0));\
    						  key.set_d(1,threadIdx.y + blockIdx.y * blockDim.y + ite_gpu.start.get(1));\
										 \
										 if (key.get(0) > ite_gpu.stop.get(0) || key.get(1) > ite_gpu.stop.get(1))\
    									 {return;}

#ifdef __NVCC__

template<unsigned int dim, typename ids_type = int>
struct grid_p
{
	__device__ static inline grid_key_dx<dim,ids_type> get_grid_point(const grid_sm<dim,void> & g)
	{
		grid_key_dx<dim,ids_type> key;

		key.set_d(0,blockIdx.x * blockDim.x + threadIdx.x);
		key.set_d(1,blockIdx.y * blockDim.y + threadIdx.y);

		unsigned int bz = blockIdx.z * blockDim.z + threadIdx.z;
		key.set_d(2,bz % g.size(2));

		for (unsigned int i = 3 ; i < dim ; i++)
		{
			bz /= g.size(i);
			key.set_d(i,bz % g.size(i));
		}

		return key;
	}

	__device__ static inline grid_key_dx<dim,ids_type> get_grid_point(const openfpm::array<ids_type,dim,unsigned int> & g)
	{
		grid_key_dx<dim,ids_type> key;

		key.set_d(0,blockIdx.x * blockDim.x + threadIdx.x);
		key.set_d(1,blockIdx.y * blockDim.y + threadIdx.y);

		unsigned int bz = blockIdx.z * blockDim.z + threadIdx.z;
		key.set_d(2,bz % g[2]);

		for (unsigned int i = 3 ; i < dim ; i++)
		{
			bz /= g[i];
			key.set_d(i,bz % g[i]);
		}

		return key;
	}
};

template<typename ids_type>
struct grid_p<3,ids_type>
{
	__device__ static inline grid_key_dx<3,ids_type> get_grid_point(const grid_sm<3,void> & g)
	{
		grid_key_dx<3,unsigned int> key;

		key.set_d(0,blockIdx.x * blockDim.x + threadIdx.x);
		key.set_d(1,blockIdx.y * blockDim.y + threadIdx.y);
		key.set_d(2,blockIdx.z * blockDim.z + threadIdx.z);

		return key;
	}

	__device__ static inline grid_key_dx<3,ids_type> get_grid_point(const openfpm::array<ids_type,3,unsigned int> & g)
	{
		grid_key_dx<3,ids_type> key;

		key.set_d(0,blockIdx.x * blockDim.x + threadIdx.x);
		key.set_d(1,blockIdx.y * blockDim.y + threadIdx.y);
		key.set_d(2,blockIdx.z * blockDim.z + threadIdx.z);

		return key;
	}
};

template<typename ids_type>
struct grid_p<2,ids_type>
{
	__device__ static inline grid_key_dx<2,ids_type> get_grid_point(const grid_sm<2,void> & g)
	{
		grid_key_dx<2,ids_type> key;

		key.set_d(0,blockIdx.x * blockDim.x + threadIdx.x);
		key.set_d(1,blockIdx.y * blockDim.y + threadIdx.y);

		return key;
	}

	__device__ static inline grid_key_dx<2,ids_type> get_grid_point(const openfpm::array<ids_type,2,unsigned int> & g)
	{
		grid_key_dx<2,ids_type> key;

		key.set_d(0,blockIdx.x * blockDim.x + threadIdx.x);
		key.set_d(1,blockIdx.y * blockDim.y + threadIdx.y);

		return key;
	}
};

template<typename ids_type>
struct grid_p<1,ids_type>
{
	__device__ static inline grid_key_dx<1,unsigned int> get_grid_point(const grid_sm<1,void> & g)
	{
		grid_key_dx<1,unsigned int> key;

		key.set_d(0,blockIdx.x * blockDim.x + threadIdx.x);

		return key;
	}
};

#endif


template<unsigned int dim>
void move_work_to_blocks(ite_gpu<dim> & ite)
{
	if (dim == 1)
	{
		ite.wthr.x = ite.wthr.x * ite.thr.x;
		ite.thr.x = 1;
	}
	else if(dim == 2)
	{
		ite.wthr.x = ite.wthr.x * ite.thr.x;
		ite.wthr.y = ite.wthr.y * ite.thr.y;
		ite.thr.x = 1;
		ite.thr.y = 1;

	}
	else
	{
		ite.wthr.x = ite.wthr.x * ite.thr.x;
		ite.wthr.x = ite.wthr.y * ite.thr.y;
		ite.wthr.x = ite.wthr.z * ite.thr.z;
		ite.thr.x = 1;
		ite.thr.y = 1;
		ite.thr.z = 1;
	}
}

#endif

#include "copy_grid_fast.hpp"

template<typename T>
struct copy_grid_fast_caller
{
	template<typename grid_type>
	static void call(grid_type & gd, const grid_type & gs, const Box<grid_type::dims,size_t> & box_src, const Box<grid_type::dims,size_t> & box_dst)
	{
		std::cout << "Error: " << __FILE__ << ":" << __LINE__ << " copy_grid_fast_caller failure" << std::endl;
	}
};

template<int ... prp>
struct copy_grid_fast_caller<index_tuple_sq<prp ...>>
{
	template<typename grid_type>
	static void call(grid_type & gd, const grid_type & gs, const Box<grid_type::dims,size_t> & box_src, const Box<grid_type::dims,size_t> & box_dst)
	{
        grid_key_dx<grid_type::dims> cnt[1];
        cnt[0].zero();

        typedef typename std::remove_reference<decltype(gd)>::type grid_cp;
        typedef typename std::remove_reference<decltype(gd.getGrid())>::type grid_info_cp;

        typedef typename to_int_sequence<0,grid_type::value_type::max_prop>::type result;

        copy_grid_fast<!is_contiguos<prp...>::type::value || has_pack_gen<typename grid_type::value_type>::value,
                                   grid_type::dims,
                                   grid_cp,
                                   grid_info_cp>::copy(gs.getGrid(),
                                   gd.getGrid(),
                                   box_src,
                                   box_dst,
                                   gs,gd,cnt);
	}
};

/*! \brief
 *
 * Implementation of a N-dimensional grid
 *
 * \tparam dim dimansionality of the grid
 * \tparam T type store by the grid
 * \tparam S Memory pool from where to take the memory
 * \tparam layout_base layout memory meta-function (the meta-function used to construct layout_)
 *
 */
template<unsigned int dim,
         typename T,
         typename S,
         template<typename> class layout_base,
         typename ord_type = grid_sm<dim,void> >
class grid_base_impl
{
	//! memory layout
	typedef typename layout_base<T>::type layout;

public:

	//! memory layout
	typedef layout layout_type;

	//! expose the dimansionality as a static const
	static constexpr unsigned int dims = dim;

	//! Access key
	typedef grid_key_dx<dim> access_key;

	//! boost::vector that describe the data type
	typedef typename T::type T_type;

	//! base layout type
	typedef layout_base<T> layout_base_;

	typedef ord_type linearizer_type;

	typedef T background_type;

protected:

	//! Memory layout specification + memory chunk pointer
	layout data_;

	//! This is a structure that store all information related to the grid and how indexes are linearized
	ord_type g1;

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
#ifndef __CUDA_ARCH__
		if (is_mem_init == false)
		{
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " you must call SetMemory before access the grid\n";
			ACTION_ON_ERROR(GRID_ERROR_OBJECT);
		}
#endif
	}

	/*! \brief Check that the key is inside the grid
	 *
	 * \param key
	 *
	 */
	inline void check_bound(const grid_key_dx<dim> & v1) const
	{
#ifndef __CUDA_ARCH__
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
#endif
	}

	/*! \brief Check that the key is inside the grid
	 *
	 * \param key
	 *
	 */
	inline void check_bound(size_t v1) const
	{
#ifndef __CUDA_ARCH__
		if (v1 >= getGrid().size())
		{
			std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << v1<< " >= " << getGrid().size() << "\n";
			ACTION_ON_ERROR(GRID_ERROR_OBJECT);
		}
#endif
	}

	/*! \brief Check that the key is inside the grid
	 *
	 * check if key2 is inside the g grid boundary
	 *
	 * \param g grid
	 * \param key2
	 *
	 */
	template<typename Mem> inline void check_bound(const grid_base_impl<dim,T,Mem,layout_base, ord_type> & g,const grid_key_dx<dim> & key2) const
	{
#ifndef __CUDA_ARCH__
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
#endif
	}

	/*! \brief Check that the key is inside the grid
	 *
	 * check if key2 is inside the g grid boundary
	 *
	 * \param g grid
	 * \param key2
	 *
	 */
	template<typename Mem, template <typename> class layout_base2> 
	inline void check_bound(const grid_base_impl<dim,T,Mem,layout_base2,ord_type> & g,const grid_key_dx<dim> & key2) const
	{
#ifndef __CUDA_ARCH__
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
#endif
	}

	/*! \brief Check that the key is inside the grid
	 *
	 * check if key2 is inside the g grid boundary
	 *
	 * \param g grid
	 * \param key2
	 *
	 */
	template<typename Mem> inline void check_bound(const grid_base_impl<dim,T,Mem,layout_base> & g,const size_t & key2) const
	{
#ifndef __CUDA_ARCH__
		if (key2 >= g.getGrid().size())
		{
			std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << key2 << " >= " << getGrid().size() << "\n";
			ACTION_ON_ERROR(GRID_ERROR_OBJECT);
		}
#endif
	}

#endif

	void resize_impl_device(const size_t (& sz)[dim],grid_base_impl<dim,T,S,layout_base,ord_type> & grid_new, unsigned int blockSize = 1)
	{
#if defined(CUDA_GPU) && defined(__NVCC__)

			// Compile time-cheking that make sense to call a GPU kernel to copy.
			

			grid_key_dx<dim> start;
			grid_key_dx<dim> stop;

			for (size_t i = 0 ; i < dim ; i++)
			{
				start.set_d(i,0);

				// take the smallest
				if (grid_new.g1.size(i) < sz[i])
				{stop.set_d(i,grid_new.g1.size(i)-1);}
				else
				{stop.set_d(i,sz[i]-1);}
			}

//			if (dim == 1)
//			{
//				copy_fast_1d_device_memory<is_layout_mlin<layout_base<T>>::value,decltype(grid_new.data_),S> cf1dm(data_,grid_new.data_);

//				boost::mpl::for_each_ref<boost::mpl::range_c<int,0,T::max_prop>>(cf1dm);
//			}
			if (dim <= 3)
			{
				auto ite = this->getGPUIterator(start,stop);
				bool has_work = has_work_gpu(ite);

				if (has_work == true)
				{
					if (blockSize == 1)
					{
						CUDA_LAUNCH((copy_ndim_grid_device<dim,decltype(grid_new.toKernel())>),ite,this->toKernel(),grid_new.toKernel());
					}
					else
					{
						move_work_to_blocks(ite);

						ite.thr.x = blockSize;

						CUDA_LAUNCH((copy_ndim_grid_block_device<dim,decltype(grid_new.toKernel())>),ite,this->toKernel(),grid_new.toKernel());
					}
				}
			}
			else
			{
				grid_key_dx<1> start;
				start.set_d(0,0);
				grid_key_dx<1> stop({});
				stop.set_d(0,this->g1.size());

				size_t sz[1];
				sz[0]= this->g1.size();

				grid_sm<1,void> g_sm_copy(sz);

				auto ite = getGPUIterator_impl<1>(g_sm_copy,start,stop);

				CUDA_LAUNCH((copy_ndim_grid_device<dim,decltype(grid_new.toKernel())>),ite,this->toKernel(),grid_new.toKernel());
			}
#else

			std::cout << __FILE__ << ":" << __LINE__ << " error: the function resize require the launch of a kernel, but it seem that this" <<
					                                    " file (grid_base_implementation.hpp) has not been compiled with NVCC  " << std::endl;

#endif
	}

	void resize_impl_host(const size_t (& sz)[dim], grid_base_impl<dim,T,S,layout_base,ord_type> & grid_new)
	{
		size_t sz_c[dim];
		for (size_t i = 0 ; i < dim ; i++)
		{sz_c[i] = (g1.size(i) < sz[i])?g1.size(i):sz[i];}

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
	}

	void resize_impl_memset(grid_base_impl<dim,T,S,layout_base,ord_type> & grid_new)
	{
		//! Set the allocator and allocate the memory
		if (isExternal == true)
		{
			mem_setext<typename std::remove_reference<decltype(grid_new)>::type,S,layout_base<T>,decltype(data_)>::set(grid_new,*this,this->data_);
		}
		else
			grid_new.setMemory();
	}

public:

	// Implementation of packer and unpacker for grid
	#include "grid_pack_unpack.ipp"

	//! it define that this data-structure is a grid
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
	}

	/*! \brief Constructor
	 *
	 * It construct a grid of specified size
	 *
	 * \param sz array that indicate the size of the grid in each dimension
	 *
	 */
	grid_base_impl(const size_t (& sz)[dim]) THROW
	:g1(sz),isExternal(false)
	{
		// Add this pointer
	}

	//! Destructor
	~grid_base_impl() THROW
	{
		// delete this pointer
	}

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 * \return itself
	 *
	 */
	grid_base_impl<dim,T,S,layout_base> & operator=(const grid_base_impl<dim,T,S,layout_base> & g)
	{
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
	grid_base_impl<dim,T,S,layout_base> & operator=(grid_base_impl<dim,T,S,layout_base> && g)
	{
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
	bool operator==(const grid_base_impl<dim,T,S,layout_base> & g)
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
	grid_base_impl<dim,T,S,layout_base> duplicate() const THROW
	{
		//! Create a completely new grid with sz

		grid_base_impl<dim,T,S,layout_base> grid_new(g1.getSize());

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

#ifdef CUDA_GPU

	/*! \brief Get an iterator for the GPU
	 *
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	struct ite_gpu<dim> getGPUIterator(grid_key_dx<dim,long int> & key1, grid_key_dx<dim,long int> & key2, size_t n_thr = default_kernel_wg_threads_) const
	{
		return getGPUIterator_impl<dim>(g1,key1,key2,n_thr);
	}
#endif

	/*! \brief Return the internal grid information
	 *
	 * Return the internal grid information
	 *
	 * \return the internal grid
	 *
	 */

	const ord_type & getGrid() const
	{
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
		mem_setm<S,layout_base<T>,decltype(this->data_),decltype(this->g1)>::setMemory(data_,g1,is_mem_init);
	}

	/*! \brief Set the object that provide memory from outside
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
	template<unsigned int p = 0> void setMemory(S & m)
	{
		//! Is external
		isExternal = true;

		//! Create and set the memory allocator
//		data_.setMemory(m);

		//! Allocate the memory and create the reppresentation
//		if (g1.size() != 0) data_.allocate(g1.size());

		bool skip_ini = skip_init<has_noPointers<T>::value,T>::skip_();

		mem_setmemory<decltype(data_),S,layout_base<T>>::template setMemory<p>(data_,m,g1.size(),skip_ini);

		is_mem_init = true;
	}

	/*! \brief Set the object that provide memory from outside
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
	void setMemoryArray(S * m)
	{
		//! Is external
		isExternal = true;

		bool skip_ini = skip_init<has_noPointers<T>::value,T>::skip_();

		mem_setmemory<decltype(data_),S,layout_base<T>>::template setMemoryArray(*this,m,g1.size(),skip_ini);

		is_mem_init = true;
	}

	/*! \brief Return a plain pointer to the internal data
	 *
	 * Return a plain pointer to the internal data
	 *
	 * \return plain data pointer
	 *
	 */

	template<unsigned int p = 0> void * getPointer()
	{
		return mem_getpointer<decltype(data_),layout_base_>::template getPointer<p>(data_);
	}

	/*! \brief Return a plain pointer to the internal data
	 *
	 * Return a plain pointer to the internal data
	 *
	 * \return plain data pointer
	 *
	 */

	template<unsigned int p = 0> const void * getPointer() const
	{
		return mem_getpointer<decltype(data_),layout_base_>::template getPointer<p>(data_);
	}


	/*! \brief In this case insert is equivalent to get
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(layout_base<T>::template get<p>(data_,g1,grid_key_dx<dim>()))>
	inline r_type insert(const grid_key_dx<dim> & v1)
	{
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return this->get<p>(v1);
	}


	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(layout_base<T>::template get<p>(data_,g1,grid_key_dx<dim>()))>
	__device__ __host__ inline r_type get_usafe(const grid_key_dx<dim> & v1)
	{
#ifdef SE_CLASS1
		check_init();
#endif
		return layout_base<T>::template get<p>(data_,g1,v1);
	}

	/*! \brief Get the const reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the const reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(layout_base<T>::template get_c<p>(data_,g1,grid_key_dx<dim>()))>
	__device__ __host__ inline r_type get_unsafe(const grid_key_dx<dim> & v1) const
	{
#ifdef SE_CLASS1
		check_init();
#endif
		return layout_base<T>::template get_c<p>(data_,g1,v1);
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(layout_base<T>::template get<p>(data_,g1,grid_key_dx<dim>()))>
	__device__ __host__ inline r_type get(const grid_key_dx<dim> & v1)
	{
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return layout_base<T>::template get<p>(data_,g1,v1);
	}

	/*! \brief Get the point flag (in this case just return 0)
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return zero
	 *
	 */
	__device__ __host__ inline unsigned char getFlag(const grid_key_dx<dim> & v1) const
	{
		return 0;
	}

	/*! \brief Get the const reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the const reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(layout_base<T>::template get_c<p>(data_,g1,grid_key_dx<dim>()))>
	__device__ __host__ inline r_type get(const grid_key_dx<dim> & v1) const
	{
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return layout_base<T>::template get_c<p>(data_,g1,v1);
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param lin_id linearized element that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(layout_base<T>::template get_lin<p>(data_,g1,0))>
	__device__ __host__ inline r_type get(const size_t lin_id)
	{
#ifdef SE_CLASS1
		check_init();
		check_bound(lin_id);
#endif
		return layout_base<T>::template get_lin<p>(data_,g1,lin_id);
	}

	/*! \brief Get the const reference of the selected element
	 *
	 * \param lin_id linearized element that identify the element in the grid
	 *
	 * \return the const reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(layout_base<T>::template get_lin<p>(data_,g1,0))>
	__device__ __host__ inline const r_type get(size_t lin_id) const
	{
#ifdef SE_CLASS1
		check_init();
		check_bound(lin_id);
#endif
		return layout_base<T>::template get_lin_const(data_,g1,lin_id);
	}


	/*! \brief Get the of the selected element as a boost::fusion::vector
	 *
	 * Get the selected element as a boost::fusion::vector
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \see encap_c
	 *
	 * \return an encap_c that is the representation of the object (careful is not the object)
	 *
	 */
	inline encapc<dim,T,layout> get_o(const grid_key_dx<dim> & v1)
	{
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
	 * \see encap_c
	 *
	 * \return an encap_c that is the representation of the object (careful is not the object)
	 *
	 */
	inline const encapc<dim,T,layout> get_o(const grid_key_dx<dim> & v1) const
	{
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return mem_geto<dim,T,layout_base<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>
			   ::get(const_cast<typename std::add_lvalue_reference<decltype(this->data_)>::type>(data_),
					 g1,v1);
	}


	/*! \brief Get the of the selected element as a boost::fusion::vector
	 *
	 * Get the selected element as a boost::fusion::vector
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \see encap_c
	 *
	 * \return an encap_c that is the representation of the object (careful is not the object)
	 *
	 */
	inline encapc<dim,T,layout> insert_o(const grid_key_dx<dim> & v1)
	{
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
	 * \param v1 linearized id that identify the element in the grid
	 *
	 * \see encap_c
	 *
	 * \return an encap_c that is the representation of the object (careful is not the object)
	 *
	 */
	inline encapc<dim,T,layout> get_o(size_t v1)
	{
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return mem_geto<dim,T,layout_base<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get_lin(data_,v1);
	}

	/*! \brief Get the of the selected element as a boost::fusion::vector
	 *
	 * Get the selected element as a boost::fusion::vector
	 *
	 * \param v1 linearized id that identify the element in the grid
	 *
	 * \see encap_c
	 *
	 * \return an encap_c that is the representation of the object (careful is not the object)
	 *
	 */
	inline const encapc<dim,T,layout> get_o(size_t v1) const
	{
#ifdef SE_CLASS1
		check_init();
		check_bound(v1);
#endif
		return mem_geto<dim,T,layout_base<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>
		       ::get_lin(const_cast<typename std::add_lvalue_reference<decltype(this->data_)>::type>(data_),v1);
	}

	/*! \brief Fill the memory with the selected byte
	 *
	 * \warning It is a low level memory operation it ignore any type and semantic safety
	 *
	 * \param fl byte pattern to fill
	 *
	 */
	template<int prp>
	void fill(unsigned char fl)
	{
		if (prp != 0 || is_layout_mlin<layout_base<T>>::type::value == false)
		{
			std::cout << "Error: " << __FILE__ << ":" << __LINE__ << " unsupported fill operation " << std::endl;
		}

		memset(getPointer(),fl,size() * sizeof(T));
	}

	/*! \brief Remove all the points in this region
	 *
	 * For this case this function does nothing
	 *
	 * \param box_src box to kill the points
	 *
	 */
	void remove(Box<dim,long int> & section_to_delete)
	{}

	/*! \brief Reset the queue to remove and copy section of grids
	 *
	 * \note for this particular implementation it does nothing
	 *
	 */
	void copyRemoveReset()
	{
	}

    /*! \brief This function check if keep geometry is possible for this grid
     *
     * \return false it does not exist this feature on this type of grids
     *
     */
    bool isSkipLabellingPossible()
    {
    	return false;
    }

	/*! \brief copy an external grid into a specific place into this grid
	 *
	 * It copy the area indicated by the  box_src from grid_src into this grid
	 * at the place box_dst. The volume of box_src and box_dst
	 *
	 * box_dst and box_src are adjusted if box_dst overflow the destination grid
	 *
	 * \param grid_src source grid
	 * \param box_src source box
	 * \param box_dst destination box
	 *
	 */
	void copy_to(const grid_base_impl<dim,T,S,layout_base> & grid_src,
			     const Box<dim,long int> & box_src,
				 const Box<dim,long int> & box_dst)
	{
		// fix box_dst

	     Box<dim,size_t> box_src_;
		 Box<dim,size_t> box_dst_;

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (box_dst.getHigh(i) >= (long int)g1.size(i))
			{
				long int shift = box_dst.getHigh(i) - g1.size(i) + 1;
				box_dst_.setHigh(i,box_dst.getHigh(i) - shift);
				box_src_.setHigh(i,box_src.getHigh(i) - shift);
			}
			else
			{
				box_dst_.setHigh(i,box_dst.getHigh(i));
				box_src_.setHigh(i,box_src.getHigh(i));
			}

			if (box_dst.getLow(i) < 0)
			{
				long int shift = -box_dst.getLow(i);
				box_dst_.setLow(i,box_dst.getLow(i) - shift);
				box_src_.setLow(i,box_src.getLow(i) - shift);
			}
			else
			{
				box_dst_.setLow(i,box_dst.getLow(i));
				box_src_.setLow(i,box_src.getLow(i));
			}
		}

        typedef typename to_int_sequence<0,T::max_prop>::type result;

        copy_grid_fast_caller<result>::call(*this,grid_src,box_src_,box_dst_);

/*        copy_grid_fast<!is_contiguos<prp...>::type::value || has_pack_gen<typename device_grid::value_type>::value,
                                   dim,
                                   grid_cp,
                                   grid_info_cp>::copy(grid_src.getGrid(),
                                   this->getGrid(),
                                   box_src,
                                   box_dst,
                                   grid_src,*this,cnt);*/

		/////////////////////////////////////////
	}

	/*! \brief copy an external grid into a specific place into this grid
	 *
	 * It copy the area indicated by the  box_src from grid_src into this grid
	 * at the place box_dst. The volume of box_src and box_dst
	 *
	 * \param grid_src source grid
	 * \param box_src source box
	 * \param box_dst destination box
	 *
	 */
	template<unsigned int ... prp>
	void copy_to_prp(const grid_base_impl<dim,T,S,layout_base> & grid_src,
			     const Box<dim,size_t> & box_src,
				 const Box<dim,size_t> & box_dst)
	{
        typedef typename std::remove_reference<decltype(grid_src)>::type grid_cp;
        typedef typename std::remove_reference<decltype(grid_src.getGrid())>::type grid_info_cp;

        grid_key_dx<dim> cnt[1];
        cnt[0].zero();

        copy_grid_fast<!is_contiguos<prp...>::type::value || has_pack_gen<T>::value,
                                   dim,
                                   grid_cp,
                                   grid_info_cp>::copy(grid_src.getGrid(),
                                   this->getGrid(),
                                   box_src,
                                   box_dst,
                                   grid_src,*this,cnt);
	}

	/*! \brief It does nothing
	 *
	 *
	 *
	 */
	void clear()
	{}

	/*! \brief copy an external grid into a specific place into this grid
	 *
	 * It copy the area indicated by the  box_src from grid_src into this grid
	 * at the place box_dst. The volume of box_src and box_dst
	 *
	 * \param grid_src source grid
	 * \param box_src source box
	 * \param box_dst destination box
	 *
	 */
	template<template<typename,typename> class op, unsigned int ... prp>
	void copy_to_op(const grid_base_impl<dim,T,S,layout_base> & gs,
			     const Box<dim,size_t> & bx_src,
				 const Box<dim,size_t> & bx_dst)
	{
		grid_key_dx_iterator_sub<dim> sub_src(gs.getGrid(),bx_src.getKP1(),bx_src.getKP2());
		grid_key_dx_iterator_sub<dim> sub_dst(this->getGrid(),bx_dst.getKP1(),bx_dst.getKP2());

		while (sub_src.isNext())
		{
			// write the object in the last element
			object_si_di_op<op,decltype(gs.get_o(sub_src.get())),decltype(this->get_o(sub_dst.get())),OBJ_ENCAP,prp...>(gs.get_o(sub_src.get()),this->get_o(sub_dst.get()));

			++sub_src;
			++sub_dst;
		}
	}

	/*! \brief Indicate that unpacking the header is supported
     *
	 * \return false
	 * 
	 */
	static bool is_unpack_header_supported()
	{return false;}

	/*! \brief Resize the grid
	 *
	 * Resize the grid to the old information is retained on the new grid,
	 * if the new grid is bigger. if is smaller the data are cropped
	 *
	 * \param sz reference to an array of dimension dim
	 * \param opt options for resize. In case we know that the data are only on device memory we can use DATA_ONLY_DEVICE,
	 *                                In case we know that the data are only on host memory we can use DATA_ONLY_HOST
	 *
	 * \param blockSize The default is equal to 1. In case of accelerator buffer resize indicate the size of the block of
	 *        threads ( this is used in case of a vector of blocks where the block object override to operator=
	 *                  to distribute threads on each block element )
	 *
	 */
	void resize(const size_t (& sz)[dim], size_t opt = DATA_ON_HOST | DATA_ON_DEVICE, unsigned int blockSize = 1)
	{
		//! Create a completely new grid with sz

		grid_base_impl<dim,T,S,layout_base,ord_type> grid_new(sz);

		resize_impl_memset(grid_new);


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

		if (opt & DATA_ON_HOST)
		{resize_impl_host(sz,grid_new);}

		if (opt & DATA_ON_DEVICE && S::isDeviceHostSame() == false)
		{resize_impl_device(sz,grid_new,blockSize);}

//		}

		// copy grid_new to the base

		this->swap(grid_new);
	}

	/*! \brief Resize the space
	 *
	 * Resize the space to a new grid, the element are retained on the new grid,
	 * if the new grid is bigger the new element are now initialized, if is smaller
	 * the data are cropped
	 *
	 * \param sz reference to an array of dimension dim
	 * \param opt options for resize. In case we know that the data are only on device memory we can use DATA_ONLY_DEVICE,
	 *                                In case we know that the data are only on host memory we can use DATA_ONLY_HOST
	 *
	 */
	void resize_no_device(const size_t (& sz)[dim])
	{
		//! Create a completely new grid with sz

		grid_base_impl<dim,T,S,layout_base> grid_new(sz);

		resize_impl_memset(grid_new);
		resize_impl_host(sz,grid_new);

		this->swap(grid_new);
	}

	/*! \brief Remove one element valid only on 1D
	 *
	 * \param key element to remove
	 *
	 */
	void remove(size_t key)
	{
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

	/*! \brief It swap the objects A become B and B become A using A.swap(B);
	 *
	 * This is a different from the standard swap and require long explanation.
	 *
	 * This object by default when it construct after we call setMemory() it create an internal memory object
	 * and use it allocate memory internally.
	 *
	 * If instead we use setMemory(external_mem) this object does not create an internal memory object but use
	 * the passed object to allocate memory. Because the external memory can already have a pool of memory preallocated
	 * we can re-use the memory.
	 *
	 * Despite this setMemory can be used to do memory retaining/re-use  and/or garbage collection.
	 * It can be seen from a different prospective of making the data structures act like a representation of external
	 *  memory. De facto we are giving meaning to the external memory so we are shaping or re-shaping pre-existing
	 *   external memory.
	 *
	 * In the following I will call these two mode Mode1 and Mode2
	 *
	 * Using the structure in this way has consequences, because now in Mode2 the memory (and so its life-span) is disentangled
	 *  by its structure.
	 *
	 *
	 * The main difference comes when we swap object in which one of both are in Mode2
	 *
	 * Let's suppose object A is in Mode1 and object B is is Mode2. The normal swap, fully swap the objects
	 *
	 * A.swap(B) A become B (in mode 2) and B become A (in mode 1)
	 *
	 * A.swap_nomode(B) In this case the mode is not swapped  A become B (in mode 1) and B become A (in mode 2).
	 *                  So the mode is not swapped and remain the original
	 *
	 * \param grid object B
	 *
	 */

	void swap_nomode(grid_base_impl<dim,T,S,layout_base> & grid)
	{
		mem_swap<T,layout_base<T>,decltype(data_),decltype(grid)>::template swap_nomode<S>(data_,grid.data_);

		// exchange the grid info
		g1.swap(grid.g1);

		// exchange the init status
		bool exg = is_mem_init;
		is_mem_init = grid.is_mem_init;
		grid.is_mem_init = exg;
	}

	/*! \brief It swap the objects A become B and B become A using A.swap(B);
	 *
	 * \param grid to swap with
	 *
	 */

	void swap(grid_base_impl<dim,T,S,layout_base,ord_type> & grid)
	{
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

	void swap(grid_base_impl<dim,T,S,layout_base> && grid)
	{
		swap(grid);
	}

	/*! \brief set only some properties
	 *
	 * \param key1 destination point
	 * \param g source
	 * \param key2 source point
	 */
	template<unsigned int ... prp>
	__device__ __host__ inline void set(const grid_key_dx<dim> & key1,const grid_base_impl & g, const grid_key_dx<dim> & key2)
	{
		auto edest = this->get_o(key1);
		auto esrc = g.get_o(key2);

		copy_cpu_encap_encap_prp<decltype(g.get_o(key2)),decltype(this->get_o(key1)),prp...> ec(esrc,edest);

		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,sizeof...(prp)>>(ec);
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
#ifdef SE_CLASS1
		check_init();
		check_bound(dx);
#endif

		// create the object to copy the properties
		copy_cpu_encap<dim,grid_base_impl<dim,T,S,layout_base>,layout> cp(dx,*this,obj);

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

	inline void set(const grid_key_dx<dim> & key1,
			        const grid_base_impl<dim,T,S,layout_base> & g,
					const grid_key_dx<dim> & key2)
	{
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

	inline void set(const size_t key1,
			        const grid_base_impl<dim,T,S,layout_base> & g,
					const size_t key2)
	{
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

	template<typename Mem> inline void set(const grid_key_dx<dim> & key1,const grid_base_impl<dim,T,Mem,layout_base> & g, const grid_key_dx<dim> & key2)
	{
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

	template<typename Mem, template <typename> class layout_base2> inline void set_general(const grid_key_dx<dim> & key1,
												   const grid_base_impl<dim,T,Mem,layout_base2> & g,
												   const grid_key_dx<dim> & key2)
	{
#ifdef SE_CLASS1
		check_init();
		check_bound(key1);
		check_bound(g,key2);
#endif

		this->get_o(key1) = g.get_o(key2);
	}

	/*! \brief return the size of the grid
	 *
	 * \return Return the size of the grid
	 *
	 */
	inline size_t size() const
	{
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
	inline grid_key_dx_iterator_sub<dim> getSubIterator(const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop) const
	{
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
		size_t sz[dim];

		for (int i = 0 ; i < dim ; i++)
		{sz[i] = g1.size(i);}

		grid_sm<dim,void> gvoid(sz);

		return grid_key_dx_iterator<dim>(gvoid);
	}

#ifdef CUDA_GPU

	/*! \brief Convert the grid into a data-structure compatible for computing into GPU
	 *
	 *  The object created can be considered like a reference of the original
	 *
	 */
	grid_gpu_ker<dim,T,layout_base> toKernel()
	{
		return grid_toKernelImpl<is_layout_inte<layout_base<T>>::value,dim,T>::toKernel(*this);
	}

	/*! \brief Convert the grid into a data-structure compatible for computing into GPU
	 *
	 *  The object created can be considered like a reference of the original
	 *
	 */
	const grid_gpu_ker<dim,T,layout_base> toKernel() const
	{
		return grid_toKernelImpl<is_layout_inte<layout_base<T>>::value,dim,T>::toKernel(*this);
	}

#endif

	/*! \brief Return a grid iterator
	 *
	 * Return a grid iterator, to iterate through the grid with stencil calculation
	 *
	 * \return a grid iterator with stencil calculation
	 *
	 */
	template<unsigned int Np>
	inline grid_key_dx_iterator<dim,stencil_offset_compute<dim,Np>>
	getIteratorStencil(const grid_key_dx<dim> (& stencil_pnt)[Np]) const
	{
		return grid_key_dx_iterator<dim,stencil_offset_compute<dim,Np>>(g1,stencil_pnt);
	}

	/*! \brief Return a grid iterator over all points included between start and stop point
	 *
	 * Return a grid iterator over all the point with the exception of the
	 * ghost part
	 *
	 * \param start point
	 * \param stop point
	 * \param to_init unused bool
	 *
	 * \return a sub-grid iterator
	 *
	 */
	inline grid_key_dx_iterator_sub<dim> getIterator(const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop, bool to_init = false) const
	{
		size_t sz[dim];

		for (int i = 0 ; i < dim ; i++)
		{sz[i] = g1.size(i);}

		grid_sm<dim,void> gvoid(sz);

		// get the starting point and the end point of the real domain
		return grid_key_dx_iterator_sub<dim>(gvoid,start,stop);
	}

	/*! \brief return the internal data_
	 *
	 * return the internal data_
	 *
	 */
	layout & get_internal_data_()
	{
		return data_;
	}

	/*! \brief return the internal data_
	 *
	 * return the internal data_
	 *
	 */
	const layout & get_internal_data_() const
	{
		return data_;
	}

	/*! \brief In this case it does nothing
	 *
	 * \note this function exist to respect the interface to work as distributed
	 *
	 */
	void removeAddUnpackReset()
	{}

	/*! \brief In this case it does nothing
	 *
	 * \note this function exist to respect the interface to work as distributed
	 *
	 * \param ctx context
	 *
	 */
	template<unsigned int ... prp, typename context_type>
	void removeAddUnpackFinalize(const context_type & ctx, int opt)
	{}

	/*! \brief In this case it does nothing
	 *
	 * \note this function exist to respect the interface to work as distributed
	 *
	 * \param ctx context
	 *
	 */
	template<unsigned int ... prp, typename context_type>
	void removeCopyToFinalize(const context_type & ctx, int opt)
	{}

	/*! \brief It does nothing
	 *
	 *
	 */
	void resetFlush()
	{}
};


#endif /* OPENFPM_DATA_SRC_GRID_GRID_BASE_IMPLEMENTATION_HPP_ */
