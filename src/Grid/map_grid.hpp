#ifndef MAP_HPP_
#define MAP_HPP_


#include "config.h"
#include "util/cuda_launch.hpp"
#include "util/object_util.hpp"
#include "Grid/util.hpp"
#include "Vector/vect_isel.hpp"
#include "Vector/util.hpp"
#include "Vector/map_vector_grow_p.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "util/util_debug.hpp"
#include "util/Pack_stat.hpp"
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include <boost/type_traits.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>
#include "memory_ly/memory_conf.hpp"
#include "util/copy_compare/meta_copy.hpp"
#include "util/for_each_ref.hpp"
#include "util.hpp"
#include <utility>
#ifdef CUDA_GPU
#include "memory/CudaMemory.cuh"
#endif
#include "grid_sm.hpp"
#include "grid_zm.hpp"
#include "memory_ly/Encap.hpp"
#include "memory_ly/memory_array.hpp"
#include "memory_ly/memory_c.hpp"
#include <vector>
#include "se_grid.hpp"
#include "memory/HeapMemory.hpp"
#include "memory/PtrMemory.hpp"
#include "grid_common.hpp"
#include "util/se_util.hpp"
#include "iterators/grid_key_dx_iterator.hpp"
#include "iterators/grid_key_dx_iterator_sub.hpp"
#include "iterators/grid_key_dx_iterator_sp.hpp"
#include "iterators/grid_key_dx_iterator_sub_bc.hpp"
#include "Packer_Unpacker/Packer_util.hpp"
#include "Packer_Unpacker/has_pack_agg.hpp"
#include "cuda/cuda_grid_gpu_funcs.cuh"
#include "grid_base_implementation.hpp"
#include "util/for_each_ref.hpp"
#include "Geometry/grid_smb.hpp"
#include "Geometry/grid_zmb.hpp"

#ifndef CUDA_GPU
typedef HeapMemory CudaMemory;
#endif

/*! \brief get the type of the SetBlock
 *
 *
 */
template<typename SGridGpu>
struct GetSetBlockType
{
	typedef typename SGridGpu::device_grid_type::container type;
};

/*! Stub grid class
 *
 */
template<unsigned int dim, typename T, typename S=HeapMemory, typename layout = typename memory_traits_lin<T>::type, typename linearizer = grid_sm<dim,void>  >
class grid_base
{
};

/*!
 *
 * \brief This is an N-dimensional grid or an N-dimensional array with memory_traits_lin layout
 *
 * it is basically an N-dimensional Cartesian grid
 *
 *	\tparam dim Dimensionality of the grid
 *	\tparam T type of object the grid store
 *	\tparam S type of memory HeapMemory CudaMemory
 *	\tparam layout memory layout
 *
 * ### Defining the grid size on each dimension
 *
 * \code{.cpp}
 *  size_t sz[3] = {16,16,16};
 * \endcode
 *
 * ### Definition and allocation of a 3D grid on CPU memory
 * \snippet grid_unit_tests.hpp Definition and allocation of a 3D grid on CPU memory
 * ### Access a grid c3 of size sz on each direction
 * \snippet grid_unit_tests.hpp Access a grid c3 of size sz on each direction
 * ### Access an N-dimensional grid with an iterator
 * \snippet grid_unit_tests.hpp Access to an N-dimensional grid with an iterator
 * ### Iterate only on a sub-set of the grid
 * \snippet grid_unit_tests.hpp Sub-grid iterator test usage
 * ### Get the full-object in an N-dimensional grid
 * \snippet grid_unit_tests.hpp Get the object in an N-dimensional grid with an iterator
 * ### Create a grid g1 and copy into another g2
 * \snippet grid_unit_tests.hpp Create a grid g1 and copy into another g2
 *
 */
template<unsigned int dim, typename T, typename S, typename linearizer>
class grid_base<dim,T,S,typename memory_traits_lin<T>::type, linearizer> : public grid_base_impl<dim,T,S, memory_traits_lin,linearizer>
{
	typedef typename apply_transform<memory_traits_lin,T>::type T_;

	T background;

public:

	//! type of layout of the structure
	typedef typename memory_traits_lin<T>::type layout;

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef typename grid_base_impl<dim,T,S, memory_traits_lin>::container container;

	//! grid_base has no grow policy
	typedef void grow_policy;

	//! type that identify one point in the grid
	typedef grid_key_dx<dim> base_key;

	//! sub-grid iterator type
	typedef grid_key_dx_iterator_sub<dim> sub_grid_iterator_type;

	//! linearizer type Z-morton Hilbert curve , normal striding
	typedef typename grid_base_impl<dim,T,S, memory_traits_lin>::linearizer_type linearizer_type;

	//! Default constructor
	inline grid_base() THROW
	:grid_base_impl<dim,T,S,memory_traits_lin, linearizer>()
	{}

	/*! \brief create a grid from another grid
	 *
	 * \tparam S memory type for allocation
	 *
	 * \param g the grid to copy
	 * \param mem memory object (only used for template deduction)
	 *
	 */
	inline grid_base(const grid_base<dim,T,S,typename memory_traits_lin<T>::type> & g) THROW
	:grid_base_impl<dim,T,S,memory_traits_lin, linearizer>(g)
	{
	}

	/*! \brief create a grid of size sz on each direction
	 *
	 * \param sz size if the grid on each directions
	 *
	 */
	inline grid_base(const size_t & sz) THROW
	:grid_base_impl<dim,T,S,memory_traits_lin,linearizer>(sz)
	{
	}

	/*! \brief Constructor allocate memory
	 *
	 * \param sz size of the grid in each dimension
	 *
	 */
	inline grid_base(const size_t (& sz)[dim]) THROW
	:grid_base_impl<dim,T,S,memory_traits_lin,linearizer>(sz)
	{
	}

	/*! \brief Stub does not do anything
	*
	*/
	template<typename pointers_type, 
			 typename headers_type, 
			 typename result_type, 
			 unsigned int ... prp >
	static void unpack_headers(pointers_type & pointers, headers_type & headers, result_type & result, int n_slot)
	{}

	template<unsigned int ... prp, typename S2, typename header_type, typename ite_type, typename context_type>
	void unpack_with_headers(ExtPreAlloc<S2> & mem,
				ite_type & sub_it,
				header_type & headers,
				int ih,
				Unpack_stat & ps,
				context_type &gpuContext,
				rem_copy_opt opt = rem_copy_opt::NONE_OPT)
	{}

#if defined(__HIP__)

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 */
	__device__ grid_base<dim,T,S> & operator=(const grid_base<dim,T,S> & g)
	{
		printf("Error grid_base operator= is not defined in device code\n");

		return *this;
	}

#endif

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 */
	__host__ grid_base<dim,T,S> & operator=(const grid_base<dim,T,S> & g)
	{
		(static_cast<grid_base_impl<dim,T,S, memory_traits_lin> *>(this))->swap(g.duplicate());

		meta_copy<T>::meta_copy_(g.background,background);

		return *this;
	}

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 */
	grid_base<dim,T,S,typename memory_traits_lin<T>::type> & operator=(grid_base<dim,T,S,typename memory_traits_lin<T>::type> && g)
	{
		(static_cast<grid_base_impl<dim,T,S, memory_traits_lin> *>(this))->swap(g);

		meta_copy<T>::meta_copy_(g.background,background);

		return *this;
	}

	/*! \brief This structure has pointers
	 *
	 * \return false
	 *
	 */
	static bool noPointers()
	{
		return false;
	}

	/*! \brief It return the properties arrays.
	 *
	 * In case of Cuda memory it return the device pointers to pass to the kernels
	 *
	 * This variant does not copy the host memory to the device memory
	 *
	 */
	template<unsigned int id> void * getDeviceBuffer()
	{
		return ((S*)this->data_.mem)->getDevicePointer();
	}

	/*! \brief This is a meta-function return which type of sub iterator a grid produce
	 *
	 * \return the type of the sub-grid iterator
	 *
	 */
	template <typename stencil = no_stencil>
	static grid_key_dx_iterator_sub<dim, stencil> type_of_subiterator()
	{
		return grid_key_dx_iterator_sub<dim, stencil>();
	}

	/*! \brief Return if in this representation data are stored is a compressed way
	 *
	 * \return false this is a normal grid no compression
	 *
	 */
	static constexpr bool isCompressed()
	{
		return false;
	}

	/*! \brief This is a meta-function return which type of iterator a grid produce
	 *
	 * \return the type of the sub-grid iterator
	 *
	 */
	static grid_key_dx_iterator<dim> type_of_iterator()
	{
		return grid_key_dx_iterator<dim>();
	}

	/*! \brief In this case it just copy the key_in in key_out
	 *
	 * \param key_out output key
	 * \param key_in input key
	 *
	 */
	void convert_key(grid_key_dx<dim> & key_out, const grid_key_dx<dim> & key_in) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{key_out.set_d(i,key_in.get(i));}
	}

	/*! \brief Get the background value
	 *
	 * For dense grid this function is useless
	 *
	 * \return background value
	 *
	 */
	T & getBackgroundValue()
	{
		return background;
	}

	/*! \brief Get the background value
	 *
	 * For dense grid this function is useless
	 *
	 * \return background value
	 *
	 */
	T & getBackgroundValueAggr()
	{
		return background;
	}

	/*! \brief Set the background value
	 *
	 * \tparam p property to set
	 *
	 */
	template<unsigned int p>
	void setBackgroundValue(const typename boost::mpl::at<typename T::type,boost::mpl::int_<p>>::type & val)
	{
		meta_copy<typename boost::mpl::at<typename T::type,boost::mpl::int_<p>>::type>::meta_copy_(val,background.template get<p>());
	}


	/*! \brief assign operator
	 *
	 * \return itself
	 *
	 */
	grid_base<dim,T,S,typename memory_traits_lin<T>::type> & operator=(const grid_base_impl<dim,T,S, memory_traits_lin> & base)
	{
		grid_base_impl<dim,T,S, memory_traits_inte>::operator=(base);

		return *this;
	}

	/*! \brief assign operator
	 *
	 * \return itself
	 *
	 */
	grid_base<dim,T,S,typename memory_traits_lin<T>::type> & operator=(grid_base_impl<dim,T,S, memory_traits_lin> && base)
	{
		grid_base_impl<dim,T,S, memory_traits_lin>::operator=((grid_base_impl<dim,T,S, memory_traits_lin> &&)base);

		return *this;
	}
};



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
template<typename T_type, unsigned int ... prp>
struct switch_copy_host_to_device
{
	//! encapsulated destination object
	typename memory_traits_inte<T_type>::type & dst;

	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline switch_copy_host_to_device(typename memory_traits_inte<T_type>::type & dst)
	:dst(dst)
	{
	};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).switchToDevicePtr();
	}
};


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
template<typename T_type, template<typename> class layout_base , typename Memory>
struct deconstruct_impl
{
	//! object to destruct
	typename memory_traits_inte<T_type>::type & dst;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline deconstruct_impl(typename memory_traits_inte<T_type>::type & dst)
	:dst(dst)
	{};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef decltype(boost::fusion::at_c<T::value>(dst).mem_r) mem_r_type;

		typedef typename boost::mpl::at<typename T_type::type,T>::type type_prp;

		typedef typename toKernel_transform<layout_base,typename mem_r_type::value_type>::type kernel_type;

		typedef boost::mpl::int_<(is_vector<typename mem_r_type::value_type>::value ||
								  is_vector_dist<typename mem_r_type::value_type>::value ||
								  is_gpu_celllist<typename mem_r_type::value_type>::value) + 2*std::is_array<type_prp>::value + std::rank<type_prp>::value> crh_cond;

		call_recursive_destructor_if_vector<typename mem_r_type::value_type,
											 kernel_type,
											 type_prp,
											 layout_base,
											 crh_cond::value>
		::template destruct<Memory,mem_r_type>(static_cast<Memory *>(boost::fusion::at_c<T::value>(dst).mem),
									 boost::fusion::at_c<T::value>(dst).mem_r);
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
template<unsigned int dim, typename T, typename S, typename linearizer>
class grid_base<dim,T,S,typename memory_traits_inte<T>::type,linearizer> : public grid_base_impl<dim,T,S, memory_traits_inte,linearizer>
{
	typedef typename apply_transform<memory_traits_inte,T>::type T_;

	T background;

public:

	//! grid layout
	typedef typename memory_traits_inte<T>::type layout;

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef typename grid_base_impl<dim,T,S, memory_traits_inte,linearizer>::container container;

	//! linearizer type Z-morton Hilbert curve , normal striding
	typedef typename grid_base_impl<dim,T,S, memory_traits_inte,linearizer>::linearizer_type linearizer_type;

	//! Default constructor
	inline grid_base() THROW
	:grid_base_impl<dim,T,S,memory_traits_inte,linearizer>()
	{
	}

	/*! \brief create a grid from another grid
	 *
	 * \param g the grid to copy
	 *
	 */
	inline grid_base(const grid_base & g) THROW
	:grid_base_impl<dim,T,S,memory_traits_inte,linearizer>(g)
	{
	}

	/*! \brief create a grid from another grid
	 *
	 * \param g the grid to copy
	 *
	 */
	inline grid_base(grid_base && g) THROW
	:grid_base_impl<dim,T,S,memory_traits_inte,linearizer>(g)
	{
	}

	/*! \brief create a grid of size sz on each direction
	 *
	 * \param sz grid size in each direction
	 *
	 */
	inline grid_base(const size_t & sz) THROW
	:grid_base_impl<dim,T,S,memory_traits_inte,linearizer>(sz)
	{
	}

	//! Constructor allocate memory and give them a representation
	inline grid_base(const size_t (& sz)[dim]) THROW
	:grid_base_impl<dim,T,S,memory_traits_inte,linearizer>(sz)
	{
	}

	/*! \brief Stub does not do anything
	*
	*/
	static void unpack_headers()
	{}

	/*! \brief Fill the memory with a byte
	 *
	 */
	template<unsigned int id> void fill(unsigned char c)
	{
		boost::fusion::at_c<id>(this->data_).mem->fill(c);
	}

	/*! \brief It return the properties arrays.
	 *
	 * In case of Cuda memory it return the device pointers to pass to the kernels
	 *
	 * This variant does not copy the host memory to the device memory
	 *
	 */
	template<unsigned int id> void * getDeviceBuffer()
	{
		return ((S*)boost::fusion::at_c<id>(this->data_).mem)->getDevicePointer();
	}

	/*! \brief This is a meta-function return which type of sub iterator a grid produce
	 *
	 * \return the type of the sub-grid iterator
	 *
	 */
	template <typename stencil = no_stencil>
	static grid_key_dx_iterator_sub<dim, stencil> type_of_subiterator()
	{
		return grid_key_dx_iterator_sub<dim, stencil>();
	}

	/*! \brief Return if in this representation data are stored is a compressed way
	 *
	 * \return false this is a normal grid no compression
	 *
	 */
	static constexpr bool isCompressed()
	{
		return false;
	}

	/*! \brief This is a meta-function return which type of iterator a grid produce
	 *
	 * \return the type of the sub-grid iterator
	 *
	 */
	static grid_key_dx_iterator<dim> type_of_iterator()
	{
		return grid_key_dx_iterator<dim>();
	}

	/*! \brief In this case it just copy the key_in in key_out
	 *
	 * \param key_out output key
	 * \param key_in input key
	 *
	 */
	void convert_key(grid_key_dx<dim> & key_out, const grid_key_dx<dim> & key_in) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{key_out.set_d(i,key_in.get(i));}
	}

	/*! \brief Get the background value
	 *
	 * For dense grid this function is useless
	 *
	 * \return background value
	 *
	 */
	T & getBackgroundValue()
	{
		return background;
	}

	/*! \brief Get the background value
	 *
	 * For dense grid this function is useless
	 *
	 * \return background value
	 *
	 */
	T & getBackgroundValueAggr()
	{
		return background;
	}

	/*! \brief assign operator
	 *
	 * \return itself
	 *
	 */
	grid_base<dim,T,S,typename memory_traits_inte<T>::type,linearizer> & operator=(const grid_base_impl<dim,T,S, memory_traits_inte,linearizer> & base)
	{
		grid_base_impl<dim,T,S, memory_traits_inte,linearizer>::operator=(base);

		return *this;
	}

	/*! \brief assign operator
	 *
	 * \return itself
	 *
	 */
	grid_base<dim,T,S,typename memory_traits_inte<T>::type,linearizer> & operator=(grid_base_impl<dim,T,S, memory_traits_inte,linearizer> && base)
	{
		grid_base_impl<dim,T,S, memory_traits_inte,linearizer>::operator=(base);

		return *this;
	}

	~grid_base()
	{
		deconstruct_impl<T,memory_traits_inte,S> dth(this->data_);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(dth);
	}
};

//! short formula for a grid on gpu
template <unsigned int dim, typename T, typename linearizer = grid_sm<dim,void> > using grid_gpu = grid_base<dim,T,CudaMemory,typename memory_traits_inte<T>::type,linearizer>;

//! short formula for a grid on gpu
template <unsigned int dim, typename T, typename linearizer = grid_sm<dim,void> > using grid_cpu = grid_base<dim,T,HeapMemory,typename memory_traits_lin<T>::type,linearizer>;


#endif


