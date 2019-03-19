#ifndef MAP_HPP_
#define MAP_HPP_

#include "config.h"

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
#ifdef SE_CLASS2
#include "Memleak_check.hpp"
#endif
#include "util/for_each_ref.hpp"
#include "util.hpp"
#include <utility>
#ifdef CUDA_GPU
#include "memory/CudaMemory.cuh"
#endif
#include "grid_sm.hpp"
#include "Encap.hpp"
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

#ifndef CUDA_GPU
typedef HeapMemory CudaMemory;
#endif


/*! Stub grid class
 *
 */
template<unsigned int dim, typename T, typename S=HeapMemory, typename layout = typename memory_traits_lin<T>::type >
class grid_cpu
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
template<unsigned int dim, typename T, typename S>
class grid_cpu<dim,T,S,typename memory_traits_lin<T>::type> : public grid_base_impl<dim,T,S,typename memory_traits_lin<T>::type, memory_traits_lin>
{
	typedef typename apply_transform<memory_traits_inte,T>::type T_;

public:

	//! type of layout of the structure
	typedef typename memory_traits_lin<T>::type layout;

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef typename grid_base_impl<dim,T,S,typename memory_traits_lin<T>::type, memory_traits_lin>::container container;

	//! Grid_cpu has no grow policy
	typedef void grow_policy;


	//! Default constructor
	inline grid_cpu() THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_lin>()
	{}

	/*! \brief create a grid from another grid
	 *
	 * \tparam S memory type for allocation
	 *
	 * \param g the grid to copy
	 * \param mem memory object (only used for template deduction)
	 *
	 */
	inline grid_cpu(const grid_cpu<dim,T,S,typename memory_traits_lin<T>::type> & g) THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_lin>(g)
	{
	}

	/*! \brief create a grid of size sz on each direction
	 *
	 * \param sz size if the grid on each directions
	 *
	 */
	inline grid_cpu(const size_t & sz) THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_lin>(sz)
	{
	}

	/*! \brief Constructor allocate memory
	 *
	 * \param sz size of the grid in each dimension
	 *
	 */
	inline grid_cpu(const size_t (& sz)[dim]) THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_lin>(sz)
	{
	}

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 */
	grid_cpu<dim,T,S,typename memory_traits_lin<T>::type> & operator=(const grid_base_impl<dim,T,S,layout,memory_traits_lin> & g)
	{
		(static_cast<grid_base_impl<dim,T,S,typename memory_traits_lin<T>::type, memory_traits_lin> *>(this))->swap(g.duplicate());

		return *this;
	}

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 */
	grid_cpu<dim,T,S,typename memory_traits_lin<T>::type> & operator=(grid_base_impl<dim,T,S,layout,memory_traits_lin> && g)
	{
		(static_cast<grid_base_impl<dim,T,S,typename memory_traits_lin<T>::type, memory_traits_lin> *>(this))->swap(g);

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

	/*! \brief Copy the memory from host to device
	 *
	 * \tparam (all properties are copied to prp is useless in this case)
	 *
	 */
	template<unsigned int ... prp> void hostToDevice()
	{
		this->data_.mem->hostToDevice();
	}

	/*! \brief Copy the memory from host to device
	 *
	 * \tparam (all properties are copied to prp is useless in this case)
	 *
	 * \param start start point
	 * \param stop stop point
	 *
	 */
	template<unsigned int ... prp> void hostToDevice(size_t start, size_t stop)
	{
		this->data_.mem->hostToDevice(start*sizeof(T),(stop+1)*sizeof(T));
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
		return this->data_.mem->getDevicePointer();
	}

	/*! \brief Synchronize the memory buffer in the device with the memory in the host
	 *
	 * \tparam ingored
	 *
	 * All properties are transfered
	 *
	 */
	template<unsigned int ... prp> void deviceToHost()
	{
		this->data_.mem->deviceToHost();
	}

	/*! \brief Synchronize the memory buffer in the device with the memory in the host
	 *
	 * \param start starting element to transfer
	 * \param stop stop element to transfer
	 *
	 * \tparam properties to transfer (ignored all properties are trasfert)
	 *
	 */
	template<unsigned int ... prp> void deviceToHost(size_t start, size_t stop)
	{
		this->data_.mem->deviceToHost(start*sizeof(T),(stop+1)*sizeof(T));
	}

#ifdef CUDA_GPU

	/*! \brief Convert the grid into a data-structure compatible for computing into GPU
	 *
	 *  The object created can be considered like a reference of the original
	 *
	 */
	grid_gpu_ker<dim,T_,memory_traits_lin> toKernel()
	{
		return grid_toKernelImpl<is_layout_inte<memory_traits_lin<T_>>::value,dim,T_>::toKernel(*this);
	}

	/*! \brief Convert the grid into a data-structure compatible for computing into GPU
	 *
	 *  The object created can be considered like a reference of the original
	 *
	 */
	const grid_gpu_ker<dim,T_,memory_traits_lin> toKernel() const
	{
		return grid_toKernelImpl<is_layout_inte<memory_traits_lin<T_>>::value,dim,T_>::toKernel(*this);
	}

#endif
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
template<typename T_type, template<typename> class layout_base , typename Memory, unsigned int ... prp>
struct host_to_device_impl
{
	//! encapsulated destination object
	typename memory_traits_inte<T_type>::type & dst;

	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	//! starting element
	size_t start;

	//! stop element
	size_t stop;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline host_to_device_impl(typename memory_traits_inte<T_type>::type & dst,size_t start, size_t stop)
	:dst(dst),start(start),stop(stop)
	{};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type ele_type;

		typedef decltype(boost::fusion::at_c<ele_type::value>(dst).mem_r) mem_r_type;

		typedef typename boost::mpl::at<typename T_type::type,ele_type>::type type_prp;

		typedef typename toKernel_transform<layout_base,typename mem_r_type::value_type>::type kernel_type;

		typedef boost::mpl::int_<(is_vector<typename mem_r_type::value_type>::value ||
								  is_vector_dist<typename mem_r_type::value_type>::value ||
								  is_gpu_celllist<typename mem_r_type::value_type>::value) + 2*std::is_array<type_prp>::value + std::rank<type_prp>::value> crh_cond;

		call_recursive_host_device_if_vector<typename mem_r_type::value_type,
											 kernel_type,
											 type_prp,
											 layout_base,
											 crh_cond::value>
		::template transform<Memory,mem_r_type>(static_cast<Memory *>(boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).mem),
									 boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).mem_r,
				                       start*sizeof(type_prp),
				                       (stop+1)*sizeof(type_prp));

		// here we have to recursively call hostToDevice for each nested vector
		call_recursive_host_device_if_vector<typename mem_r_type::value_type,
											 kernel_type,
											 type_prp,
											 layout_base,
											 0>
		::call(boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).mem_r,start,stop);
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
struct device_to_host_impl
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
	inline device_to_host_impl(typename memory_traits_inte<T_type>::type & dst)
	:dst(dst)
	{
	};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).mem->deviceToHost();
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
struct device_to_host_start_stop_impl
{
	//! encapsulated destination object
	typename memory_traits_inte<T_type>::type & dst;

	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	//! start
	size_t start;

	//! stop
	size_t stop;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline device_to_host_start_stop_impl(typename memory_traits_inte<T_type>::type & dst,size_t start,size_t stop)
	:dst(dst),start(start),stop(stop)
	{
	};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<typename T_type::type,T>::type p_type;

		boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).mem->deviceToHost(start*sizeof(p_type),(stop+1)*sizeof(p_type));
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
	typedef typename apply_transform<memory_traits_inte,T>::type T_;

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

	/*! \brief Fill the memory with a byte
	 *
	 */
	template<unsigned int id> void fill(unsigned char c)
	{
		boost::fusion::at_c<id>(this->data_).mem->fill(c);
	}

	/*! \brief Copy the memory from host to device
	 *
	 */
	template<unsigned int ... prp> void hostToDevice()
	{
		host_to_device_impl<T,memory_traits_inte,S,prp ...> htd(this->data_,0,this->getGrid().size()-1);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(htd);
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
		return boost::fusion::at_c<id>(this->data_).mem->getDevicePointer();
	}

	/*! \brief Synchronize the memory buffer in the device with the memory in the host
	 *
	 *
	 */
	template<unsigned int ... prp> void deviceToHost()
	{
		device_to_host_impl<T, prp ...> dth(this->data_);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(dth);
	}

	/*! \brief Synchronize the memory buffer in the device with the memory in the host
	 *
	 * \param start starting element to transfer
	 * \param stop stop element to transfer
	 *
	 * \tparam properties to transfer
	 *
	 */
	template<unsigned int ... prp> void deviceToHost(size_t start, size_t stop)
	{
		device_to_host_start_stop_impl<T, prp ...> dth(this->data_,start,stop);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(dth);
	}

	/*! \brief Synchronize the memory buffer in the device with the memory in the host
	 *
	 * \param start starting element to transfer
	 * \param stop stop element to transfer
	 *
	 * \tparam properties to transfer
	 *
	 */
	template<unsigned int ... prp> void hostToDevice(size_t start, size_t stop)
	{
		host_to_device_impl<T,memory_traits_inte,S, prp ...> dth(this->data_,start,stop);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(dth);
	}

#ifdef CUDA_GPU

	/*! \brief Convert the grid into a data-structure compatible for computing into GPU
	 *
	 *  The object created can be considered like a reference of the original
	 *
	 */
	grid_gpu_ker<dim,T_,memory_traits_inte> toKernel()
	{
		return grid_toKernelImpl<is_layout_inte<memory_traits_inte<T_>>::value,dim,T_>::toKernel(*this);
	}

	/*! \brief Convert the grid into a data-structure compatible for computing into GPU
	 *
	 *  The object created can be considered like a reference of the original
	 *
	 */
	const grid_gpu_ker<dim,T_,memory_traits_inte> toKernel() const
	{
		return grid_toKernelImpl<is_layout_inte<memory_traits_inte<T>>::value,dim,T_>::toKernel(*this);
	}

#endif
};

//! short formula for a grid on gpu
template <unsigned int dim, typename T> using grid_gpu = grid_cpu<dim,T,CudaMemory,typename memory_traits_inte<T>::type>;

#endif


