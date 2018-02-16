#ifndef MAP_HPP_
#define MAP_HPP_

#include "config.h"

//! Warning: apparently you cannot used nested boost::mpl with boost::fusion
//! can create template circularity, this include avoid the problem
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

	T background;

public:

	//! type of layout of the structure
	typedef typename memory_traits_lin<T>::type layout;

	//! Object container for T, it is the return type of get_o it return a object type trough
	// you can access all the properties of T
	typedef typename grid_base_impl<dim,T,S,typename memory_traits_lin<T>::type, memory_traits_lin>::container container;

	//! type that identify one point in the grid
	typedef grid_key_dx<dim> base_key;

	//! sub-grid iterator type
	typedef grid_key_dx_iterator_sub<dim> sub_grid_iterator_type;

	//! Default constructor
	inline grid_cpu() THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_lin>()
	{
	}

	/*! \brief create a grid from another grid
	 *
	 * \tparam S memory type for allocation
	 *
	 * \param g the grid to copy
	 * \param mem memory object (only used for template deduction)
	 *
	 */
	inline grid_cpu(const grid_cpu & g) THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_lin>(g)
	{
	}

	/*! \brief create a grid of size sz on each direction
	 *
	 * \tparam S memory type for allocation
	 *
	 * \param sz size if the grid on each directions
	 *
	 */
	inline grid_cpu(const size_t & sz) THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_lin>(sz)
	{
	}

	//! Constructor allocate memory and give them a representation
	inline grid_cpu(const size_t (& sz)[dim]) THROW
	:grid_base_impl<dim,T,S,layout,memory_traits_lin>(sz)
	{
	}

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 */
	grid_cpu<dim,T,S,typename memory_traits_lin<T>::type> & operator=(const grid_cpu<dim,T,S,typename memory_traits_lin<T>::type> & g)
	{
		(static_cast<grid_base_impl<dim,T,S,typename memory_traits_lin<T>::type, memory_traits_lin> *>(this))->swap(g.duplicate());

		meta_copy<T>::meta_copy_(g.background,background);

		return *this;
	}

	/*! \brief It copy a grid
	 *
	 * \param g grid to copy
	 *
	 */
	grid_cpu<dim,T,S,typename memory_traits_lin<T>::type> & operator=(grid_cpu<dim,T,S,typename memory_traits_lin<T>::type> && g)
	{
		(static_cast<grid_base_impl<dim,T,S,typename memory_traits_lin<T>::type, memory_traits_lin> *>(this))->swap(g);

		meta_copy<T>::meta_copy_(g.background,background);

		return *this;
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
};

#include "grid_gpu.hpp"

#endif


