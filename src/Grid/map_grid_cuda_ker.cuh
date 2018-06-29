/*
 * map_grid_cuda_ker.hpp
 *
 *  Created on: Jun 28, 2018
 *      Author: i-bird
 */

#ifndef MAP_GRID_CUDA_KER_HPP_
#define MAP_GRID_CUDA_KER_HPP_


/*! \brief grid interface available when on gpu
 *
 * \tparam n_buf number of template buffers
 *
 */
template<unsigned int dim, typename T>
struct grid_gpu_ker
{
	//! grid information
	grid_sm<dim,void> g1;

	//! type of layout of the structure
	typedef typename memory_traits_inte<T>::type layout;

	//! layout data
	layout data_;

	grid_gpu_ker()
	{}

	grid_gpu_ker(const grid_sm<dim,void> & g1)
	:g1(g1)
	{
	}

	grid_gpu_ker(const grid_gpu_ker & cpy)
	:g1(cpy.g1)
	{
		copy_switch_memory_c_no_cpy<T> bp_mc(cpy.data_,this->data_);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(bp_mc);
	}

	/*! \brief Return the internal grid information
	 *
	 * Return the internal grid information
	 *
	 * \return the internal grid
	 *
	 */
	__device__ __host__ const grid_sm<dim,void> & getGrid() const
	{
		return g1;
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(mem_get<p,memory_traits_inte<T>,layout,grid_sm<dim,T>,grid_key_dx<dim>>::get(data_,g1,grid_key_dx<dim>()))>
	__device__ __host__ inline r_type get(const grid_key_dx<dim> & v1)
	{
		return mem_get<p,memory_traits_inte<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(data_,g1,v1);
	}

	/*! \brief Get the const reference of the selected element
	 *
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the const reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(mem_get<p,memory_traits_inte<T>,layout,grid_sm<dim,T>,grid_key_dx<dim>>::get(data_,g1,grid_key_dx<dim>()))>
	__device__ __host__ inline const r_type get(const grid_key_dx<dim> & v1) const
	{
		return mem_get<p,memory_traits_inte<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(data_,g1,v1);
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \param lin_id linearized element that identify the element in the grid
	 *
	 * \return the reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(mem_get<p,memory_traits_inte<T>,layout,grid_sm<dim,T>,grid_key_dx<dim>>::get_lin(data_,g1,0))>
	__device__ __host__ inline r_type get(const size_t lin_id)
	{
		return mem_get<p,memory_traits_inte<T>,decltype(this->data_),decltype(this->g1),grid_key_dx<dim>>::get_lin(data_,g1,lin_id);
	}

	/*! \brief Get the const reference of the selected element
	 *
	 * \param lin_id linearized element that identify the element in the grid
	 *
	 * \return the const reference of the element
	 *
	 */
	template <unsigned int p, typename r_type=decltype(mem_get<p,memory_traits_inte<T>,layout,grid_sm<dim,T>,grid_key_dx<dim>>::get_lin(data_,g1,0))>
	__device__ __host__ inline const r_type get(size_t lin_id) const
	{
		return mem_get<p,memory_traits_inte<T>,decltype(this->data_),decltype(this->g1),grid_key_dx<dim>>::get_lin(data_,g1,lin_id);
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
	__device__ inline encapc<dim,T,layout> get_o(const grid_key_dx<dim> & v1)
	{
		return mem_geto<dim,T,memory_traits_inte<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(data_,g1,v1);
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
	__device__ inline const encapc<dim,T,layout> get_o(const grid_key_dx<dim> & v1) const
	{
		return mem_geto<dim,T,memory_traits_inte<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(const_cast<decltype(this->data_) &>(data_),g1,v1);
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
//	__device__ inline encapc<dim,T,layout> get_o(int v1)
//	{
//		return mem_geto<dim,T,memory_traits_inte<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(data_,g1,v1);
//	}

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
//	__device__ inline const encapc<dim,T,layout> get_o(int v1) const
//	{
//		return mem_geto<dim,T,memory_traits_inte<T>,decltype(this->data_),decltype(this->g1),decltype(v1)>::get(const_cast<decltype(this->data_) &>(data_),g1,v1);
//	}

	/*! \brief Set an element of the grid from another element of another grid
	 *
	 * \param key1 element of the grid to set
	 * \param g source grid
	 * \param key2 element of the source grid to copy
	 *
	 */

	__device__ inline void set(const grid_key_dx<dim> & key1,const grid_gpu_ker<dim,T> & g, const grid_key_dx<dim> & key2)
	{
		this->get_o(key1) = g.get_o(key2);
	}
};


#endif /* MAP_GRID_CUDA_KER_HPP_ */
