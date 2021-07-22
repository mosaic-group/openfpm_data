/*
 * grid_zm.hpp
 *
 *  Created on: Mar 25, 2020
 *      Author: i-bird
 */

#ifndef GRID_ZM_HPP_
#define GRID_ZM_HPP_

#include "util/zmorton.hpp"

/*! \brief class that store the information of the grid like number of point on each direction and
 *  define the index linearization by stride
 *
 * \param N dimensionality
 * \param T type of object is going to store the grid
 *
 */
template<unsigned int N, typename T>
class grid_zm : private grid_sm<N,T>
{

public:


	/*! \brief Reset the dimension of the grid
	 *
	 * \param dims store on each dimension the size of the grid
	 *
	 */
	inline void setDimensions(const size_t  (& dims)[N])
	{
		((grid_sm<N,T> *)this)->setDimensions(dims);
	}

	grid_zm(){};

	/*! \brief construct a grid from another grid
	 *
	 * \param g grid info
	 *
	 * construct a grid from another grid, type can be different
	 *
	 */

	template<typename S> inline grid_zm(const grid_zm<N,S> & g)
	{
		this->setDimensions(this->getSize());
	}


	/*! \brief Construct a grid of a specified size
	 *
	 * Construct a grid of a specified size
	 *
	 * \param sz is an array that contain the size of the grid on each dimension
	 *
	 */

	inline grid_zm(const size_t & sz)
	:grid_sm<N,T>(sz)
	{}

	/*! \brief Construct a grid of a specified size
	 *
	 * Construct a grid of a specified size
	 *
	 * \param sz is an array that contain the size of the grid on each dimension
	 *
	 */

	inline grid_zm(const size_t (& sz)[N])
	{
		this->setDimensions(sz);
	}

	//! Destructor
	~grid_zm() {};

	/*! \brief Linearization of the grid_key_dx
	 *
	 * Linearization of the grid_key_dx given a key, it spit out a number that is just the 1D linearization
	 * of the key. In this case is the linearization of N index
	 *
	 * \param gk grid key to access the element of the grid
	 *
	 */
	template<typename ids_type> __device__ __host__ inline mem_id LinId(const grid_key_dx<N,ids_type> & gk) const
	{
		return lin_zid(gk);
	}


	/*! \brief Copy the grid from another grid
	 *
	 * \param g grid from witch to copy
	 *
	 */

	__device__ __host__ inline grid_zm<N,T> & operator=(const grid_zm<N,T> & g)
	{
		((grid_sm<N,T> *)this)->operator=(g);

		return *this;
	}

	/*! \brief Check if the two grid_sm are the same
	 *
	 * \param g element to check
	 *
	 * \return true if they are the same
	 *
	 */

	inline bool operator==(const grid_zm<N,T> & g)
	{
		return ((grid_sm<N,T> *)this)->operator==(g);
	}

	/*! \brief Check if the two grid_sm are the same
	 *
	 * \param g element to check
	 *
	 */

	inline bool operator!=(const grid_zm<N,T> & g)
	{
		return ((grid_sm<N,T> *)this)->operator!=(g);
	}

	/*! \brief swap the grid_sm informations
	 *
	 * \param g grid to swap
	 *
	 */
	inline void swap(grid_zm<N,T> & g)
	{
		((grid_sm<N,T> *)this)->swap(g);
	}

	/**
	 *
	 * Get the size of the grid on the direction i
	 *
	 * \param i direction
	 * \return the size on the direction i
	 *
	 */
	inline size_t size(unsigned int i) const
	{
		return ((grid_sm<N,T> *)this)->size(i);
	}

	/**
	 *
	 * Get the total size of the grid
	 *
	 * \return the total size on the grid
	 *
	 */
	inline size_t size() const
	{
		return ((grid_sm<N,T> *)this)->size();
	}

	//!  It simply mean that all the classes grid are friend of all its specialization
	template <unsigned int,typename> friend class grid_zm;
};


#endif /* GRID_ZM_HPP_ */
