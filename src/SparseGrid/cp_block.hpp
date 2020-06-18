/*
 * cp_block.hpp
 *
 *  Created on: Apr 20, 2020
 *      Author: i-bird
 */

#ifndef CP_BLOCK_HPP_
#define CP_BLOCK_HPP_

#include "util/create_vmpl_sequence.hpp"

template<typename T, unsigned int stencil_size, typename vector_vmpl, unsigned int dim>
class cp_block
{
	// we create first a vector with

	typedef typename vmpl_sum_constant<2*stencil_size,vector_vmpl>::type stop_border_vmpl;
	typedef typename vmpl_create_constant<dim,stencil_size>::type start_border_vmpl;

	static const int sizeBlock = vmpl_reduce_prod<vector_vmpl>::type::value;
	static const int sizeBlockBord = vmpl_reduce_prod<stop_border_vmpl>::type::value;

	T (& ptr)[sizeBlock];

public:

	__device__ __host__ cp_block(T (& ptr)[sizeBlock])
	:ptr(ptr)
	{}

	template<typename ... ArgsT>
	__device__ __host__  T & operator()(ArgsT ... args)
	{
		return ptr[Lin_vmpl_off<vector_vmpl,start_border_vmpl>(args ...)];
	}
};

template<typename T, unsigned int stencil_size, typename vector_vmpl>
class cp_block<T,stencil_size,vector_vmpl,2>
{
	// we create first a vector with

	typedef typename vmpl_sum_constant<2*stencil_size,vector_vmpl>::type stop_border_vmpl;
	typedef typename vmpl_create_constant<2,stencil_size>::type start_border_vmpl;

	static const int sizeBlock = vmpl_reduce_prod<vector_vmpl>::type::value;
	static const int sizeBlockBord = vmpl_reduce_prod<stop_border_vmpl>::type::value;

	T (& ptr)[sizeBlock];

public:

	__device__ __host__ cp_block(T (& ptr)[sizeBlock])
	:ptr(ptr)
	{}

	__device__ __host__ int LinId(int i , int j)
	{
		return Lin_vmpl_off<vector_vmpl,start_border_vmpl>(i,j);
	}

	__device__ __host__ T & operator()( int i , int j)
	{
		return ptr[Lin_vmpl_off<vector_vmpl,start_border_vmpl>(i,j)];
	}
};

template<typename T, unsigned int stencil_size, typename vector_vmpl>
class cp_block<T,stencil_size,vector_vmpl,3>
{
	// we create first a vector with

	typedef typename vmpl_sum_constant<2*stencil_size,vector_vmpl>::type stop_border_vmpl;
	typedef typename vmpl_create_constant<3,stencil_size>::type start_border_vmpl;

	static const int sizeBlock = vmpl_reduce_prod<vector_vmpl>::type::value;
	static const int sizeBlockBord = vmpl_reduce_prod<stop_border_vmpl>::type::value;

	T (& ptr)[sizeBlock];

public:

	__device__ __host__ cp_block(T (& ptr)[sizeBlock])
	:ptr(ptr)
	{}

	__device__ __host__ int LinId(int i , int j, int k)
	{
		return Lin_vmpl_off<vector_vmpl,start_border_vmpl>(i,j,k);
	}

	__device__ __host__ T & operator()(int i , int j, int k)
	{
		return ptr[Lin_vmpl_off<vector_vmpl,start_border_vmpl>(i,j,k)];
	}
};

#endif /* CP_BLOCK_HPP_ */
