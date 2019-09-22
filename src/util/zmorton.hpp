/*
 * zmorton.hpp
 *
 *  Created on: Jul 31, 2019
 *      Author: i-bird
 */

#ifndef ZMORTON_HPP_
#define ZMORTON_HPP_

#include "Grid/grid_key.hpp"

template<typename T>
inline __device__ __host__ size_t lin_zid(const grid_key_dx<1,T> & key)
{
	return key.get(0);
}

template<typename T>
inline __device__ __host__ void invlin_zid(size_t lin, grid_key_dx<1,T> & key)
{
	return key.set_d(0,lin);
}


template<typename T>
inline __device__ __host__  size_t lin_zid(const grid_key_dx<2,T> & key)
{
	size_t x = key.get(0);
	size_t y = key.get(1);


	x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
	x = (x | (x << 8)) & 0x00FF00FF00FF00FF;
	x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0F;
	x = (x | (x << 2)) & 0x3333333333333333;
	x = (x | (x << 1)) & 0x5555555555555555;

	y = (y | (y << 16)) & 0x0000FFFF0000FFFF;
	y = (y | (y << 8)) & 0x00FF00FF00FF00FF;
	y = (y | (y << 4)) & 0x0F0F0F0F0F0F0F0F;
	y = (y | (y << 2)) & 0x3333333333333333;
	y = (y | (y << 1)) & 0x5555555555555555;

	return x | (y << 1);
}

template<typename T>
inline __device__ __host__  void invlin_zid(size_t lin, grid_key_dx<2,T> & key)
{
	size_t x = lin & 0x5555555555555555;
	size_t y = (lin & 0xAAAAAAAAAAAAAAAA) >> 1;

	x = (x | (x >> 1)) & 0x3333333333333333;
	x = (x | (x >> 2)) & 0x0F0F0F0F0F0F0F0F;
	x = (x | (x >> 4)) & 0x00FF00FF00FF00FF;
	x = (x | (x >> 8)) & 0x0000FFFF0000FFFF;
	x = (x | (x >> 16)) & 0x00000000FFFFFFFF;

	y = (y | (y >> 1)) & 0x3333333333333333;
	y = (y | (y >> 2)) & 0x0F0F0F0F0F0F0F0F;
	y = (y | (y >> 4)) & 0x00FF00FF00FF00FF;
	y = (y | (y >> 8)) & 0x0000FFFF0000FFFF;
	y = (y | (y >> 16)) & 0x00000000FFFFFFFF;

	key.set_d(0,x);
	key.set_d(1,y);
}


static const size_t S3[] = {2, 4, 8, 16, 32};

template<typename T>
inline __device__ __host__  size_t lin_zid(const grid_key_dx<3,T> & key)
{
	size_t x = key.get(0);
	size_t z = key.get(2);
	size_t y = key.get(1);

	x = (x | (x << 32)) & 0xFFFF0000FFFFFFFF;
	x = (x | (x << 16)) & 0x0FFF000FFF000FFF;
	x = (x | (x << 8)) & 0xF00F00F00F00F00F;
	x = (x | (x << 4)) & 0x30C30C30C30C30C3;
	x = (x | (x << 2)) & 0x9249249249249249;

	y = (y | (y << 32)) & 0xFFFF0000FFFFFFFF;
	y = (y | (y << 16)) & 0x0FFF000FFF000FFF;
	y = (y | (y << 8)) & 0xF00F00F00F00F00F;
	y = (y | (y << 4)) & 0x30C30C30C30C30C3;
	y = (y | (y << 2)) & 0x9249249249249249;

	z = (z | (z << 32)) & 0xFFFF0000FFFFFFFF;
	z = (z | (z << 16)) & 0x0FFF000FFF000FFF;
	z = (z | (z << 8)) & 0xF00F00F00F00F00F;
	z = (z | (z << 4)) & 0x30C30C30C30C30C3;
	z = (z | (z << 2)) & 0x9249249249249249;

	return x | (y << 1) | (z << 2);
}

template<typename T>
inline __device__ __host__  void invlin_zid(size_t lin, grid_key_dx<3,T> & key)
{
	size_t x = lin & 0x9249249249249249;
	size_t y = (lin >> 1) & 0x9249249249249249;
	size_t z = (lin >> 2) & 0x9249249249249249;

	x = (x | (x >> 2)) & 0x30C30C30C30C30C3;
	x = (x | (x >> 4)) & 0xF00F00F00F00F00F;
	x = (x | (x >> 8)) & 0x00FF0000FF0000FF;
	x = (x | (x >> 16)) & 0x00000FF0000FFFF;
	x = (x | x >> 16) & 0xFFFFFF;

	y = (y | (y >> 2)) & 0x30C30C30C30C30C3;
	y = (y | (y >> 4)) & 0xF00F00F00F00F00F;
	y = (y | (y >> 8)) & 0x00FF0000FF0000FF;
	y = (y | (y >> 16)) & 0x00000FF0000FFFF;
	y = (y | y >> 16) & 0xFFFFFF;

	z = (z | (z >> 2)) & 0x30C30C30C30C30C3;
	z = (z | (z >> 4)) & 0xF00F00F00F00F00F;
	z = (z | (z >> 8)) & 0x00FF0000FF0000FF;
	z = (z | (z >> 16)) & 0x00000FF0000FFFF;
	z = (z | z >> 16) & 0xFFFFFF;

	key.set_d(0,x);
	key.set_d(1,y);
	key.set_d(2,z);
}

#endif /* ZMORTON_HPP_ */
