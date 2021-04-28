/*
 * compute_optimal_device_grid.hpp
 *
 *  Created on: Oct 1, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_COMPUTE_OPTIMAL_DEVICE_GRID_HPP_
#define OPENFPM_DATA_SRC_UTIL_COMPUTE_OPTIMAL_DEVICE_GRID_HPP_


/*! \brief Calculate an optimal decomposition in grids and threads for GPU or parallel accelerators
 *
 * First it factorize the grid size in each direction in prime factors. It try to get
 * a multiple of two for the block size that does not pass the maximum block size.
 * If the block size result too small (smaller than min_block_size) it try to construct a block
 * non multiple of two but bigger than min_block_size and smaller than max_block_size.
 * The min_block_size is constrain weakly imposed. When impossible (or difficult) to find
 * a block size in the range the min_block_size constrain can be broken, but the max_block_size
 * is never broken
 *
 * the function guarantee that
 *
 * sz[0] = dg.threads.x*dg.grids.x
 * sz[1] = dg.threads.y*dg.grids.y
 * sz[2]*...sz[dim-1] = dg.threads.z*dg.grids.z
 *
 *
 *
 * \param dg device grid
 * \param sz size of the grid
 * \param max_block_size maximum block size
 * \param min_block_size minimum block size
 *
 */
template<unsigned int dim>
void calculate_optimal_device_grid(device_grid<dim> & dg,
								   size_t (& sz)[dim],
								   size_t max_block_size,
								   size_t min_block_size)
{

	if (dim == 0)	{return ;}

	size_t tot_block = 1;

	// Here we calculate the factors for each grid dimension and prioritize the
	// factors by 2 for the blocks

	// Get the factors for x
	std::vector<size_t> x;
	openfpm::math::getFactorization(sz[0],x);

	dg.threads.x = 1;
	size_t jx = 0;

	while(jx < x.size() && x[jx] == 2 && tot_block < max_block_size)
	{
		if (tot_block * 2 > max_block_size)
		{break;}

		dg.threads.x *= 2;
		tot_block *= 2;

		jx++;
	}

	// if we already reach the maximum block size we finished
	if (tot_block * 2 > max_block_size)
	{
		dg.threads.y = 1;
		dg.threads.z = 1;

		dg.grids.x = 1;
		for (; jx < x.size() ; jx++)
		{dg.grids.x *= x[jx];}

		if (dim >= 2)
		{dg.grids.y = sz[1];}
		else
		{dg.grids.y = 1;}

		dg.grids.z = 1;
		for (size_t k = 2 ; k < dim ; k++)
		// coverty[dead_error_line]
		{dg.grids.z *= sz[k];}

		return;
	}


	// Get the factors for y
	std::vector<size_t> y;
	size_t jy = 0;
	dg.threads.y = 1;

	if (dim >= 2)
	{
		openfpm::math::getFactorization(sz[1],y);

		while(jy < y.size() && y[jy] == 2 && tot_block < max_block_size)
		{
			if (tot_block * 2 > max_block_size)
			{break;}

			dg.threads.y *= 2;
			tot_block *= 2;

			jy++;
		}

		// if we already reach the maximum block size we finished
		if (tot_block * 2 > max_block_size)
		{
			dg.threads.z = 1;

			dg.grids.x = 1;
			for (; jx < x.size() ; jx++)
			{dg.grids.x *= x[jx];}

			dg.grids.y = 1;
			for (; jy < y.size() ; jy++)
			{dg.grids.y *= y[jy];}

			dg.grids.z = 1;
			for (size_t k = 2 ; k < dim ; k++)
			{dg.grids.z *= sz[k];}

			return;
		}
	}

	// Get the factors for z
	std::vector<size_t> z;

	size_t jz = 0;
	dg.threads.z = 1;

	if (dim >= 3)
	{
		openfpm::math::getFactorization(sz[2],z);

		while(jz < z.size() && z[jz] == 2 && tot_block < max_block_size)
		{
			if (tot_block * 2 > max_block_size)
			{break;}

			dg.threads.z *= 2;
			tot_block *= 2;

			jz++;
		}

		// if we already reach the maximum block size we finished
		if (tot_block * 2 > max_block_size)
		{
			dg.grids.x = 1;
			for (; jx < x.size() ; jx++)
			{dg.grids.x *= x[jx];}

			dg.grids.y = 1;
			for (; jy < y.size() ; jy++)
			{dg.grids.y *= y[jy];}

			dg.grids.z = 1;
			for (; jz < z.size() ; jz++)
			{dg.grids.z *= z[jz];}

			for (size_t k = 3 ; k < dim ; k++)
			// coverty[dead_error_line]
			{dg.grids.z *= sz[k];}

			return;
		}
	}

	if (tot_block >= min_block_size)
	{return;}

	// Calculate the grids from the threads configuration

	dg.grids.x = 1;
	for (size_t k =  jx ; k < x.size() ; k++)
	{dg.grids.x *= x[k];}

	dg.grids.y = 1;
	for (size_t k = jy ; k < y.size() ; k++)
	{dg.grids.y *= y[k];}

	dg.grids.z = 1;
	for (size_t k = jz ; k < z.size() ; k++)
	{dg.grids.z *= z[k];}

	std::vector<size_t> * ptr_xyz[3];
	ptr_xyz[0] = &x;
	ptr_xyz[1] = &y;
	ptr_xyz[2] = &z;

	size_t  * jj[3];
	jj[0] = &jx;
	jj[1] = &jy;
	jj[2] = &jz;

	while (tot_block < min_block_size && (jx < x.size() || jy < y.size() || jz < z.size() ) )
	{
		size_t best_fact = std::numeric_limits<size_t>::max();
		size_t k_best = 0;

		for (size_t k = 0 ; k < dim ; k++)
		{
			if (*jj[k] < ptr_xyz[k]->size() && ptr_xyz[k]->operator[](*jj[k]) < best_fact )
			{
				best_fact = ptr_xyz[k]->operator[](*jj[k]);
				k_best = k;
			}
		}

		// The maximum block size cannot be passed
		if (tot_block * best_fact > max_block_size)
		{break;}

		if (k_best == 0)
		{
			dg.threads.x *= best_fact;
			dg.grids.x /= best_fact;
		}
		else if (k_best == 1)
		{
			dg.threads.y *= best_fact;
			dg.grids.y /= best_fact;
		}
		/* coverty[dead_error_line] */
		else if (k_best == 2)
		{
			dg.threads.z *= best_fact;
			dg.grids.z /= best_fact;
		}

		tot_block *= best_fact;

		(*jj[k_best])++;
	}
}


#endif /* OPENFPM_DATA_SRC_UTIL_COMPUTE_OPTIMAL_DEVICE_GRID_HPP_ */
