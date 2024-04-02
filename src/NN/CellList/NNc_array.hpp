/*
 * NNc_array.hpp
 *
 *  Created on: Feb 6, 2018
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_NNC_ARRAY_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_NNC_ARRAY_HPP_

#include "Grid/grid_sm.hpp"
#include <boost/mpl/bool.hpp>

/*! \brief Set a dimension threshold
 *
 * \param dim dimension
 *
 */
template<unsigned int dim>
struct as_array_nnc
{
	//! true only if dimension is smaller than 10
	typedef boost::mpl::bool_< dim < 10 > type;
};

/* \brief NNc_array
 *
 * \param size
 *
 */
template<unsigned int dim, unsigned int size, bool thr = as_array_nnc<dim>::type::value>
class NNc_array
{
	//! NNc_array
	long int NNc_arr[size];

	//! size of the cell array on each dimension
	grid_sm<dim,void> gs;

public:

	/*! \brief Set the size in each
	 *
	 * \param sz size og the cell grid in each dimensions
	 *
	 */
	void set_size(const size_t (& sz)[dim])
	{
		gs.setDimensions(sz);
	}

	/*! \brief return the element i
	 *
	 * \return element i
	 *
	 */
	inline const long int & operator[](size_t i) const
	{
		return NNc_arr[i];
	}

	/*! \brief return the element i
	 *
	 * \return element i
	 *
	 */
	inline long int & operator[](size_t i)
	{
		return NNc_arr[i];
	}

	/*! \brief Initialize the NNc array with full neighborhood cells indexes
	 *
	 *
	 */
	void init_full()
	{
		typedef typename generate_array<size_t,dim, Fill_zero>::result NNzero;
		typedef typename generate_array<size_t,dim, Fill_two>::result NNtwo;
		typedef typename generate_array<size_t,dim, Fill_one>::result NNone;

		// Generate the sub-grid iterator

		grid_key_dx_iterator_sub<dim> gr_sub3(gs,NNzero::data,NNtwo::data);

		// Calculate the NNc array

		size_t middle = gs.LinId(NNone::data);
		size_t i = 0;
		while (gr_sub3.isNext())
		{
			NNc_arr[i] = (long int)gs.LinId(gr_sub3.get()) - middle;

			++gr_sub3;
			i++;
		}
	}


	/*! \brief Initialize the NNc array with symmetric neighborhood cells indexes
	 *
	 *
	 */
	void init_sym()
	{
		// compile-time array {0,0,0,....}  {2,2,2,...} {1,1,1,...}

		typedef typename generate_array<size_t,dim, Fill_zero>::result NNzero;
		typedef typename generate_array<size_t,dim, Fill_two>::result NNtwo;
		typedef typename generate_array<size_t,dim, Fill_one>::result NNone;

		size_t middle = gs.LinId(NNone::data);

		// Generate the sub-grid iterator

		grid_key_dx_iterator_sub<dim> gr_sub3(gs,NNzero::data,NNtwo::data);

		// Calculate the NNc_sym array

		size_t i = 0;
		while (gr_sub3.isNext())
		{
			auto key = gr_sub3.get();

			size_t lin = gs.LinId(key);

			// Only the first half is considered
			if (lin < middle)
			{
				++gr_sub3;
				continue;
			}

			NNc_arr[i] = lin - middle;

			++gr_sub3;
			i++;
		}
	}

	/*! \brief Initialize the NNc array with symmetric local neighborhood cells indexes
	 *
	 * The neighborhood is the same as with the full iteration scheme
	 * to account for interactions with ghost particles in bottom left
	 * neighboring cells. In the symmetric case those are completely skipped
	 *
	 */
	void init_sym_local()
	{
		this->init_full();
	}

	/*! \brief return the pointer to the array
	 *
	 * \return the pointer
	 *
	 */
	const long int * getPointer() const
	{
		return NNc_arr;
	}

	/*! \brief Copy the NNc_array
	 *
	 * \param nnc NNc_array to copy
	 *
	 */
	NNc_array<dim,size,thr> & operator=(const NNc_array<dim,size,thr> & nnc)
	{
		std::copy(&nnc.NNc_arr[0],&nnc.NNc_arr[size],&NNc_arr[0]);
		gs = nnc.gs;

		return *this;
	}

	/*! \brief swap NNc_array
	 *
	 * \param nnc NNc_array to copy
	 *
	 */
	void swap(NNc_array<dim,size,thr> & nnc)
	{
		gs.swap(nnc.gs);

		long int NNc_full_tmp[openfpm::math::pow(3,dim)];

		std::copy(&nnc.NNc_arr[0],&nnc.NNc_arr[size],&NNc_full_tmp[0]);
		std::copy(&NNc_arr[0],&NNc_arr[size],&nnc.NNc_arr[0]);
		std::copy(&NNc_full_tmp[0],&NNc_full_tmp[size],&NNc_arr[0]);
	}
};

/* \brief NNc_array
 *
 * It encapsulate a 3^{dim} array containing the neighborhood cells-id
 *
 * \tparam dim dimensionality
 * \tparam size total number of neighborhood cells
 *
 */
template<unsigned int dim, unsigned int size>
class NNc_array<dim,size,false>
{
	//! NNc_array is a grid in general 3^{dim}, this object contain the information
	//! about this grid
	grid_sm<dim,void> gs;

	//! Information about the grid in the reduced space
	grid_sm<dim,void> gs_base;

	//!
	size_t sub_off;
	size_t sym_mid;

	bool full_or_sym;

public:

	/*! \brief set the size of the cell grid
	 *
	 * \param sz[dim] size of the cell grid in each dimension
	 *
	 */
	void set_size(const size_t (& sz)[dim])
	{
		typedef typename generate_array<size_t,dim, Fill_three>::result NNthree;

		gs.setDimensions(sz);
		gs_base.setDimensions(NNthree::data);

		typedef typename generate_array<size_t,dim, Fill_one>::result NNone;
		sub_off = gs.LinId(NNone::data);
		sym_mid = gs_base.LinId(NNone::data);
	}

	/*! \brief return the element i
	 *
	 * \return element i
	 *
	 */
	long int operator[](size_t i) const
	{
		if (full_or_sym == true)
		{
			grid_key_dx<dim> key = gs_base.InvLinId(i);
			return gs.LinId(key) - sub_off;
		}

		grid_key_dx<dim> key = gs_base.InvLinId(i + sym_mid);
		return gs.LinId(key) - sub_off;
	}

	void init_full()
	{
		full_or_sym = true;
	}

	void init_sym()
	{
		full_or_sym = false;
	}

	void init_sym_local()
	{
		full_or_sym = false;
	}

	/*! \brief return the pointer to the array
	 *
	 * \return the pointer
	 *
	 */
	const long int * getPointer() const
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " error dimension is too high to use this type of neighborhood" << std::endl;
		return NULL;
	}

	/*! \brief Copy the NNc_array
	 *
	 * \param nnc NNc_array to copy
	 *
	 */
	NNc_array<dim,size,false> & operator=(const NNc_array<dim,size,false> & nnc)
	{
		gs = nnc.gs;
		gs_base = nnc.gs_base;
		sub_off = nnc.sub_off;
		sym_mid = nnc.sym_mid;

		full_or_sym = nnc.full_or_sym;

		return *this;
	}

	/*! \brief swap the NNc_array
	 *
	 * \param nnc NNc_array to copy
	 *
	 */
	void swap(NNc_array<dim,size,false> & nnc)
	{
		gs.swap(nnc.gs);
		gs_base.swap(nnc.gs_base);

		size_t sub_off_tmp = sub_off;
		sub_off = nnc.sub_off;
		nnc.sub_off = sub_off_tmp;

		size_t sym_mid_tmp = sym_mid;
		sym_mid = nnc.sym_mid;
		nnc.sym_mid = sym_mid_tmp;

		bool full_or_sym_tmp = full_or_sym;
		full_or_sym = nnc.full_or_sym;
		nnc.full_or_sym = full_or_sym_tmp;
	}
};


#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_NNC_ARRAY_HPP_ */
