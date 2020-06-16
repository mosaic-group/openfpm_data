/*
 * grid_key_dx_iterator_sub_p.hpp
 *
 *  Created on: Dec 15, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SUB_BC_HPP_
#define OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SUB_BC_HPP_

#include "grid_key_dx_iterator_sub.hpp"

/*! \brief The same as grid_key_dx_iterator_sub_p but with periodic boundary
 *
 * In this case the boundaries are periodic
 *
 */
template<unsigned int dim, typename stencil=no_stencil, typename linearizer = grid_sm<dim,void>, typename warn=print_warning_on_adjustment<dim,linearizer>>
class grid_key_dx_iterator_sub_bc : public grid_key_dx_iterator_sub<dim,stencil,linearizer,warn>
{
	//! Boundary conditions
	size_t bc[dim];

	//! actual iterator box
	size_t act;

	//! Here we have all the boxes that this iterator produce
	std::vector<Box<dim,size_t>> boxes;

	/*! \brief Check if the position is valid with the actual boundary conditions
	 *
	 * \param key key to check
	 *
	 */
	bool inline check_invalid(const grid_key_dx<dim> & key, const size_t (& bc)[dim])
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == NON_PERIODIC && key.get(i) != 1)
			{
				return true;
			}
		}

		return false;
	}

public:


	grid_key_dx_iterator_sub_bc(const grid_key_dx_iterator_sub_bc & tmp)
	{
		this->operator=(tmp);

	}

	grid_key_dx_iterator_sub_bc(grid_key_dx_iterator_sub_bc && tmp)
	{
		this->operator=(tmp);

	}

	grid_key_dx_iterator_sub_bc & operator=(const grid_key_dx_iterator_sub_bc & tmp)
	{
		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = tmp.bc[i];

		act = tmp.act;
		boxes = tmp.boxes;

		return *this;
	}

	grid_key_dx_iterator_sub_bc & operator=(grid_key_dx_iterator_sub_bc && tmp)
	{
		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = tmp.bc[i];

		act = tmp.act;
		boxes.swap(tmp.boxes);

		return *this;
	}

	/*! \brief Constructor
	 *
	 *
	 */
	grid_key_dx_iterator_sub_bc()
	:act(0)
	{

	}


	/*! \brief Constructor
	 *
	 * \param g Grid information
	 * \param start starting point
	 * \param stop stop point
	 * \param bc boundary conditions
	 *
	 */
	grid_key_dx_iterator_sub_bc(const linearizer & g,
								const grid_key_dx<dim> & start,
								const grid_key_dx<dim> & stop,
								const size_t (& bc)[dim])
	:act(0)
	{
		Initialize(g,start,stop,bc);
	}

	/*! \brief Initialize the iterator
	 *
	 * \param g Grid information
	 * \param start starting point
	 * \param stop stop point
	 * \param bc boundary conditions
	 *
	 */
	void Initialize(const linearizer & g,
										 const grid_key_dx<dim> & start ,
										 const grid_key_dx<dim> & stop,
										 const size_t (& bc)[dim])
	{
		// copy the boundary conditions

		for (size_t i = 0 ; i < dim ; i++)
		{this->bc[i] = bc[i];}

		// compile-time array {3,3,3,....}

		typedef typename generate_array<size_t,dim, Fill_three>::result NNthree;

		// In case of high dimension grid_key_dx_iterator_sub can become exponentially
		// expensive, on the other hand we are expecting that some of the dimensions
		// are cropped by non periodicity, so we adjust the iteration

		size_t end_p[dim];
		size_t start_p[dim];

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (this->bc[i] == NON_PERIODIC && g.size(i) == 1)
			{start_p[i] = 1; end_p[i] = 1;}
			else
			{start_p[i] = 0; end_p[i] = 2;}
		}

		// Generate the sub-grid iterator

		grid_sm<dim,void> nn(NNthree::data);
		grid_key_dx_iterator_sub<dim,stencil> it(nn,start_p,end_p);

		// Box base
		Box<dim,long int> base_b(start,stop);

		// intersect with all the boxes
		while (it.isNext())
		{
			auto key = it.get();

			if (check_invalid(key,bc) == true)
			{
				++it;
				continue;
			}

			bool intersect;

			// intersection box
			Box<dim,long int> b_int;
			Box<dim,long int> b_out;

			for (size_t i = 0 ; i < dim ; i++)
			{
				b_int.setLow(i,(key.get(i)-1)*g.getSize()[i]);
				b_int.setHigh(i,key.get(i)*g.getSize()[i]-1);
			}

			intersect = base_b.Intersect(b_int,b_out);

			// Bring to 0 and size[i]
			for (size_t i = 0 ; i < dim ; i++)
			{
				if (bc[i] == PERIODIC)
				{
					b_out.setLow(i,openfpm::math::positive_modulo(b_out.getLow(i),g.size(i)));
					b_out.setHigh(i,openfpm::math::positive_modulo(b_out.getHigh(i),g.size(i)));
				}
			}

			// if intersect add in the box list
			if (intersect == true)
			{boxes.push_back(b_out);}

			++it;
		}

		// initialize the first iterator
		if (boxes.size() > 0)
		{grid_key_dx_iterator_sub<dim,stencil,linearizer,warn>::reinitialize(grid_key_dx_iterator_sub<dim>(g,boxes[0].getKP1(),boxes[0].getKP2()));}
	}

	/*! \brief Get the next element
	 *
	 * Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	grid_key_dx_iterator_sub_bc<dim,stencil,linearizer,warn> & operator++()
	{
		grid_key_dx_iterator_sub<dim,stencil,linearizer,warn>::operator++();
		if (grid_key_dx_iterator_sub<dim,stencil,linearizer,warn>::isNext() == true)
		{
			return *this;
		}
		else
		{
			act++;
			if (act < boxes.size())
			{grid_key_dx_iterator_sub<dim,stencil,linearizer,warn>::reinitialize(grid_key_dx_iterator_sub<dim>(this->getGridInfo(),boxes[act].getKP1(),boxes[act].getKP2()));}
		}

		return *this;
	}

	/*! \brief Check if there is the next element
	 *
	 * Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */

	inline bool isNext()
	{
		return act < boxes.size();
	}

	/*! \brief Return the actual grid key iterator
	 *
	 * \return the actual key
	 *
	 */
	inline const grid_key_dx<dim> get() const
	{
		return grid_key_dx_iterator_sub<dim,stencil,linearizer,warn>::get();
	}

	/*! \brief Reset the iterator (it restart from the beginning)
	 *
	 */
	inline void reset()
	{
		act = 0;
		// initialize the first iterator
		reinitialize(boxes.get(0).getKP1,boxes.get(0).getKP2);
	}
};


#endif /* OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SUB_BC_HPP_ */
