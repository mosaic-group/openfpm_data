#ifndef GRID_KEY_DX
#define GRID_KEY_DX

#include "Grid/comb.hpp"
#include "Grid/grid_key_expression.hpp"

/*! \brief In this header file the main grid accessor are defined
 *
 * Given a grid a in general a set of indexes define one element in the grid
 * In this header file several classes that encapsulate this information are
 * defined
 *
 */

/*! \brief grid_key_dx is the key to access any element in the grid
 *
 * grid_key_dx is the key to access any element in the grid
 *
 * \param dim dimensionality of the grid
 *
 */

template<unsigned int dim>
class grid_key_dx
{
public:

	// Constructor from expression
	template<typename exp1> grid_key_dx(const grid_key_dx_expression<exp1> & exp)
	{
		for (int i = 0 ; i < dim ; i++)
			this->k[i] = exp.value(i);
	}

	//! Constructor
	grid_key_dx()
	{}

	//! Constructor from buffer reference
	grid_key_dx(const size_t (&k)[dim])
	{
		for (int i = 0 ; i < dim ; i++)
			this->k[i] = k[i];
	}

	//! Construct a grid key from a list of numbers
	template<typename ...T> grid_key_dx(const size_t v,const T...t)
	{
		k[dim-1] = v;
		invert_assign(t...);
	}

	/* \brief Set to zero the key
	 *
	 *
	 */
	inline void zero()
	{
		for (int i = 0 ; i < dim ; i++)
		{
			k[i] = 0;
		}
	}

	/* \brief Set to invalid the key
	 *
	 *
	 */
	inline void invalid()
	{
		for (int i = 0 ; i < dim ; i++)
		{
			k[i] = -1;
		}
	}

	/* \brief sum an a combination to the grid_key
	 *
	 * \param comb combination (or relative movement)
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	inline grid_key_dx_sum<dim,grid_key_dx<dim>,comb<dim>> operator+(comb<dim> & cmb) const
	{
		grid_key_dx_sum<dim,grid_key_dx<dim>,comb<dim>> exp_sum(*this,cmb);
//		grid_key_dx_expression<dim,grid_key_dx_sum<dim,grid_key_dx<dim>,comb<dim>>> exp(exp_sum);

		return exp_sum;
	}

	/* \brief Check if two key are the same
	 *
	 * \param key_t key to check
	 *
	 * \return true if the two key are identical
	 *
	 */

	template<unsigned int dim_t> bool operator==(grid_key_dx<dim_t> & key_t)
		  {
		if (dim != dim_t)
		{
			return false;
		}

		// Check the two key index by index

		for (int i = 0 ; i < dim ; i++)
		{
			if (k[i] != key_t.k[i])
			{
				return false;
			}
		}

		// identical key
		return true;
		  }

	//! set the grid key from a list of numbers
	template<typename a, typename ...T>void set(a v, T...t)
	{
		k[dim-1] = v;
		invert_assign(t...);
	}

	/*! \brief Get the i index
	 *
	 * \param i index to get
	 *
	 * \return the index value
	 *
	 */
	mem_id value(size_t i) const
	{
		return k[i];
	}

	/*! \brief Get the i index
	 *
	 *
	 * \param i index to get
	 *
	 * \return the index value
	 *
	 */
	mem_id get(size_t i) const
	{
		return k[i];
	}

	/*! \brief Set the i index
	 *
	 * Set the i index
	 *
	 * \param i index to set
	 * \param id value to set
	 *
	 */
	void set_d(size_t i, mem_id id)
	{
#ifdef DEBUG

		if (i >= dim)
			std::cerr << "grid_key_dx error: " << __FILE__ << " " << __LINE__ << " try to access dimension " << i << " on a grid_key_dx of size " << dim << "\n";

#endif
		k[i] = id;
	}

	//! structure that store all the index
	mem_id k[dim];

private:

	/*! \brief Recursively invert the assignment
	 *
	 * Recursively invert the assignment at compile-time
	 *
	 */
	template<typename a, typename ...T>void invert_assign(a v,T...t)
	{
		k[sizeof...(T)] = v;
		invert_assign(t...);
	}

	template<typename a, typename ...T>void invert_assign(a v)
	{
		k[0] = v;
	}

	void invert_assign()
	{
	}

};


/*! \brief grid_key_d is the key to access any element in the grid
 *
 * grid_key_d is the key to access any element in the grid
 *
 * \param dim dimensionality of the grid
 * \param p object property to get from the element of the grid
 *
 */

template<unsigned int dim, unsigned int p>
class grid_key_d
{
public:


	template<typename a, typename ...T>grid_key_d(a v,T...t)
	{
		k[dim-1] = v;
		invert_assign(t...);
	}

	template<typename a, typename ...T>void invert_assign(a v,T...t)
	{
		k[sizeof...(T)] = v;
		invert_assign(t...);
	}

	template<typename a, typename ...T>void invert_assign(a v)
	{
		k[0] = v;
	}

	mem_id k[dim];
};


/*! \brief grid_key is the key to access any element in the grid
 *
 * grid_key is the key to access any element in the grid
 *
 * \param p dimensionality of the grid
 *
 */

template<unsigned int p>
class grid_key
{
public:

	template<unsigned int s>grid_key(grid_key<s> & key)
	{
		id = key.id;
	}

	grid_key(size_t sz)
	:id(new mem_id[sz])
	{
	}

	/*! \brief set the dimension d of the grid_key
	 *
	 * set the dimension d of the grid_key
	 *
	 * \param id_ d dimension
	 * \param val value to set
	 *
	 */

	void set_d(size_t id_, mem_id val)
	{
		id[id_] = val;
	}

	void set(mem_id val1)
	{
		id[0] = val1;
	}

	void set(mem_id val1, mem_id val2)
	{
		id[0] = val2;
		id[1] = val1;
	}

	void set(mem_id val1, mem_id val2, mem_id val3)
	{
		id[0] = val3;
		id[1] = val2;
		id[2] = val1;
	}

	void set(mem_id val1, mem_id val2, mem_id val3, mem_id val4)
	{
		id[0] = val4;
		id[1] = val3;
		id[2] = val2;
		id[3] = val1;
	}

	mem_id * getId()
	{
		return id.get();
	}

	boost::shared_array<mem_id> id;
};

#endif
