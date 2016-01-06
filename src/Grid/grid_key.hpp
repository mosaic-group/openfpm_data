#ifndef GRID_KEY_DX
#define GRID_KEY_DX

#include "Grid/comb.hpp"
#include "Grid/grid_key_expression.hpp"
#include "Space/Shape/Point.hpp"

/*! \brief grid_key_dx is the key to access any element in the grid
 *
 * Given a grid in general a set of indexes define one element in the grid
 * For example 2 indexes identify one element on a two dimensional grid,
 * 3 indexes on a 3 dimensional grid ...
 *
 * \tparam dim dimensionality of the grid
 *
 */
template<unsigned int dim>
class grid_key_dx
{
public:

	// Constructor from expression
	template<typename exp1> inline grid_key_dx(const grid_key_dx_expression<dim,exp1> & exp)
	{
		for (size_t i = 0 ; i < dim ; i++)
			this->k[i] = exp.value(i);
	}

	//! Constructor
	inline grid_key_dx()
	{}

	//! Constructor from an other key
	inline grid_key_dx(const grid_key_dx<dim> & key)
	:grid_key_dx(key.k)
	{
	}

	//! Constructor from buffer reference
	inline grid_key_dx(const size_t (&k)[dim])
	{
		for (size_t i = 0 ; i < dim ; i++)
			this->k[i] = k[i];
	}

	//! Constructor from buffer reference
	inline grid_key_dx(const long int (&k)[dim])
	{
		for (size_t i = 0 ; i < dim ; i++)
			this->k[i] = k[i];
	}

	//! Construct a grid key from a list of numbers
	template<typename ...T> inline grid_key_dx(const size_t v,const T...t)
	{
#ifdef DEBUG
		if (sizeof...(t) != dim -1)
			std::cerr << "Error grid_key: " << __FILE__ << " " << __LINE__ << "creating a key of dimension " << dim << " require " << dim << " numbers " << sizeof...(t) + 1 << " provided" << "\n";
#endif
		k[dim-1] = v;
		invert_assign(t...);
	}

	/* \brief Set to zero the key
	 *
	 */
	inline void zero()
	{
		for (size_t i = 0 ; i < dim ; i++)
			k[i] = 0;
	}

	/* \brief Set to one the key
	 *
	 */
	inline void one()
	{
		for (size_t i = 0 ; i < dim ; i++)
			k[i] = 1;
	}

	/* \brief Set to invalid the key
	 *
	 */
	inline void invalid()
	{
		for (size_t i = 0 ; i < dim ; i++)
			k[i] = -1;
	}

	/* \brief sum a grid_key
	 *
	 * \param comb combination (or relative movement)
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	inline grid_key_dx<dim> & operator+=(const grid_key_dx<dim> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
			k[i] += p.k[i];

		return *this;
	}

	/* \brief sum a grid_key
	 *
	 * \param comb combination (or relative movement)
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	inline grid_key_dx<dim> & operator-=(const grid_key_dx<dim> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
			k[i] -= p.k[i];

		return *this;
	}

	/* \brief sum a grid_key to the grid_key
	 *
	 * \param p grid_key to sum
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	inline grid_key_dx_sum<dim,grid_key_dx<dim>,grid_key_dx<dim>> operator+(const grid_key_dx<dim> & p) const
	{
		grid_key_dx_sum<dim,grid_key_dx<dim>,grid_key_dx<dim>> exp_sum(*this,p);

		return exp_sum;
	}

	/* \brief sum an a combination to the grid_key
	 *
	 * \param comb combination (or relative movement)
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	inline grid_key_dx_sum<dim,grid_key_dx<dim>,Point<dim,long int>> operator+(const Point<dim,long int> & p) const
	{
		grid_key_dx_sum<dim,grid_key_dx<dim>,Point<dim,long int>> exp_sum(*this,p);

		return exp_sum;
	}

	/* \brief sum an a combination to the grid_key
	 *
	 * \param comb combination (or relative movement)
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	inline grid_key_dx_sum<dim,grid_key_dx<dim>,comb<dim>> operator+(const comb<dim> & cmb) const
	{
		grid_key_dx_sum<dim,grid_key_dx<dim>,comb<dim>> exp_sum(*this,cmb);

		return exp_sum;
	}

	/* \brief sum an a combination to the grid_key
	 *
	 * \param comb combination (or relative movement)
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	inline grid_key_dx_sub<dim,grid_key_dx<dim>,grid_key_dx<dim>> operator-(const grid_key_dx<dim> & cmb) const
	{
		grid_key_dx_sub<dim,grid_key_dx<dim>,grid_key_dx<dim>> exp_sum(*this,cmb);

		return exp_sum;
	}

	/* \brief sum this key to another grid expression
	 *
	 * \param cmb expression
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	template <typename T> inline grid_key_dx_sub<dim,grid_key_dx<dim>,grid_key_dx_expression<dim,T>> operator-(const grid_key_dx_expression<dim,T> & cmb) const
	{
		grid_key_dx_sub<dim,grid_key_dx<dim>,grid_key_dx_expression<dim,T>> exp_sum(*this,cmb);

		return exp_sum;
	}

	/* \brief Check if two key are the same
	 *
	 * \param key_t key to check
	 *
	 * \return true if the two key are equal
	 *
	 */
	template<unsigned int dim_t> bool operator==(const grid_key_dx<dim_t> & key_t)
	{
		if (dim != dim_t)
		{
			return false;
		}

		// Check the two key index by index

		for (size_t i = 0 ; i < dim ; i++)
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
#ifdef DEBUG
		if (sizeof...(t) != dim -1)
			std::cerr << "Error grid_key: " << __FILE__ << " " << __LINE__ << "setting a key of dimension " << dim << " require " << dim << " numbers " << sizeof...(t) + 1 << " provided" << "\n";
#endif
		k[dim-1] = v;
		invert_assign(t...);
	}

	/*! \brief Convert to a point the grid_key_dx
	 *
	 * \see toPoint
	 *
	 */
	Point<dim,long int> toPointS() const
	{
		Point<dim,long int> p;

		for (size_t i = 0; i < dim ; i++)
		{
			p.get(i) = get(i);
		}

		return p;
	}

	/*! \brief convert the information into a string
	 *
	 * \return a string
	 *
	 */
	std::string to_string()
	{
		return this->toPointS().toString();
	}

	/*! \brief Convert to a point the grid_key_dx
	 *
	 * \see toPointS
	 *
	 */
	Point<dim,size_t> toPoint() const
	{
		Point<dim,size_t> p;

		for (size_t i = 0; i < dim ; i++)
		{
			p.get(i) = get(i);
		}

		return p;
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


#endif
