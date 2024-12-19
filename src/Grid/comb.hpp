#ifndef COMB_HPP
#define COMB_HPP

#define COMB_ERROR 1001lu

#include "util/se_util.hpp"

/*! \brief Position of the element of dimension d in the hyper-cube of dimension dim
 *
 * These objects are returned by the Hyper-cube static functions
 * The number of non-zero d define the dimensionality of the object ( +1 or -1 its position in the hypercube)
 *
 * ## Example
 *
  \verbatim

                           (0,1)
				(-1,1)  +---------+ (1,1)
					    |         |
					    |         |
		          (-1,0)|  (0,0)  | (1,0)
					    |         |
					    |         |
			    (-1,-1) +---------+ (1,-1)
			               (0,-1)
  \endverbatim
 *
 *
 * \tparam dim of the hyper-cube
 *
 *
 */
template<unsigned int dim>
struct comb
{
	//! Array that store the combination
	signed char c[dim];

	/*! \brief check if it is a valid combination
	 *
	 * \return true if it is valid
	 *
	 */

	inline bool isValid()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (c[i] != 1 && c[i] != -1 && c[i] != 0)
			{
				return false;
			}
		}

		return true;
	}

	/*! \brief Check if the combination is a sub-element
	 *
	 * \param cmb combination to check if it is sub-element
	 *
	 * \return true if cmb is a sub-element of this combination
	 *
	 */

	inline bool isSub(comb<dim> cmb)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (c[i] != 0 && c[i] != cmb.c[i])
			{
				return false;
			}
		}

		return true;
	}

	/*! \brief Set all the elements to zero
	 *
	 */

	inline void zero()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			c[i] = 0;
		}
	}

	/*! \brief Set all the elements to -1
	 *
	 */

	inline void mone()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			c[i] = -1;
		}
	}

	/*! \brief Bitwise operator &
	 *
	 * \return Result combination
	 *
	 */
	inline comb<dim> operator&(signed char c_)
	{
		comb<dim> ret;

		for (size_t i = 0 ; i < dim ; i++)
		{
			ret.c[i] = c[i] & c_;
		}

		return ret;
	}

	/*! \brief Flip the sign of the combination
	 *
	 *
	 */
	inline void sign_flip()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			c[i] = -c[i];
		}
	}

	/*! \brief Subtract the combinations and return the result
	 *
	 * \return Result combination
	 *
	 */
	inline comb<dim> operator-(const comb<dim> & t)
	{
		comb<dim> ret;

		for (size_t i = 0 ; i < dim ; i++)
		{
			ret.c[i] = c[i] - t.c[i];
		}

		return ret;
	}

	/*! \brief flip the coefficent of the combination to be < 0
	 *
	 *
	 * * 0 --> 0
	 * * 1 --> -1
	 * * 2 --> -2
	 * * -1 --> -1
	 *
	 * \return flipped combination
	 *
	 */
	inline comb<dim> flip()
	{
		comb<dim> ret;

		for (size_t i = 0 ; i < dim ; i++)
		{
			ret.c[i] = -abs(c[i]);
		}

		return ret;
	}

	/*! \brief Subtract the combinations and return the result
	 *
	 * \return Result combination
	 *
	 */
	inline comb<dim> operator-()
	{
		comb<dim> ret;

		for (size_t i = 0 ; i < dim ; i++)
		{
			ret.c[i] = -c[i];
		}

		return ret;
	}

	/*! \brief sum the combinations and return the result
	 *
	 * \return Result combination
	 *
	 */
	inline comb<dim> operator+(const comb<dim> & t)
	{
		comb<dim> ret;

		for (size_t i = 0 ; i < dim ; i++)
		{
			ret.c[i] = c[i] + t.c[i];
		}

		return ret;
	}

	/*! \brief Compare two combination
	 *
	 * check if they match
	 *
	 * \param t combination to check
	 * \return true if the two combination match, false otherwise
	 *
	 */
	inline bool operator!=(const comb<dim> & t) const
	{
		// Check if the two combination match

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (c[i] != t.c[i])
				return true;
		}

		// They match

		return false;
	}

	/*! \brief Compare two combination
	 *
	 * check if they match
	 *
	 * \param t combination to check
	 * \return true if the two combination match, false otherwise
	 *
	 */
	inline bool operator==(const comb<dim> & t) const
	{
		return !this->operator!=(t);
	}

	/*! \brief Get the i combination coefficient
	 *
	 * \param i coefficent
	 * \return the coefficent of the combination at position i
	 *
	 */

	inline signed char operator[](int i) const
	{
		return c[i];
	}

	/*! \brief get the combination array pointer
	 *
	 * \return an array of char representing the combination
	 *
	 */

	inline signed char * getComb()
	{
		return c;
	}

	/*! \brief get the combination array pointer
	 *
	 * \return an array of char representing the combination
	 *
	 */

	inline const signed char * getComb() const
	{
		return c;
	}

	/*! \brief get the index i of the combination
	 *
	 * NOTE: used on expression template
	 *
	 * \param i index
	 *
	 * \return value of the i index
	 *
	 */
	inline signed char value(int i) const
	{
		return c[i];
	}


	/* \brief It return the number of zero in the combination
	 *
	 * \return number of zero
	 *
	 */
	inline int n_zero() const
	{
		int zero = 0;

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (c[i] == 0) zero++;
		}

		return zero;
	}

	//! Default constructor
	comb()
	{}

	/*! \brief Constructor from a list of numbers
	 *
	 * \param c list of numbers
	 *
	 */
	comb(std::initializer_list<signed char> c)
	{
		size_t i = 0;
	    for(signed char x : c)
	    {this->c[c.size() - i - 1] = x;i++;}
	}

	/*! \brief Check if any alement in the combination is <= 0
	 *
	 *    +----#----+
	 *    |         |
	 *    |         |
	 *    #    *    #
	 *    |         |
	 *    |         |
	 *    +----#----+
	 *
	 *  Is negative return true when the combination indicate * the down-left vertex + and down and left edge #
	 *
	 */
	bool isNegative()
	{
		for (size_t i = 0; i < dim ; i++)
		{
			if (c[i] > 0)
				return false;
		}

		return true;
	}

	/*! \brief Convert a combination into string
	 *
	 * \return a string reppresentation of the combination
	 *
	 */
	std::string to_string() const
	{
		size_t i;

		std::stringstream str;
		str << "(";

		// convert to string
		for (i = 0 ; i < dim - 1 ; i++)
			str << std::to_string(c[i]) << ",";

		str << std::to_string(c[i]) << ")";

		return str.str();
	}

	/*! \brief Linearization
	 *
	 * From the combination produce a number
	 *
	 * \return the number
	 *
	 */
	size_t inline lin() const
	{
		size_t ret = 0;
		size_t accu = 1;

		for (size_t i = 0 ; i < dim ; i++)
		{
			ret += (c[i] + 1) * accu;
			accu *= 3;
		}

		return ret;
	}
};

/*! brief specialization of comb in case of dim 0
 *
 *
 */

template<>
struct comb<0>
{
	//! FIX
	signed char c[0];

	/*! \brief check if it is a valid combination
	 *
	 * \return true if it is valid
	 *
	 */

	inline bool isValid()
	{
		return true;
	}

	/*! \brief Check if the combination is a sub-element
	 *
	 * \param cmb combination to check if it is sub-element
	 *
	 * \return true if cmb is a sub-element of this combination
	 *
	 */

	inline bool isSub(comb<0> cmb)
	{
		return true;
	}

	/*! \brief Set all the elements to zero
	 *
	 */

	inline void zero()
	{
	}


	/*! \brief Compare two combination
	 *
	 * check if they match
	 *
	 * \param t combination to check
	 * \return true if the two combination match, false otherwise
	 *
	 */
	inline bool operator!=(const comb<0> & t) const
	{
		// They match

		return true;
	}

	/*! \brief Compare two combination
	 *
	 * check if they match
	 *
	 * \param t combination to check
	 * \return true if the two combination match, false otherwise
	 *
	 */
	inline bool operator==(const comb<0> & t) const
	{
		return true;
	}

	/*! \brief Get the i combination coefficient
	 *
	 * \param i coefficent
	 * \return the coefficent of the combination at position i
	 *
	 */

	inline signed char operator[](int i)
	{
		return 0;
	}

	/*! \brief get the combination array pointer
	 *
	 * \return an array of char representing the combination
	 *
	 */

	inline signed char * getComb()
	{
		return c;
	}

	/*! \brief get the index i of the combination
	 *
	 * NOTE: used on expression template
	 *
	 * \param i index
	 *
	 * \return value of the i index
	 *
	 */
	inline signed char value(int i) const
	{
		return c[i];
	}


	/* \brief It return the number of zero in the combination
	 *
	 * \return number of zero
	 *
	 */
	inline int n_zero()
	{
		return 0;
	}

	/* \brief produce an unique number from the combination
	 *
	 *
	 *
	 */
	inline size_t lin()
	{
		return 0;
	}
};

// create a specialization of std::vector<comb<0>>

namespace std
{
	/*! \brief Stub vector specialization
	 *
	 *  On compiler previous 4.9.2 the code compile, and provide trivial functionality,
     *(complex operations produce crash at runtime)
	 *
	 * On 4.9.2 does not even compile
     * So we create a particular case that pass compilation give basic functionality, and
     * crash in case of complex usage
     *
    */

	template<>
	class vector<comb<0>>
	{
		//! Pointer to nothing
		comb<0> * ptr;

	public:

		/*! \brief Return 0
		 *
		 * \return 0
		 *
		 */
		size_t size()
		{
			return 0;
		}

		/*! \brief Do nothing
		 *
		 * \param i ignored
		 *
		 * \return nothing
		 *
		 */
		comb<0> & operator[](size_t i)
		{
			return ptr[i];
		}

		/*! \brief Do nothing
		 *
		 * \param obj ignored
		 *
		 */
		void push_back(comb<0> & obj)
		{
		}
	};
}

#endif
