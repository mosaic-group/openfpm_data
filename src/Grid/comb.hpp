#ifndef COMB_HPP
#define COMB_HPP


/*! brief Position of the element of dimension d in the hyper-cube of dimension dim
 *
 * Position of the element of dimension d in the hyper-cube of dimension dim
 *
 * \tparam dim of the hyper-cube
 *
 *
 * hyper-cube define only the features of an N-dimensional hyper-cube, does not define
 * where is is located and its size, use Box for that purpose
 *
 */

template<unsigned int dim>
struct comb
{
	//! Array that store the combination
	char c[dim];

	/*! \brief check if it is a valid combination
	 *
	 * \return true if it is valid
	 *
	 */

	inline bool isValid()
	{
		for (int i = 0 ; i < dim ; i++)
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
		for (int i = 0 ; i < dim ; i++)
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
		for (int i = 0 ; i < dim ; i++)
		{
			c[i] = 0;
		}
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

		for (int i = 0 ; i < dim ; i++)
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

	inline char operator[](int i)
	{
		return c[i];
	}

	/*! \brief get the combination array pointer
	 *
	 * \return an array of char representing the combination
	 *
	 */

	inline char * getComb()
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
	inline char value(int i) const
	{
		return c[i];
	}
};

#endif
