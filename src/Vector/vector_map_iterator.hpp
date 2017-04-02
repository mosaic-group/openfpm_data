/*
 * vector_map_iterator.hpp
 *
 *  Created on: Apr 1, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_VECTOR_VECTOR_MAP_ITERATOR_HPP_
#define OPENFPM_DATA_SRC_VECTOR_VECTOR_MAP_ITERATOR_HPP_

namespace openfpm
{
	/*! \brief Vector iterator
	 *
	 */

	class vector_key_iterator
	{
		//! Linearized end element
		size_t end;

	protected:

		//! Actual key
		size_t gk;

	public:

		/*! \brief Constructor require the size of the vector
		 *
		 * \param end point
		 * \param start starting point
		 *
		 */
		vector_key_iterator(size_t end, size_t start = 0)
		: end(end),gk(start)
		{}


		/*! \brief Get the next element
		 *
		 * Get the next element
		 *
		 * \return the next grid_key
		 *
		 */
		vector_key_iterator operator++()
		{
			//! increment the first index

			gk++;

			return *this;
		}

		/*! \brief Set the dimension
		 *
		 * \param d is the dimension (IGNORED is by default 0)
		 * \param sz set the counter to sz
		 *
		 */
		void set(int d, size_t sz)
		{
			// set the counter dim to sz

			gk = sz;
		}

		/*! \brief Check if there is the next element
		 *
		 * Check if there is the next element
		 *
		 * \return true if there is the next, false otherwise
		 *
		 */
		bool isNext() const
		{
			if (gk < end)
			{
				//! we did not reach the end of the grid

				return true;
			}

			//! we reach the end of the grid
			return false;
		}

		/*! \brief Get the actual key
		 *
		 * Get the actual key
		 *
		 * \return the actual key
		 *
		 */
		size_t get() const
		{
			return gk;
		}
	};

	/*! \brief Vector iterator
	 *
	 * Vector iterator over a predefined sequence
	 *
	 */

	template<typename lid>
	class vector_key_iterator_seq
	{
		openfpm::vector<lid> & dp;

	protected:

		//! Actual key
		size_t gk;

	public:

		/*! \brief Constructor require the sequence
		 *
		 * \param dp
		 *
		 */
		vector_key_iterator_seq(openfpm::vector<lid> & dp)
		:dp(dp),gk(0)
		{}


		/*! \brief Get the next element
		 *
		 * Get the next element
		 *
		 * \return the next grid_key
		 *
		 */
		vector_key_iterator_seq<lid> operator++()
		{
			//! increment the first index

			gk++;

			return *this;
		}

		/*! \brief Set the dimension
		 *
		 * \param d is the dimension (IGNORED is by default 0)
		 * \param sz set the counter to sz
		 *
		 */
		void set(int d, size_t sz)
		{
			// set the counter dim to sz

			gk = sz;
		}

		/*! \brief Check if there is the next element
		 *
		 * Check if there is the next element
		 *
		 * \return true if there is the next, false otherwise
		 *
		 */
		bool isNext() const
		{
			if (gk < dp.size())
			{
				//! we did not reach the end of the grid

				return true;
			}

			//! we reach the end of the grid
			return false;
		}

		/*! \brief Get the actual key
		 *
		 * Get the actual key
		 *
		 * \return the actual key
		 *
		 */
		size_t get() const
		{
			return dp.get(gk);
		}
	};
}


#endif /* OPENFPM_DATA_SRC_VECTOR_VECTOR_MAP_ITERATOR_HPP_ */
