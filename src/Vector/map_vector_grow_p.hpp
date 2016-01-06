/*
 * Packer_cls.hpp
 *
 *  Created on: Jul 15, 2015
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_GROW_POLICY_HPP_
#define MAP_VECTOR_GROW_POLICY_HPP_

#define PAGE_ALLOC 1

namespace openfpm
{

	/*! \brief Grow policy define how the vector should grow every time we exceed the size
	 *
	 * In this case it return the requested size
	 *
	 */

	class grow_policy_identity
	{
	public:

		/*! \brief It say how much the vector must grow
		 *
		 * \param original size
		 * \param requested size
		 *
		 * \return how much to grow
		 *
		 */
		static size_t grow(size_t original, size_t requested)
		{
			return requested;
		}
	};

	/*! \brief Grow policy define how the vector should grow every time we exceed the size
	 *
	 * In this case it double up the size
	 *
	 */

	class grow_policy_double
	{
	public:

		/*! \brief It say how much the vector must grow
		 *
		 * \param original size
		 * \param requested size
		 *
		 * \return how much to grow
		 *
		 */
		static size_t grow(size_t original, size_t requested)
		{
			size_t grow = (original == 0)?1:original;
			while (grow < requested)	{grow *= 2;}
			return grow;
		}
	};

	//! default grow policy
	typedef grow_policy_double vector_grow_policy_default;

	/*! \brief Grow policy define how the vector should grow every time we exceed the size
	 *
	 * In this case it increase of 4096 elements
	 *
	 */

	class grow_policy_page
	{
	public:

		/*! \brief It say how much the vector must grow
		 *
		 * \param original size
		 * \param requested size
		 *
		 * \return how much to grow
		 *
		 */
		static size_t grow(size_t original, size_t requested)
		{
			return (requested / PAGE_ALLOC) * PAGE_ALLOC + PAGE_ALLOC;
		}
	};
}

#endif /* MAP_VECTOR_GROW_POLICY_HPP_ */
