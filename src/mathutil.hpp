#ifndef MATHUTIL_HPP
#define MATHUTIL_HPP

/*! \brief calculate the factorial
 *
 * calculate the factorial of a number
 *
 * \param f number
 * \return the factorial
 *
 */

static size_t factorial(size_t f)
{
	size_t fc = 1;

	for (size_t s = 2 ; s <= f ; s++)
	{
		fc *= s;
	}

	return fc;
}

/*! \brief C(n,k) Combination of n objects taken on group of k elements
 *
 * C(n,k) Combination of n objects taken on group of k elements, defined as
 *
 * n!/(k!(n-k)!)
 *
 */

static size_t C(size_t n, size_t k)
{
	return factorial(n)/(factorial(k)*factorial(n-k));
}

/*! \brief Round to the nearest bigger power of 2 number
 *
 * \param number
 *
 * \return nearest bigger power of 2 number
 *
 */

static size_t round_big_2(size_t n)
{
	n--;
	n |= n >> 1;   // Divide by 2^k for consecutive doublings of k up to 32,
	n |= n >> 2;   // and then or the results.
	n |= n >> 4;
	n |= n >> 8;
	n |= n >> 16;
	n++;

	return n;
}


#endif
