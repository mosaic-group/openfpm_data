#ifndef MATHUTIL_HPP
#define MATHUTIL_HPP

#ifdef HAVE_LIBQUADMATH
#include <boost/multiprecision/float128.hpp>
#endif

namespace openfpm
{
	namespace math
	{
		/*! \brief calculate the factorial
		 *
		 * calculate the factorial of a number
		 *
		 * # Example
		 * \snippet mathutil_unit_test.hpp factorial usage
		 *
		 * \param f number
		 * \return the factorial
		 *
		 */
		static inline size_t factorial(size_t f)
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
		 * # Example
		 * \snippet mathutil_unit_test.hpp Combination usage
		 *
		 * \param n
		 * \param k
		 *
		 */
		static inline size_t C(size_t n, size_t k)
		{
			return factorial(n)/(factorial(k)*factorial(n-k));
		}

		/*! \brief Round to the nearest bigger power of 2 number
		 *
		 * \param n number
		 * \return nearest bigger power of 2 number
		 *
		 * # Example
		 * \snippet mathutil_unit_test.hpp round to big pow
		 *
		 *
		 */
		static inline size_t round_big_2(size_t n)
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


		/*! \brief It calculate at compile-time and runtime the power with recursion
		 *
		 * # Example
		 * \snippet mathutil_unit_test.hpp pow
		 *
		 * \tparam type of the pow expression
		 *
		 * \param base
		 * \param exponent
		 *
		 */

		template<class T>
		inline constexpr size_t pow(const T base, unsigned const exponent)
		{
			// (parentheses not required in next line)
			return (exponent == 0) ? 1 : (base * pow(base, exponent-1));
		}

		/* \brief Return the positive modulo of a number
		 *
		 * # Example
		 * \snippet mathutil_unit_test.hpp positive modulo
		 *
		 * \param i number
		 * \param n modulo
		 *
		 */
		static inline long int positive_modulo(long int i, long int n)
		{
		    return (i % n + n) % n;
		}

		/*! \brief Bound the position to be inside p2 and p1
		 *
		 * Given pos = 10.9, p2 = 1.0 and p1 = 0.1 and l = p2 - p1 = 1.0,
		 * it return 10.9 - (integer)( (10.9 - 0.1) / 1.0) * 1.0
		 *
		 * # Example
		 * \snippet mathutil_unit_test.hpp periodic
		 *
		 * \param pos
		 * \param l
		 * \param b
		 *
		 * \return the bound number
		 *
		 */
		template<typename T> static inline T periodic(const T & pos, const T & p2, const T & p1)
		{
			T pos_tmp;

			pos_tmp = pos - (p2 - p1) * (long int)( (pos -p1) / (p2 - p1));
			pos_tmp += (pos < p1)?(p2 - p1):0;

			return pos_tmp;
		}

		/*! \brief Bound the position to be inside p2 and p1
		 *
		 * It is like periodic but faster
		 *
		 * \warning pos should not overshoot the domain more than one time, for example
		 *          if p2 = 1.1 and p1 = 0.1, pos can be between 2.1 and -0.9
		 *
		 * # Example
		 * \snippet mathutil_unit_test.hpp periodic_l
		 *
		 * \param pos
		 * \param l
		 * \param b
		 *
		 * \return the bound number
		 *
		 */
		template<typename T> static inline T periodic_l(const T & pos, const T & p2, const T & p1)
		{
			T pos_tmp = pos;

			if (pos >= p2)
			{
				pos_tmp = p1 + (pos - p2);
			}
			else if (pos < p1)
			{
				pos_tmp = p2 - (p1 - pos);

				// This is a round off error fix
				// if the shift bring exactly on p2 p2 we back the particle to p1
				if (pos_tmp == p2)
					pos_tmp = p1;
			}

			return pos_tmp;
		}

		/*! \brief floor math function
		 *
		 *
		 */
		inline long int size_t_floor(double x)
		{
		  size_t i = (long int)x; /* truncate */
		  return i - ( i > x ); /* convert trunc to floor */
		}

		/*! \brief floor math function
		 *
		 *
		 */
		inline long int size_t_floor(float x)
		{
		  size_t i = (long int)x; /* truncate */
		  return i - ( i > x ); /* convert trunc to floor */
		}

		/*! \brief floor math function
		 *
		 *
		 */
		__device__ __host__ inline int uint_floor(double x)
		{
		 unsigned  int i = (int)x; /* truncate */
		  return i - ( i > x ); /* convert trunc to floor */
		}

		/*! \brief floor math function
		 *
		 *
		 */
		__device__ __host__ inline int uint_floor(float x)
		{
		  unsigned int i = (int)x; /* truncate */
		  return i - ( i > x ); /* convert trunc to floor */
		}

		const int tab64[64] = {
		    63,  0, 58,  1, 59, 47, 53,  2,
		    60, 39, 48, 27, 54, 33, 42,  3,
		    61, 51, 37, 40, 49, 18, 28, 20,
		    55, 30, 34, 11, 43, 14, 22,  4,
		    62, 57, 46, 52, 38, 26, 32, 41,
		    50, 36, 17, 19, 29, 10, 13, 21,
		    56, 45, 25, 31, 35, 16,  9, 12,
		    44, 24, 15,  8, 23,  7,  6,  5};

		/*! \brief Calculate the logarith 2 of a 64 bit integer
		 *
		 * \param value
		 *
		 */
		inline int log2_64 (uint64_t value)
		{
		    value |= value >> 1;
		    value |= value >> 2;
		    value |= value >> 4;
		    value |= value >> 8;
		    value |= value >> 16;
		    value |= value >> 32;
		    return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
		}


#ifdef HAVE_LIBQUADMATH

		/*! \brief floor math function
		 *
		 *
		 */
		inline long int size_t_floor(boost::multiprecision::float128 x)
		{
		  size_t i = (long int)x; /* truncate */
		  return i - ( i > x ); /* convert trunc to floor */
		}

#endif
	}
}

#endif
