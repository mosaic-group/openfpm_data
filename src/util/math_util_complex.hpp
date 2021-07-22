/*
 * math_util_complex.hpp
 *
 *  Created on: Oct 1, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_UTIL_MATH_UTIL_COMPLEX_HPP_
#define OPENFPM_DATA_SRC_UTIL_MATH_UTIL_COMPLEX_HPP_

//! extern vector for getFactorization
extern std::vector<int> sieve_spf;

namespace openfpm
{
	namespace math
	{

		#define SIEVE_MAXN   4096

		/*! \brief Precompute SPF (Smallest Prime Factor) for every
		 * number till MAXN.
		 * Time Complexity : O(nloglogn)
		 *
		 */
		void inline init_getFactorization()
		{
			sieve_spf.resize(SIEVE_MAXN);

			sieve_spf[1] = 1;
			for (int i = 2; i < SIEVE_MAXN; i++)
			{
				// marking smallest prime factor for every
				// number to be itself.
				sieve_spf[i] = i;
			}

			// separately marking spf for every even
			// number as 2
			for (int i = 4; i < SIEVE_MAXN; i += 2)
			{sieve_spf[i] = 2;}

			for (int i = 3; i*i < SIEVE_MAXN; i++)
			{
				// checking if i is prime
				if (sieve_spf[i] == i)
				{
					// marking SPF for all numbers divisible by i
					for (int j=i*i; j < SIEVE_MAXN; j+=i)
					{
						// marking spf[j] if it is not
						// previously marked
						if (sieve_spf[j] == j)
						{sieve_spf[j] = i;}
					}
				}
			}
		}

		/*! \brief A O(log n) function returning prime factorization
		 * by dividing by smallest prime factor at every step
		 *
		 * \param x number to factor
		 * \param factorization output
		 *
		 */
		inline void getFactorization(int x, std::vector<size_t> & ret)
		{
			ret.clear();
			while (x != 1)
			{
				ret.push_back(sieve_spf[x]);
				x = x / sieve_spf[x];
			}
		}

	}
}


#endif /* OPENFPM_DATA_SRC_UTIL_MATH_UTIL_COMPLEX_HPP_ */
