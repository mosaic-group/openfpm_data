//
// Created by tommaso on 28/06/19.
//

#ifndef OPENFPM_PDATA_MATHUTILS_HPP
#define OPENFPM_PDATA_MATHUTILS_HPP

#include <cstdlib>

template<unsigned int base, unsigned int exponent>

struct IntPow
{
    constexpr static size_t value = base * IntPow<base, exponent - 1>::value;
};

template<unsigned int base>
struct IntPow<base, 0>
{
    constexpr static size_t value = 1;
};

template <unsigned int numerator, unsigned int denominator>
struct UIntDivCeil
{
    constexpr static unsigned int value = numerator / denominator + (numerator%denominator!=0);
};

#endif //OPENFPM_PDATA_MATHUTILS_HPP
