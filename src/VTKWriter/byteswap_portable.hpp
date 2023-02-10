/*
 * byteswap_portable.hpp
 *
 *  Created on: Feb 20, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_VTKWRITER_BYTESWAP_PORTABLE_HPP_
#define OPENFPM_IO_SRC_VTKWRITER_BYTESWAP_PORTABLE_HPP_

#include <climits>

/*! \brief This function swap byte from little endian to big endian format
 *
 * \warning in the case of big-endian machine this function should do nothing.
 *          Unfortunately this is not the case because I never had the bad luck
 *          of getting one
 *
 * \param T value to convert
 *
 */
template <typename T>
T swap_endian_lt(T u)
{
    static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");

    union
    {
        T u;
        unsigned char u8[sizeof(T)];
    } source, dest;

    source.u = u;

    for (size_t k = 0; k < sizeof(T); k++)
    {dest.u8[k] = source.u8[sizeof(T) - k - 1];}

    return dest.u;
}

#endif /* OPENFPM_IO_SRC_VTKWRITER_BYTESWAP_PORTABLE_HPP_ */
