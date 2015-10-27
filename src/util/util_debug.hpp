/*
 * util.hpp
 *
 *  Created on: Nov 20, 2014
 *      Author: Pietro Incardona
 */

#ifndef UTIL_DEBUG_HPP
#define UTIL_DEBUG_HPP

/*! \brief check that two array match
 *
 * check that two array match
 *
 * \param ptr1 Point to array 1
 * \param ptr2 Pointer to array 2
 * \param sz Size of the array
 *
 */

template<typename T> void boost_check_array(const T * ptr1, const T * ptr2, size_t sz)
{
	// Check the array
	for (size_t i = 0 ; i < sz ; i++)
	{
		BOOST_REQUIRE_EQUAL(ptr1[i],ptr2[i]);
	}
}

#include <typeinfo>

#include <cxxabi.h>

/*! \brief typeid().name() return a mangled name this function demangle it creating a more readable type string
 *
 * \param mangle string
 *
 * \return demangled string
 */
static std::string demangle(const char* name) {

    int status = -4; // some arbitrary value to eliminate the compiler warning

    // enable c++11 by passing the flag -std=c++11 to g++
    std::unique_ptr<char, void(*)(void*)> res {
        abi::__cxa_demangle(name, NULL, NULL, &status),
        std::free
    };

    return (status==0) ? res.get() : name ;
}

#endif
