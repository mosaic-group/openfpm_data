/*
 * vector_impl.hpp
 *
 *  Created on: Mar 8, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_IMPL_HPP_
#define VECTOR_IMPL_HPP_

#include "common.hpp"

#define STD_VECTOR 1
#define OPENFPM_NATIVE 2

/*! \brief It analyze the type given and it select correctly the implementation
 *  for vector
 *
 * \tparam type to analyze
 *
 * [Example]
 *
 * vect_isel<T>::value
 *
 * will return 1 for std base implementation
 * will return 2 for openfpm native implementation
 *
 * Basically the openfpm native implementation require that T
 * has some specific structure, this class check for it, if T
 * does not have this structure it fall to the case 1
 *
 */

namespace openfpm
{
	template<typename T>
	struct vect_isel
	{
		enum
		{
			value = is_typedef_and_data_same<has_typedef_type<T>::value && has_data<T>::value,T>::value + 1
		};
	};
}


#endif /* VECTOR_IMPL_HPP_ */
