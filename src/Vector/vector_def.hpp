/*
 * vector_def.hpp
 *
 *  Created on: Jun 29, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_VECTOR_VECTOR_DEF_HPP_
#define OPENFPM_DATA_SRC_VECTOR_VECTOR_DEF_HPP_


namespace openfpm
{
	template<typename T, typename Memory=HeapMemory, template<typename> class layout_base=memory_traits_lin , typename grow_p=grow_policy_double, unsigned int impl=vect_isel<T>::value> class vector;
}


#endif /* OPENFPM_DATA_SRC_VECTOR_VECTOR_DEF_HPP_ */
