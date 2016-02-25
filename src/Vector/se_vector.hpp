/*
 * se_vector.hpp
 *
 *  Created on: Oct 30, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_VECTOR_SE_VECTOR_HPP_
#define OPENFPM_DATA_SRC_VECTOR_SE_VECTOR_HPP_


/*! \brief It define a set of MACRO for security enhancement class 1 check
 *
 */

#define VECTOR_ERROR 2000lu

/*#define VECTOR_OVERFLOW_NATIVE(id) if (id >= v_size)\
		{\
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " overflow id: " << id << "\n";\
			size_t * err_code_pointer = (size_t *)&this->err_code;\
			*err_code_pointer = 2001;\
			ACTION_ON_ERROR(VECTOR_ERROR);\
		}*/

#endif /* OPENFPM_DATA_SRC_VECTOR_SE_VECTOR_HPP_ */
