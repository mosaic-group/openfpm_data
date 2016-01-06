/*
 * se_grid.hpp
 *
 *  Created on: Aug 13, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_SE_GRID_HPP_
#define SRC_GRID_SE_GRID_HPP_

#include "util/se_util.hpp"

/*! \brief It define a set of MACRO for security enhancement class 1 check
 *
 */

#define GRID_ERROR 1000lu

// Macro for security enhancement
/*#define CHECK_INIT() 		if (is_mem_init == false)\
							{\
								std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " you must call SetMemory before access the grid\n";\
								size_t * err_code_pointer = (size_t *)&this->err_code;\
								*err_code_pointer = 1001;\
								ACTION_ON_ERROR(GRID_ERROR);\
							}

#define GRID_OVERFLOW(v1) for (long int i = 0 ; i < dim ; i++)\
						  {\
						  	  if (v1.get(i) >= (long int)getGrid().size(i))\
							  {\
						  		  std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << v1.get(i) << " >= " << getGrid().size(i) << "\n";\
								  size_t * err_code_pointer = (size_t *)&this->err_code;\
								  *err_code_pointer = 1002;\
						  		  ACTION_ON_ERROR(GRID_ERROR);\
							  }\
							  else if (v1.get(i) < 0)\
							  {\
						  		  std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << v1.get(i) << " is negative " << "\n";\
								  size_t * err_code_pointer = (size_t *)&this->err_code;\
								  *err_code_pointer = 1003;\
						  		  ACTION_ON_ERROR(GRID_ERROR);\
							  }\
						  }

#define GRID_OVERFLOW_EXT(g,key2) 		for (size_t i = 0 ; i < dim ; i++)\
		{\
			if (key2.get(i) >= (long int)g.g1.size(i))\
			{\
				std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << key2.get(i) << " >= " << g.g1.size(i) << "\n";\
				size_t * err_code_pointer = (size_t *)&this->err_code;\
				*err_code_pointer = 1004;\
				ACTION_ON_ERROR(GRID_ERROR);\
			}\
			else if (key2.get(i) < 0)\
			{\
				std::cerr << "Error " __FILE__ << ":" << __LINE__ <<" grid overflow " << "x=[" << i << "]=" << key2.get(i) << " is negative " << "\n";\
				size_t * err_code_pointer = (size_t *)&this->err_code;\
				*err_code_pointer = 1005;\
				ACTION_ON_ERROR(GRID_ERROR);\
			}\
		}*/

#endif /* SRC_GRID_SE_GRID_HPP_ */
