/*
 * memory.hpp
 *
 *  Created on: Aug 8, 2014
 *      Author: Pietro Incardona
 */

/*! \brief This class is an interface
 *
 * This class is an interface to allocate memory
 *
 */

typedef long int mem_id;

#include <stddef.h>

class memory
{
	/*! \brief allocate on device a buffer of
	 *
	 * Allocate on the device a buffer of memory
	 *
	 * \param sz the size of the buffer
	 *
	 */

	virtual bool allocate(size_t sz) = 0;

	/*! \brief resize on device a buffer
	 *
	 * On the device resize
	 *
	 * \param sz the size the resized buffer
	 *
	 */
	virtual bool resize(size_t sz) = 0;

	/*! \brief Destroy the buffer on the device
	 *
	 * destroy the buffer on the device
	 *
	 */
	virtual void destroy() = 0;


	/*! \brief Copy the memory from one device to another device
	 *
	 * copy the memory from one device to another device
	 *
	 * \param memory from where to copy
	 *
	 */
	virtual bool copy(memory & m) = 0;

	/*! \brief get the size of the buffer
	 *
	 * get the size of the buffer
	 *
	 */

	virtual size_t size() = 0;

	/*! \ brief destructor
	 *
	 * destructor
	 *
	 */
	~memory() {};
};


