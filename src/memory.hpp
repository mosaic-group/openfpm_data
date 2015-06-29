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

#ifndef MEMORY_HPP_
#define MEMORY_HPP_

typedef long int mem_id;

#include "config.h"
#include <stddef.h>

class memory
{
	public:

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
	 * \param m from where to copy
	 *
	 */
	virtual bool copy(memory & m) = 0;

	/*! \brief get the size of the buffer
	 *
	 */

	virtual size_t size() = 0;

	/*! \brief return a data pointer
	 *
	 * return readable pointer with the data stored
	 *
	 */

	virtual void * getPointer() = 0;

	/*! \brief destructor
	 *
	 * destructor
	 *
	 */
	virtual ~memory() {};

	/*! \brief Increment the internal reference counter
	 *
	 *
	 */
	virtual void incRef() = 0;

	/*! \brief Decrement the internal reference counter
	 *
	 *
	 */
	virtual void decRef() = 0;

	/*! \brief Return the actual reference counter
	 *
	 * \return the reference counter
	 *
	 */
	virtual long int ref() = 0;

	/*! \brief Return if the actual memory that is going to be allocated is already initialized
	 *
	 * \return true if already initialized
	 *
	 */
	virtual bool isInitialized() = 0;


};

#endif
