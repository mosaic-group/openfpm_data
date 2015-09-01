#include "config.h"
#include <iostream>
#include <map>

#ifndef MEMLEAK_CHECK_HPP
#define MEMLEAK_CHECK_HPP

typedef unsigned char * byte_ptr;

#ifdef SE_CLASS2

extern size_t new_data;
extern size_t delete_data;

extern std::map<byte_ptr,size_t> active_ptr;

extern long int process_v_cl;
extern long int process_to_print;

/*! \brief Check and remove the active pointer
 *
 * Check and remove the pointer from the active list
 *
 * \param pointer to check and remove
 *
 * \return true if the operation succeded, false if the pointer does not exist
 *
 */
static bool remove_ptr(void * ptr)
{
	// Check if the pointer exist
	std::map<byte_ptr, size_t>::iterator it = active_ptr.find((byte_ptr)ptr);

	// if the element does not exist, print that something wrong happened and return
	if ( it == active_ptr.end() )
	{
		std::cout << "Error pointer not found " << ptr << "\n";
		return false;
	}

	// erase the pointer
	active_ptr.erase((byte_ptr)ptr);

	return true;
}

/*! \brief Print all active pointer
 *
 * Print all active pointer
 *
 */
static void print_unalloc()
{
	for (std::map<byte_ptr,size_t>::iterator it = active_ptr.begin(); it != active_ptr.end(); ++it)
	{
		std::cout << "Unallocated memory " << it->first << "     size " << it->second << "\n";
	}
}


/*! \brief Add the new allocated active pointer
 *
 * Add the new allocated active pointer
 *
 * \param new data active pointer
 * \param sz size of the new allocated memory
 *
 */
static bool check_new(void * data, size_t sz)
{
	// Add a new pointer
	new_data++;
	active_ptr[(byte_ptr)data] = sz;

#ifdef SE_CLASS2_VERBOSE
	if (process_to_print < 0 || process_to_print == process_v_cl)
		std::cout << "New data: " << new_data << "   " << data << "\n";
#endif

	return true;
}

/*! \brief check and delete a pointer
 *
 * check and delete a pointer from the list of active pointers
 *
 * \param pointer data
 * \return true if the operation to delete succeed
 *
 */
static bool check_delete(void * data)
{
	if (data == NULL)	return true;
	// Delete the pointer
	delete_data++;
	bool result = remove_ptr(data);

#ifdef SE_CLASS2_VERBOSE
	if (process_to_print < 0 || process_to_print == process_v_cl)
		std::cout << "Delete data: " << delete_data << "   " << data << "\n";
#endif

	return result;
}

/*! \brief check if the access is valid
 *
 * check if the access is valid
 *
 * \param ptr pointer we are going to access
 * \param size_access is the size in byte of the data we are fetching
 *
 * \return true if the pointer is valid
 *
 */
static bool check_valid(void * ptr, size_t size_access)
{
	// get the lower bound

	std::map<byte_ptr, size_t>::iterator l_b = active_ptr.upper_bound((byte_ptr)ptr);

	if (active_ptr.size() == 0)
	{
		std::cout << "Error invalid pointer: " << __FILE__ << ":" << __LINE__ << "  " << ptr << "\n";
		return false;
	}

	//! subtract one
	l_b--;

	// if there is no memory that satisfy the request
	if (l_b == active_ptr.end())
	{
		if (process_to_print < 0 || process_to_print == process_v_cl)
			std::cout << "Error invalid pointer: " << __FILE__ << ":" << __LINE__ << "  " << ptr << "\n";
		return false;
	}

	// check if ptr is in the range

	size_t sz = l_b->second;

	if (((unsigned char *)l_b->first) + sz < ((unsigned char *)ptr) + size_access)
	{
		if (process_to_print < 0 || process_to_print == process_v_cl)
		{
			std::cout << "Error invalid pointer: " << __FILE__ << ":" << __LINE__ << "  "  << ptr << "\n";
			return false;
		}
	}

	return true;
}

/*! \brief In case of Parallel application it set the process that print the
 *
 * \param p_to_print is < 0 (Mean all)
 *
 */
static void set_process_to_print(long int p_to_print)
{
	process_to_print = p_to_print;
}

#else



#endif
#endif
