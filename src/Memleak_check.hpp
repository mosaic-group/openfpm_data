#include <iostream>
#include <map>

#ifdef MEMLEAK_CHECK

extern size_t new_data;
extern size_t delete_data;

extern std::map<void *,size_t> active_ptr;

/*! \brief Check and remove the active pointer
 *
 * Check and remove the pointer from the active list
 *
 * \param pointer to check and remove
 *
 */
static void remove_ptr(void * ptr)
{
	// Check if the pointer exist
	std::map<void *, size_t>::iterator it = active_ptr.find(ptr);

	// if the element does not exist, print that something wrong happened and return
	if ( it == active_ptr.end() )
	{
		std::cout << "ERROR POINTER NOT FOUND " << ptr << "\n";
		return;
	}

	// erase the pointer
	active_ptr.erase(ptr);
}

/*! \brief Print all active pointer
 *
 * Print all active pointer
 *
 */
static void print_unalloc()
{
	for (std::map<void *,size_t>::iterator it = active_ptr.begin(); it != active_ptr.end(); ++it)
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
static void check_new(void * data, size_t sz)
{
	// Add a new pointer
	new_data++;
	active_ptr[data] = sz;
	std::cout << "New data: " << new_data << "   " << data << "\n";
}

/*! \brief check and delete a pointer
 *
 * check and delete a pointer from the list of active pointers
 *
 * \param pointer data
 *
 */
static void check_delete(void * data)
{
	// Delete the pointer
	delete_data++;
	remove_ptr(data);
	std::cout << "Delete data: " << delete_data << "   " << data << "\n";
}

/*! \brief check if the access is valid
 *
 * check if the access is valid
 *
 * \param ptr pointer we are going to access
 * \param size_access is the size in byte of the data we are fetching
 *
 */
static void check_valid(void * ptr, size_t size_access)
{
	// get the lower bound

	std::map<void *, size_t>::iterator l_b = active_ptr.lower_bound(ptr);

	// if where is not memory that satisfy the request
	if (l_b == active_ptr.end())
	{
		std::cout << "Error invalid pointer: " << ptr << "\n";
	}

	// check if ptr is in the range

	size_t sz = l_b->second;

	if (l_b->first + sz + size_access > ptr)
	{
		std::cout << "Error invalid pointer: " << ptr << "\n";
	}
}

#else

static void remove_ptr(void * ptr)
{
}

static void print_unalloc()
{
}

static void check_new(void * ptr)
{
}

static void check_delete(void * ptr)
{
}

static void check_valid(void * ptr)
{
}

#endif
