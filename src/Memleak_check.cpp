#include "config.h"
#include "Memleak_check.hpp"

#ifdef MEMLEAK_CHECK

// counter for allocation of new memory
size_t new_data;

// counter to delete memory
size_t delete_data;

// structure that store all the active pointer
std::map<byte_ptr, size_t> active_ptr;

#endif
