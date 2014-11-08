#include "Memleak_check.hpp"

// counter for allocation of new memory
size_t new_data;

// counter to delete memory
size_t delete_data;

// structure that store all the active pointer
std::map<void *, size_t> active_ptr;
