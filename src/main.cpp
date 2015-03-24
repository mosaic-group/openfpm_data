#define BOOST_DISABLE_ASSERTS


#define BOOST_TEST_MODULE "C++ test module for OpenFPM_data project"
#include <boost/test/included/unit_test.hpp>

#include <iostream>
#include <boost/mpl/int.hpp>
#include <typeinfo>
#include <ct_array.hpp>
#include "memory/CudaMemory.cuh"
#include "memory/HeapMemory.hpp"
#include "memory_conf.hpp"
#include "Grid/map_grid.hpp"
#include "Vector/map_vector.hpp"
#include "Graph/map_graph.hpp"

// Include tests

#include "Vector/vector_unit_tests.hpp"
#include "hypercube_unit_test.hpp"
#include "Graph/graph_unit_tests.hpp"
#include "Grid/grid_unit_tests.hpp"

