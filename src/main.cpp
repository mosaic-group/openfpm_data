#define BOOST_DISABLE_ASSERTS

#include "config.h"
#define BOOST_TEST_MODULE "C++ test module for OpenFPM_data project"
#include <boost/test/included/unit_test.hpp>

#include <iostream>
#include <boost/mpl/int.hpp>
#include <typeinfo>
#include "util/ct_array.hpp"
#ifdef CUDA_GPU
#include "memory/CudaMemory.cuh"
#endif
#include "memory/HeapMemory.hpp"
#include "memory_conf.hpp"
#include "Grid/map_grid.hpp"
#include "Vector/map_vector.hpp"
#include "Graph/map_graph.hpp"

// Include tests

#include "util/util_test.hpp"
#include "Space/SpaceBox_unit_tests.hpp"
#include "Space/Shape/Box_unit_tests.hpp"
#include "NN/CellList/CellList_test.hpp"
#include "Vector/vector_unit_tests.hpp"
#include "hypercube_unit_test.hpp"
#include "Graph/graph_unit_tests.hpp"
#include "Grid/grid_unit_tests.hpp"

