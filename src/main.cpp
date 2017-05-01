#include "config.h"
#undef VERSION

#include <iostream>

#define BOOST_DISABLE_ASSERTS

#define BOOST_TEST_MODULE "C++ test module for OpenFPM_io project"
#include <boost/test/included/unit_test.hpp>

#include "VCluster/VCluster.hpp"
#include "unit_test_init_cleanup.hpp"

#include "VCluster/VCluster.hpp"
#include "CSVWriter/CSVWriter_unit_tests.hpp"
#include "GraphMLWriter/GraphMLWriter_unit_tests.hpp"
#include "VTKWriter/VTKWriter_unit_tests.hpp"
//#include "HDF5_XdmfWriter/HDF5_XdmfWriter_unit_tests.hpp"
#include "Plot/Plot_unit_tests.hpp"
#include "RawReader/RawReader_unit_tests.hpp"
#include "HDF5_wr/HDF5_writer_unit_tests.hpp"
