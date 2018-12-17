#include "config.h"
#undef VERSION

#include <iostream>

#define BOOST_DISABLE_ASSERTS

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


// initialization function:
bool init_unit_test()
{
  return true;
}

// entry point:
int main(int argc, char* argv[])
{
  return boost::unit_test::unit_test_main( &init_unit_test, argc, argv );
}


#include "VCluster/VCluster.hpp"
#include "unit_test_init_cleanup_io.hpp"

#include "VCluster/VCluster.hpp"
#include "CSVWriter/CSVWriter_unit_tests.hpp"
#include "GraphMLWriter/GraphMLWriter_unit_tests.hpp"
#include "VTKWriter/VTKWriter_unit_tests.hpp"
//#include "HDF5_XdmfWriter/HDF5_XdmfWriter_unit_tests.hpp"
#include "Plot/Plot_unit_tests.hpp"
#include "RawReader/RawReader_unit_tests.hpp"
#include "HDF5_wr/HDF5_writer_unit_tests.hpp"
