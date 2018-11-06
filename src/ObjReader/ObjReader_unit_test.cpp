#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ObjReader.hpp"

#ifdef HAVE_TINYOBJLOADER

BOOST_AUTO_TEST_SUITE( obj_reader_test_suite )

BOOST_AUTO_TEST_CASE( obj_reader_test_use)
{
	ObjReader<float> reader;

	reader.read("");
}

BOOST_AUTO_TEST_SUITE_END()


#endif
