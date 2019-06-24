//
// Created by tommaso on 18/06/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "../BlockGeometry.hpp"

template <unsigned int dim, typename BGT>
void testStandardLinearizations(BGT geometry)
{
    grid_key_dx<dim> origin({0,0});
    BOOST_REQUIRE_EQUAL(geometry.LinId(origin), 0);

    grid_key_dx<dim> block0a({7,0});
    BOOST_REQUIRE_EQUAL(geometry.LinId(block0a), 7);
    grid_key_dx<dim> block0b({0,1});
    BOOST_REQUIRE_EQUAL(geometry.LinId(block0b), 8);
    grid_key_dx<dim> block0c({7,7});
    BOOST_REQUIRE_EQUAL(geometry.LinId(block0c), 63);

    grid_key_dx<dim> block1a({8+7,0});
    BOOST_REQUIRE_EQUAL(geometry.LinId(block1a), 64+7);
    grid_key_dx<dim> block1b({8+0,1});
    BOOST_REQUIRE_EQUAL(geometry.LinId(block1b), 64+8);
    grid_key_dx<dim> block1c({8+7,7});
    BOOST_REQUIRE_EQUAL(geometry.LinId(block1c), 64+63);

    grid_key_dx<dim> block3a({7,8+0});
    BOOST_REQUIRE_EQUAL(geometry.LinId(block3a), (64*3)+7);
    grid_key_dx<dim> block3b({0,8+1});
    BOOST_REQUIRE_EQUAL(geometry.LinId(block3b), (64*3)+8);
    grid_key_dx<dim> block3c({7,8+7});
    BOOST_REQUIRE_EQUAL(geometry.LinId(block3c), (64*3)+63);
}

BOOST_AUTO_TEST_SUITE(BlockGeometry_tests)
    BOOST_AUTO_TEST_CASE(testLinId)
    {
        constexpr unsigned int dim = 2;
        const size_t blockGridDim[dim] = {3,7};
        BlockGeometry<dim, 8> geometry(blockGridDim);

        testStandardLinearizations<dim>(geometry);
    }

    BOOST_AUTO_TEST_CASE(testCopyConstructor)
    {
        constexpr unsigned int dim = 2;
        const size_t blockGridDim[dim] = {3,7};
        BlockGeometry<dim, 8> geometry0(blockGridDim);

        // Here copy-construct
        BlockGeometry<dim, 8> geometry(geometry0);

        // Then test...
        testStandardLinearizations<dim>(geometry);
    }

    BOOST_AUTO_TEST_CASE(testCopyAssignOp)
    {
        constexpr unsigned int dim = 2;
        const size_t blockGridDim[dim] = {3,7};
        BlockGeometry<dim, 8> geometry0(blockGridDim);

        // Here copy-assign
        const size_t blockGridDimOther[dim] = {3,7};
        BlockGeometry<dim, 8> geometry(blockGridDimOther);
        geometry = geometry0;

        // Then test...
        testStandardLinearizations<dim>(geometry);
    }
BOOST_AUTO_TEST_SUITE_END()
