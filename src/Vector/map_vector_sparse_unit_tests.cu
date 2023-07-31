#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "map_vector_sparse.hpp"

BOOST_AUTO_TEST_SUITE( sparse_vector_test )


BOOST_AUTO_TEST_CASE ( test_sparse_vector_use )
{
	openfpm::vector_sparse<aggregate<size_t,size_t[3],size_t[3][3]>> vs;

	aggregate<size_t,size_t[3],size_t[3][3]> bck;
	vs.template setBackground<0>(0);

	size_t v[3] = {0,0,0};

	vs.template setBackground<1>(v);

	size_t t[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

	vs.template setBackground<2>(t);

	vs.template insert<0>(5) = 5;
	vs.template insert<0>(8) = 8;
	vs.template insert<0>(25) = 25;
	vs.template insert<0>(84) = 84;
	vs.template insert<0>(54) = 54;
	vs.template insert<0>(81) = 81;
	vs.template insert<0>(57) = 57;
	vs.template insert<0>(83) = 83;
	vs.template insert<0>(85) = 85;
	vs.template insert<0>(82) = 82;
	vs.template insert<0>(35) = 35;
	vs.template insert<0>(28) = 28;

	gpu::ofp_context_t ctx;
	vs.template flush<sadd_<0>>(ctx);

	BOOST_REQUIRE_EQUAL(vs.get<0>(5),5);
	BOOST_REQUIRE_EQUAL(vs.get<0>(54),54);
	BOOST_REQUIRE_EQUAL(vs.get<0>(85),85);
	BOOST_REQUIRE_EQUAL(vs.get<0>(1000),0);

	//// Add new elements

	vs.template insert<0>(4) = 4;
	vs.template insert<0>(28) = 28;
	vs.template insert<0>(45) = 45;
	vs.template insert<0>(94) = 94;
	vs.template insert<0>(88) = 88;
	vs.template insert<0>(823) = 823;

	vs.template flush<sadd_<0>>(ctx);

	BOOST_REQUIRE_EQUAL(vs.get<0>(5),5);
	BOOST_REQUIRE_EQUAL(vs.get<0>(54),54);
	BOOST_REQUIRE_EQUAL(vs.get<0>(85),85);
	BOOST_REQUIRE_EQUAL(vs.get<0>(1000),0);
	BOOST_REQUIRE_EQUAL(vs.get<0>(4),4);
	BOOST_REQUIRE_EQUAL(vs.get<0>(28),56);
	BOOST_REQUIRE_EQUAL(vs.get<0>(45),45);
	BOOST_REQUIRE_EQUAL(vs.get<0>(823),823);

	// we try to overlap some numbers

	vs.template insert<0>(45) = 450;
	vs.template insert<0>(94) = 940;
	vs.template insert<0>(88) = 880;

	vs.template flush<sadd_<0>>(ctx);

	BOOST_REQUIRE_EQUAL(vs.get<0>(45),495);
	BOOST_REQUIRE_EQUAL(vs.get<0>(94),1034);
	BOOST_REQUIRE_EQUAL(vs.get<0>(88),968);

	// try to insert the same thing multiple time

	vs.template insert<0>(101) = 1850;
	vs.template insert<0>(45) = 1450;
	vs.template insert<0>(45) = 1940;
	vs.template insert<0>(45) = 1880;
	vs.template insert<0>(1) = 2050;

	vs.template flush<sadd_<0>>(ctx);

	BOOST_REQUIRE_EQUAL(vs.get<0>(101),1850);
	BOOST_REQUIRE_EQUAL(vs.get<0>(45),5765);
	BOOST_REQUIRE_EQUAL(vs.get<0>(1),2050);
}

BOOST_AUTO_TEST_SUITE_END()
