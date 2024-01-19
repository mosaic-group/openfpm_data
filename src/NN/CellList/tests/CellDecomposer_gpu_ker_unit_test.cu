#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "NN/CellList/cuda/CellDecomposer_gpu_ker.cuh"
#include "Space/Shape/Point.hpp"
#include "Vector/map_vector.hpp"


template<typename vector_type, typename celldec>
__global__ void check(vector_type vd, celldec cd, unsigned int id, Point<3,float> p)
{
	vd.template get<0>(id) = cd.getCell(p);
}

BOOST_AUTO_TEST_SUITE( CellDecomposer_gpu_test_suite )


BOOST_AUTO_TEST_CASE( CellDecomposer_gpu_test_use )
{
	//! Spacing
	openfpm::array<float,3> spacing_c = {0.1,0.1,0.1};

	//! \brief number of sub-divisions in each direction
	openfpm::array<unsigned int,3> div_c = {10,10,10};

	//! \brief cell offset
	openfpm::array<unsigned int,3> off = {2,2,2};

	Point<3,float> trans({0.0,0.0,0.0});

	shift_only<3,float> t(Matrix<3,float>::identity(),trans);

	CellDecomposer_gpu_ker<3,float,unsigned int,shift_only<3,float>> clk(spacing_c,div_c,off,t);

	openfpm::vector_gpu<aggregate<grid_key_dx<3,unsigned int>>> output(8);

	CUDA_LAUNCH_DIM3(check,1,1,output.toKernel(),clk,0,Point<3,float>({0.2,0.2,0.2}));
	CUDA_LAUNCH_DIM3(check,1,1,output.toKernel(),clk,1,Point<3,float>({0.1,0.2,0.3}));
	CUDA_LAUNCH_DIM3(check,1,1,output.toKernel(),clk,2,Point<3,float>({0.25,0.55,0.45}));
	CUDA_LAUNCH_DIM3(check,1,1,output.toKernel(),clk,3,Point<3,float>({0.15,0.15,0.95}));
	CUDA_LAUNCH_DIM3(check,1,1,output.toKernel(),clk,4,Point<3,float>({1.05,1.05,1.05}));
	CUDA_LAUNCH_DIM3(check,1,1,output.toKernel(),clk,5,Point<3,float>({1.15,1.15,1.15}));
	CUDA_LAUNCH_DIM3(check,1,1,output.toKernel(),clk,6,Point<3,float>({-0.05,-0.05,-0.05}));
	CUDA_LAUNCH_DIM3(check,1,1,output.toKernel(),clk,7,Point<3,float>({-0.15,-0.15,-0.15}));

	output.template deviceToHost<0>();

	grid_key_dx<3,unsigned int> k = output.template get<0>(0);
	BOOST_REQUIRE_EQUAL(k.get(0),4);
	BOOST_REQUIRE_EQUAL(k.get(1),4);
	BOOST_REQUIRE_EQUAL(k.get(2),4);
	k = output.template get<0>(1);
	BOOST_REQUIRE_EQUAL(k.get(0),3);
	BOOST_REQUIRE_EQUAL(k.get(1),4);
	BOOST_REQUIRE_EQUAL(k.get(2),5);
	k = output.template get<0>(2);
	BOOST_REQUIRE_EQUAL(k.get(0),4);
	BOOST_REQUIRE_EQUAL(k.get(1),7);
	BOOST_REQUIRE_EQUAL(k.get(2),6);
	k = output.template get<0>(3);
	BOOST_REQUIRE_EQUAL(k.get(0),3);
	BOOST_REQUIRE_EQUAL(k.get(1),3);
	BOOST_REQUIRE_EQUAL(k.get(2),11);
	k = output.template get<0>(4);
	BOOST_REQUIRE_EQUAL(k.get(0),12);
	BOOST_REQUIRE_EQUAL(k.get(1),12);
	BOOST_REQUIRE_EQUAL(k.get(2),12);
	k = output.template get<0>(5);
	BOOST_REQUIRE_EQUAL(k.get(0),13);
	BOOST_REQUIRE_EQUAL(k.get(1),13);
	BOOST_REQUIRE_EQUAL(k.get(2),13);
	k = output.template get<0>(6);
	BOOST_REQUIRE_EQUAL(k.get(0),1);
	BOOST_REQUIRE_EQUAL(k.get(1),1);
	BOOST_REQUIRE_EQUAL(k.get(2),1);
	k = output.template get<0>(7);
	BOOST_REQUIRE_EQUAL(k.get(0),0);
	BOOST_REQUIRE_EQUAL(k.get(1),0);
	BOOST_REQUIRE_EQUAL(k.get(2),0);
}

BOOST_AUTO_TEST_SUITE_END()
