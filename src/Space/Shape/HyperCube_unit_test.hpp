#ifndef HYPERCUBE_UNIT_TEST_HPP
#define HYPERCUBE_UNIT_TEST_HPP

#include "Space/Shape/HyperCube.hpp"
#include "util/util_debug.hpp"

/*! \brief Check if the linearization is correct
 *
 * \param cmb vector of combination to checl
 *
 */

template<unsigned int dim> void check_lin(std::vector<comb<dim>> cmb)
{
	// Check the linearization
	for (size_t i = 0 ; i < cmb.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(HyperCube<dim>::LinId(cmb[i]),i);
	}
}

/*! \brief Check if the linearization of the permutation is correct
 *
 * \param cmb vector of combination to checl
 *
 */

template<unsigned int dim> void check_lin_perm(std::vector<comb<dim>> cmb)
{
	// Check the linearization
	for (size_t i = 0 ; i < cmb.size() ; i++)
	{
		BOOST_REQUIRE_EQUAL(HyperCube<dim>::LinPerm(cmb[i]),i);
	}
}

/*! \brief Check if the combinations are dinstict
 *
 * \param combs combinations to check
 *
 */

template<unsigned int dim> bool isDinstict(std::vector<comb<dim>> combs)
{
	// Check if the combination are dinstinct

	for (size_t i = 0 ; i < combs.size() ; i++)
	{
		for (size_t j = i+1 ; j < combs.size() ; j++)
		{
			if (combs[i] == combs[j])
				return false;
		}

	}

	return true;
}

/*! \brief isSubdecomposition check if combs are elements that exist in c as sub-elements
 *
 * Check if combs are lower dimensional elements of c
 *
 * \param combs elements to check
 * \param c element that must contain all the combs
 *
 */

template<unsigned int dim> bool isSubdecomposition(std::vector<comb<dim>> combs, comb<dim> c)
{
	for (size_t i = 0 ; i < combs.size() ; i++)
	{
		if (c.isSub(combs[i]) == false)
			return false;
	}

	return true;
}

template<unsigned int dim> bool isValid(std::vector<comb<dim>> combs)
{
	// Check if the combinations are valid

	for (size_t i = 0 ; i < combs.size() ; i++)
	{
		if (combs[i].isValid() == false)
			return false;
	}

	return true;
}

BOOST_AUTO_TEST_SUITE( Hyper_cube )

BOOST_AUTO_TEST_CASE( Hyper_cube_comb_use )
{
	{
	comb<3> c1({1,2,3});
	comb<3> c2({1,1,1});

	comb<3> c_ret = c1 - c2;

	BOOST_REQUIRE_EQUAL(c_ret.c[0],2);
	BOOST_REQUIRE_EQUAL(c_ret.c[1],1);
	BOOST_REQUIRE_EQUAL(c_ret.c[2],0);
	}

	{
	comb<3> c1({1,2,3});
	comb<3> c2;
	c2 = c1 & 0x01;

	BOOST_REQUIRE_EQUAL(c2.c[0],1);
	BOOST_REQUIRE_EQUAL(c2.c[1],0);
	BOOST_REQUIRE_EQUAL(c2.c[2],1);
	}

	{
	comb<3> c1({1,2,3});
	comb<3> c2({1,1,1});

	comb<3> c_ret = c1 + c2;

	BOOST_REQUIRE_EQUAL(c_ret.c[0],4);
	BOOST_REQUIRE_EQUAL(c_ret.c[1],3);
	BOOST_REQUIRE_EQUAL(c_ret.c[2],2);
	}

	{
	comb<3> c1({1,2,3});

	comb<3> c_ret = -c1;

	BOOST_REQUIRE_EQUAL(c_ret.c[0],-3);
	BOOST_REQUIRE_EQUAL(c_ret.c[1],-2);
	BOOST_REQUIRE_EQUAL(c_ret.c[2],-1);
	}

	{
	comb<3> c1({-1,1,0});

	comb<3> c_ret = c1.flip();

	BOOST_REQUIRE_EQUAL(c_ret.c[0],0);
	BOOST_REQUIRE_EQUAL(c_ret.c[1],-1);
	BOOST_REQUIRE_EQUAL(c_ret.c[2],-1);
	}
}

BOOST_AUTO_TEST_CASE( Hyper_cube_use)
{
	std::cout << "Hyper-cube unit test start" << "\n";
	// Test Hyper-Cube functionality

	size_t ele[8];

	//! [Get vertex and edge on a line]

	// Get the number of vertex (elements of dimension 0) of a line (dimension 1)
	ele[0] = HyperCube<1>::getNumberOfElements_R(0);
	// Get the number of edge (elements of dimension 1) of a line (dimension 1)
	ele[1] = HyperCube<1>::getNumberOfElements_R(1);

	// Get combination for each dimensions
	std::vector<comb<1>> v_c1_0 = HyperCube<1>::getCombinations_R(0);

	//! [Get vertex and edge on a line]

	// Check
	BOOST_REQUIRE_EQUAL(ele[0],2ul);
	BOOST_REQUIRE_EQUAL(ele[1],1ul);

	// Fill the expected vector

	comb<1> c[] = {{1},{-1}};

	// Check the linearization
	check_lin(v_c1_0);

	std::vector<comb<1>> v1_st;
	HyperCube<1>::BinPermutationsSt(v1_st);
	check_lin_perm(v1_st);

	boost_check_array(&c[0],&v_c1_0[0],2);

	//! [Get vertex edge and surfaces of a square]
	// Number of vertex
	ele[0] = HyperCube<2>::getNumberOfElements_R(0);
	// Number of edge
	ele[1] = HyperCube<2>::getNumberOfElements_R(1);
	// Number of faces
	ele[2] = HyperCube<2>::getNumberOfElements_R(2);

	// Get combination for vertex (1,1) (-1,1) (-1,1) (-1,-1)
	std::vector<comb<2>> v_c2_0 = HyperCube<2>::getCombinations_R(0);
	// Get combination for edges  (1,0) (-1,0) (0,1) (0,-1)
	std::vector<comb<2>> v_c2_1 = HyperCube<2>::getCombinations_R(1);

	//! [Get vertex edge and surfaces of a square]

	// Check
	BOOST_REQUIRE_EQUAL(ele[0],4ul);
	BOOST_REQUIRE_EQUAL(ele[1],4ul);
	BOOST_REQUIRE_EQUAL(ele[2],1ul);

	// Check combination

	comb<2> c2_0[] = {{1,1},{-1,1},{1,-1},{-1,-1}};
	comb<2> c2_1[] = {{0,1},{0,-1},{1,0},{-1,0}};
	check_lin(v_c2_0);
	check_lin(v_c2_1);

	std::vector<comb<2>> v2_st;
	HyperCube<2>::BinPermutationsSt(v2_st);
	check_lin_perm(v2_st);

	boost_check_array(&c2_0[0],&v_c2_0[0],4);
	boost_check_array(&c2_1[0],&v_c2_1[0],4);

	//! [Get vertex edge surfaces and volumes of a cube]
	// Number of vertex
	ele[0] = HyperCube<3>::getNumberOfElements_R(0);
	// Number of edge
	ele[1] = HyperCube<3>::getNumberOfElements_R(1);
	// Number of faces
	ele[2] = HyperCube<3>::getNumberOfElements_R(2);
	// Number of Cubes
	ele[3] = HyperCube<3>::getNumberOfElements_R(3);

	// Get combination for vertex
	std::vector<comb<3>> v_c3_0 = HyperCube<3>::getCombinations_R(0);
	// Get combinations for edge
	std::vector<comb<3>> v_c3_1 = HyperCube<3>::getCombinations_R(1);
	// Get combinations for surfaces
	std::vector<comb<3>> v_c3_2 = HyperCube<3>::getCombinations_R(2);

	//! [Get vertex edge surfaces and volumes of a cube]

	// Check
	BOOST_REQUIRE_EQUAL(ele[0],8ul);
	BOOST_REQUIRE_EQUAL(ele[1],12ul);
	BOOST_REQUIRE_EQUAL(ele[2],6ul);
	BOOST_REQUIRE_EQUAL(ele[3],1ul);

	// Check combination

	comb<3> c3_0[] = {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
	comb<3> c3_1[] = {{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0}};
	comb<3> c3_2[] = {{{0,0,1}},{{0,0,-1}},{{0,1,0}},{{0,-1,0}},{{1,0,0}},{{-1,0,0}}};
	check_lin(v_c3_0);
	check_lin(v_c3_1);
	check_lin(v_c3_2);

	std::vector<comb<3>> v3_st;
	HyperCube<3>::BinPermutationsSt(v3_st);
	check_lin_perm(v3_st);

	boost_check_array(&c3_0[0],&v_c3_0[0],8);
	boost_check_array(&c3_1[0],&v_c3_1[0],12);
	boost_check_array(&c3_2[0],&v_c3_2[0],6);

	// Tesseract
	// Number of vertex
	ele[0] = HyperCube<4>::getNumberOfElements_R(0);
	// Number of edge
	ele[1] = HyperCube<4>::getNumberOfElements_R(1);
	// Number of faces
	ele[2] = HyperCube<4>::getNumberOfElements_R(2);
	// Number of Cubes
	ele[3] = HyperCube<4>::getNumberOfElements_R(3);
	// Number of Tessaract
	ele[4] = HyperCube<4>::getNumberOfElements_R(4);

	// Get combination for each dimensions
	std::vector<comb<4>> v_c4_0 = HyperCube<4>::getCombinations_R(0);
	std::vector<comb<4>> v_c4_1 = HyperCube<4>::getCombinations_R(1);
	std::vector<comb<4>> v_c4_2 = HyperCube<4>::getCombinations_R(2);
	std::vector<comb<4>> v_c4_3 = HyperCube<4>::getCombinations_R(3);
	check_lin(v_c4_0);
	check_lin(v_c4_1);
	check_lin(v_c4_2);
	check_lin(v_c4_3);

	std::vector<comb<4>> v4_st;
	HyperCube<4>::BinPermutationsSt(v4_st);
	check_lin_perm(v1_st);

	// Check
	BOOST_REQUIRE_EQUAL(ele[0],16ul);
	BOOST_REQUIRE_EQUAL(ele[1],32ul);
	BOOST_REQUIRE_EQUAL(ele[2],24ul);
	BOOST_REQUIRE_EQUAL(ele[3],8ul);
	BOOST_REQUIRE_EQUAL(ele[4],1ul);

	// Penteract
	// Number of vertex
	ele[0] = HyperCube<5>::getNumberOfElements_R(0);
	// Number of edge
	ele[1] = HyperCube<5>::getNumberOfElements_R(1);
	// Number of faces
	ele[2] = HyperCube<5>::getNumberOfElements_R(2);
	// Number of Cubes
	ele[3] = HyperCube<5>::getNumberOfElements_R(3);
	// Number of Tessaract
	ele[4] = HyperCube<5>::getNumberOfElements_R(4);
	// Number of Penteract
	ele[5] = HyperCube<5>::getNumberOfElements_R(5);

	// Get combination for each dimensions
	std::vector<comb<5>> v_c5_0 = HyperCube<5>::getCombinations_R(0);
	std::vector<comb<5>> v_c5_1 = HyperCube<5>::getCombinations_R(1);
	std::vector<comb<5>> v_c5_2 = HyperCube<5>::getCombinations_R(2);
	std::vector<comb<5>> v_c5_3 = HyperCube<5>::getCombinations_R(3);
	std::vector<comb<5>> v_c5_4 = HyperCube<5>::getCombinations_R(4);
	check_lin(v_c5_0);
	check_lin(v_c5_1);
	check_lin(v_c5_2);
	check_lin(v_c5_3);
	check_lin(v_c5_4);

	std::vector<comb<5>> v5_st;
	HyperCube<5>::BinPermutationsSt(v5_st);
	check_lin_perm(v5_st);

	// Check
	BOOST_REQUIRE_EQUAL(ele[0],32ul);
	BOOST_REQUIRE_EQUAL(ele[1],80ul);
	BOOST_REQUIRE_EQUAL(ele[2],80ul);
	BOOST_REQUIRE_EQUAL(ele[3],40ul);
	BOOST_REQUIRE_EQUAL(ele[4],10ul);
	BOOST_REQUIRE_EQUAL(ele[5],1ul);


	// Test SubHypercube 2D and 3D

	//! [Get the vertices of a square]
	std::vector<comb<2>> sc2_0 = HyperCube<2>::getCombinations_R(0);
	//! [Get the vertices of a square]

	for (size_t i = 0 ; i < sc2_0.size() ; i++)
	{
		// Expecting one element equal to c2[i]
		//! [Getting the vertex as a sub-hypercube]
		std::vector<comb<2>> combs = SubHyperCube<2,0>::getCombinations_R(sc2_0[i],0);
		//! [Getting the vertex as a sub-hypercube]
		BOOST_REQUIRE_EQUAL(combs.size(),1ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc2_0[i]),true);
	}

	//! [Getting the edge of a square]
	std::vector<comb<2>> sc2_1 = HyperCube<2>::getCombinations_R(1);
	//! [Getting the edge of a square]

	for (size_t i = 0 ; i < sc2_1.size() ; i++)
	{
		// Expecting two elements, valid, distinct,  sub-decomposition of c2[i]

		//! [Getting the vertex of the line of the original square]
		std::vector<comb<2>> combs = SubHyperCube<2,1>::getCombinations_R(sc2_1[i],0);
		//! [Getting the vertex of the line of the original square]
		BOOST_REQUIRE_EQUAL(combs.size(),2ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc2_1[i]),true);

		// Expecting one element, valid, distinct,  sub-decomposition of c2[i]
		//! [Getting the edge(line) of the line of the original square]
		combs = SubHyperCube<2,1>::getCombinations_R(sc2_1[i],1);
		//! [Getting the edge(line) of the line of the original square]
		BOOST_REQUIRE_EQUAL(combs.size(),1ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc2_1[i]),true);
	}

	//! [Getting the square of a square]
	std::vector<comb<2>> sc2_2 = HyperCube<2>::getCombinations_R(2);
	//! [Getting the square of a square]

	for (size_t i = 0 ; i < sc2_2.size() ; i++)
	{
		// Expecting two elements, valid, distinct,  sub-decomposition of sc2_2[i]
		//! [Getting the vertices of the square of the original square]
		std::vector<comb<2>> combs = SubHyperCube<2,2>::getCombinations_R(sc2_2[i],0);
		//! [Getting the vertices of the square of the original square]
		BOOST_REQUIRE_EQUAL(combs.size(),4ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc2_2[i]),true);

		// Expecting one element, valid, distinct,  sub-decomposition of c2[i]
		//! [Getting the lines of the square of the original square]
		combs = SubHyperCube<2,2>::getCombinations_R(sc2_2[i],1);
		//! [Getting the lines of the square of the original square]
		BOOST_REQUIRE_EQUAL(combs.size(),4ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc2_2[i]),true);

		// Expecting one element, valid, distinct,  sub-decomposition of c2[i]
		//! [Getting the square of the square of the original square]
		combs = SubHyperCube<2,2>::getCombinations_R(sc2_2[i],2);
		//! [Getting the square of the square of the original square]
		BOOST_REQUIRE_EQUAL(combs.size(),1ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc2_2[i]),true);
	}

	////////////// 3D ////////////////

	//! [Getting the vertices of the cube]
	std::vector<comb<3>> sc3_0 = HyperCube<3>::getCombinations_R(0);
	//! [Getting the vertices of the cube]

	for (size_t i = 0 ; i < sc3_0.size() ; i++)
	{
		// Expecting one element equal to sc3[i]
		//! [Getting the vertices of the vertices of the cube]
		std::vector<comb<3>> combs = SubHyperCube<3,0>::getCombinations_R(sc3_0[i],0);
		//! [Getting the vertices of the vertices of the cube]
		BOOST_REQUIRE_EQUAL(combs.size(),1ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
//		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc3_0[i]),true);
	}

	//! [Getting the edges of the cube]
	std::vector<comb<3>> sc3_1 = HyperCube<3>::getCombinations_R(1);
	//! [Getting the edges of the cube]

	for (size_t i = 0 ; i < sc3_1.size() ; i++)
	{
		// Expecting two elements, valid, distinct,  sub-decomposition of sc3[i]
		//! [Getting the vertices of the edge of the cube]
		std::vector<comb<3>> combs = SubHyperCube<3,1>::getCombinations_R(sc3_1[i],0);
		//! [Getting the vertices of the vertices of the cube]
		BOOST_REQUIRE_EQUAL(combs.size(),2ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc3_1[i]),true);

		// Expecting one element, valid, distinct,  sub-decomposition of c2[i]
		//! [Getting the edges of the edges of the cube]
		combs = SubHyperCube<3,1>::getCombinations_R(sc3_1[i],1);
		//! [Getting the edges of the edges of the cube]
		BOOST_REQUIRE_EQUAL(combs.size(),1ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc3_1[i]),true);
	}

	//! [Getting the surfaces of the cube]
	std::vector<comb<3>> sc3_2 = HyperCube<3>::getCombinations_R(2);
	//! [Getting the surfaces of the cube]

	for (size_t i = 0 ; i < sc3_2.size() ; i++)
	{
		// Expecting two elements, valid, distinct,  sub-decomposition of sc3_2[i]
		//! [Getting the vertices of one surface of the cube]
		std::vector<comb<3>> combs = SubHyperCube<3,2>::getCombinations_R(sc3_2[i],0);
		//! [Getting the vertices of one surface of the cube]
		BOOST_REQUIRE_EQUAL(combs.size(),4ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc3_2[i]),true);

		// Expecting one element, valid, distinct,  sub-decomposition of c2[i]
		//! [Getting the edges of the surfaces of the cube]
		combs = SubHyperCube<3,2>::getCombinations_R(sc3_2[i],1);
		//! [Getting the edges of the surfaces of the cube]
		BOOST_REQUIRE_EQUAL(combs.size(),4ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc3_2[i]),true);

		// Expecting one element, valid, distinct,  sub-decomposition of c2[i]
		//! [Getting the surfaces of the surfaces of the cube]
		combs = SubHyperCube<3,2>::getCombinations_R(sc3_2[i],2);
		//! [Getting the surfaces of the surfaces of the cube]
		BOOST_REQUIRE_EQUAL(combs.size(),1ul);
		BOOST_REQUIRE_EQUAL(isDinstict(combs),true);
		BOOST_REQUIRE_EQUAL(isValid(combs),true);
		BOOST_REQUIRE_EQUAL(isSubdecomposition(combs,sc3_2[i]),true);
	}

	// High dimension test
	size_t bc[50];

	for (size_t i = 0 ; i < 50 ; i++)
	{bc[i] = NON_PERIODIC;}

	bc[0] = PERIODIC;
	bc[11] = PERIODIC;
	bc[17] = PERIODIC;

	std::vector<comb<50>> sc50_2 = HyperCube<50>::getCombinations_R_bc(49,bc);

	BOOST_REQUIRE_EQUAL(sc50_2.size(),6ul);

	BOOST_REQUIRE_EQUAL(sc50_2[0].c[0],1);
	BOOST_REQUIRE_EQUAL(sc50_2[0].c[11],0);
	BOOST_REQUIRE_EQUAL(sc50_2[0].c[17],0);

	BOOST_REQUIRE_EQUAL(sc50_2[1].c[0],-1);
	BOOST_REQUIRE_EQUAL(sc50_2[1].c[11],0);
	BOOST_REQUIRE_EQUAL(sc50_2[1].c[17],0);

	BOOST_REQUIRE_EQUAL(sc50_2[2].c[0],0);
	BOOST_REQUIRE_EQUAL(sc50_2[2].c[11],1);
	BOOST_REQUIRE_EQUAL(sc50_2[2].c[17],0);

	BOOST_REQUIRE_EQUAL(sc50_2[3].c[0],0);
	BOOST_REQUIRE_EQUAL(sc50_2[3].c[11],-1);
	BOOST_REQUIRE_EQUAL(sc50_2[3].c[17],0);

	BOOST_REQUIRE_EQUAL(sc50_2[4].c[0],0);
	BOOST_REQUIRE_EQUAL(sc50_2[4].c[11],0);
	BOOST_REQUIRE_EQUAL(sc50_2[4].c[17],1);

	BOOST_REQUIRE_EQUAL(sc50_2[5].c[0],0);
	BOOST_REQUIRE_EQUAL(sc50_2[5].c[11],0);
	BOOST_REQUIRE_EQUAL(sc50_2[5].c[17],-1);

	std::vector<comb<50>> sc50_3 = HyperCube<50>::getCombinations_R_bc(48,bc);

	BOOST_REQUIRE_EQUAL(sc50_3.size(),12ul);

	BOOST_REQUIRE_EQUAL(sc50_3[0].c[0],1);
	BOOST_REQUIRE_EQUAL(sc50_3[0].c[11],1);
	BOOST_REQUIRE_EQUAL(sc50_3[0].c[17],0);

	BOOST_REQUIRE_EQUAL(sc50_3[1].c[0],1);
	BOOST_REQUIRE_EQUAL(sc50_3[1].c[11],-1);
	BOOST_REQUIRE_EQUAL(sc50_3[1].c[17],0);

	BOOST_REQUIRE_EQUAL(sc50_3[2].c[0],-1);
	BOOST_REQUIRE_EQUAL(sc50_3[2].c[11],1);
	BOOST_REQUIRE_EQUAL(sc50_3[2].c[17],0);

	BOOST_REQUIRE_EQUAL(sc50_3[3].c[0],-1);
	BOOST_REQUIRE_EQUAL(sc50_3[3].c[11],-1);
	BOOST_REQUIRE_EQUAL(sc50_3[3].c[17],0);

	BOOST_REQUIRE_EQUAL(sc50_3[4].c[0],1);
	BOOST_REQUIRE_EQUAL(sc50_3[4].c[11],0);
	BOOST_REQUIRE_EQUAL(sc50_3[4].c[17],1);

	BOOST_REQUIRE_EQUAL(sc50_3[5].c[0],1);
	BOOST_REQUIRE_EQUAL(sc50_3[5].c[11],0);
	BOOST_REQUIRE_EQUAL(sc50_3[5].c[17],-1);

	BOOST_REQUIRE_EQUAL(sc50_3[6].c[0],-1);
	BOOST_REQUIRE_EQUAL(sc50_3[6].c[11],0);
	BOOST_REQUIRE_EQUAL(sc50_3[6].c[17],1);

	BOOST_REQUIRE_EQUAL(sc50_3[7].c[0],-1);
	BOOST_REQUIRE_EQUAL(sc50_3[7].c[11],0);
	BOOST_REQUIRE_EQUAL(sc50_3[7].c[17],-1);

	BOOST_REQUIRE_EQUAL(sc50_3[8].c[0],0);
	BOOST_REQUIRE_EQUAL(sc50_3[8].c[11],1);
	BOOST_REQUIRE_EQUAL(sc50_3[8].c[17],1);

	BOOST_REQUIRE_EQUAL(sc50_3[9].c[0],0);
	BOOST_REQUIRE_EQUAL(sc50_3[9].c[11],1);
	BOOST_REQUIRE_EQUAL(sc50_3[9].c[17],-1);

	BOOST_REQUIRE_EQUAL(sc50_3[10].c[0],0);
	BOOST_REQUIRE_EQUAL(sc50_3[10].c[11],-1);
	BOOST_REQUIRE_EQUAL(sc50_3[10].c[17],1);

	BOOST_REQUIRE_EQUAL(sc50_3[11].c[0],0);
	BOOST_REQUIRE_EQUAL(sc50_3[11].c[11],-1);
	BOOST_REQUIRE_EQUAL(sc50_3[11].c[17],-1);

	std::vector<comb<50>> sc50_4 = HyperCube<50>::getCombinations_R_bc(47,bc);

	BOOST_REQUIRE_EQUAL(sc50_4.size(),8ul);

	BOOST_REQUIRE_EQUAL(sc50_4[0].c[0],1);
	BOOST_REQUIRE_EQUAL(sc50_4[0].c[11],1);
	BOOST_REQUIRE_EQUAL(sc50_4[0].c[17],1);

	BOOST_REQUIRE_EQUAL(sc50_4[1].c[0],1);
	BOOST_REQUIRE_EQUAL(sc50_4[1].c[11],1);
	BOOST_REQUIRE_EQUAL(sc50_4[1].c[17],-1);

	BOOST_REQUIRE_EQUAL(sc50_4[2].c[0],1);
	BOOST_REQUIRE_EQUAL(sc50_4[2].c[11],-1);
	BOOST_REQUIRE_EQUAL(sc50_4[2].c[17],1);

	BOOST_REQUIRE_EQUAL(sc50_4[3].c[0],1);
	BOOST_REQUIRE_EQUAL(sc50_4[3].c[11],-1);
	BOOST_REQUIRE_EQUAL(sc50_4[3].c[17],-1);

	BOOST_REQUIRE_EQUAL(sc50_4[4].c[0],-1);
	BOOST_REQUIRE_EQUAL(sc50_4[4].c[11],1);
	BOOST_REQUIRE_EQUAL(sc50_4[4].c[17],1);

	BOOST_REQUIRE_EQUAL(sc50_4[5].c[0],-1);
	BOOST_REQUIRE_EQUAL(sc50_4[5].c[11],1);
	BOOST_REQUIRE_EQUAL(sc50_4[5].c[17],-1);

	BOOST_REQUIRE_EQUAL(sc50_4[6].c[0],-1);
	BOOST_REQUIRE_EQUAL(sc50_4[6].c[11],-1);
	BOOST_REQUIRE_EQUAL(sc50_4[6].c[17],1);

	BOOST_REQUIRE_EQUAL(sc50_4[7].c[0],-1);
	BOOST_REQUIRE_EQUAL(sc50_4[7].c[11],-1);
	BOOST_REQUIRE_EQUAL(sc50_4[7].c[17],-1);

	std::cout << "Hyper-cube unit test end" << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

#endif
