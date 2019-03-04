/*
 * Box_unit_tests.hpp
 *
 *  Created on: May 8, 2015
 *      Author: i-bird
 */

#ifndef BOX_UNIT_TESTS_HPP_
#define BOX_UNIT_TESTS_HPP_

BOOST_AUTO_TEST_SUITE( box_test )

BOOST_AUTO_TEST_CASE( box_use)
{
	std::cout << "Box unit test start" << "\n";

	//! [expand the box with some spacing]

	float spacing[2] = {0.1,0.1};

	Box<2,float> sp({1.0,1.0},{2.0,3.0});
	sp.expand(spacing);

	BOOST_REQUIRE_CLOSE(sp.getLow(0),1.0,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getLow(1),1.0,0.0001);

	BOOST_REQUIRE_CLOSE(sp.getHigh(0),2.1,0.0001);
	BOOST_REQUIRE_CLOSE(sp.getHigh(1),3.1,0.0001);

	//! [expand the box with some spacing]

	{
	//! [create an enclosing box]

	Box<3,float> box1({0.1,0.2,0.3},{1.0,1.1,1.3});
	Box<3,float> box2({0.5,0.6,0.7},{2.0,2.1,2.2});

	box1.enclose(box2);

	BOOST_REQUIRE_EQUAL(box1.getLow(0),0.1f);
	BOOST_REQUIRE_EQUAL(box1.getLow(1),0.2f);
	BOOST_REQUIRE_EQUAL(box1.getLow(2),0.3f);

	BOOST_REQUIRE_EQUAL(box1.getHigh(0),2.0f);
	BOOST_REQUIRE_EQUAL(box1.getHigh(1),2.1f);
	BOOST_REQUIRE_EQUAL(box1.getHigh(2),2.2f);

	//! [create an enclosing box]
	}

	// Create the smallest boxes between several boxes or
	// Create a box that is contained into more boxes centering them on p1

	{

	//! [Create the smallest boxes between several boxes]
	Box<3,float> box1({0.0,0.0,0.0},{1.0,1.1,1.3});
	Box<3,float> box2({0.5,2.0,0.5},{2.0,2.1,2.2});
	Box<3,float> box3({1.5,1.5,4.2},{5.0,5.1,5.2});

	box1.contained(box2);
	box1.contained(box3);

	BOOST_REQUIRE_CLOSE(box1.getHigh(0),1.0f,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1),0.1f,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2),1.0f,0.0001);

	//! [Create the smallest boxes between several boxes]
	}

	{
	//! [Enlarge the box]

	Box<3,float> box1({0.1,0.2,0.3},{1.0,1.1,1.3});
	Box<3,float> box2({-0.5,-0.6,-0.7},{0.5,0.6,0.7});

	box1.enlarge(box2);

	BOOST_REQUIRE_CLOSE(box1.getLow(0),-0.4,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(1),-0.4,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(2),-0.4,0.0001);

	BOOST_REQUIRE_CLOSE(box1.getHigh(0),1.5,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1),1.7,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2),2.0,0.0001);

	//! [Enlarge the box]
	}

	{
	//! [Enlarge the box with fixed P1]

	Box<3,float> box1({0.1,0.2,0.3},{1.0,1.1,1.3});
	Box<3,float> box2({-0.5,-0.6,-0.7},{0.5,0.6,0.7});

	box1.enlarge_fix_P1(box2);

	BOOST_REQUIRE_CLOSE(box1.getLow(0),0.1,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(1),0.2,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(2),0.3,0.0001);

	BOOST_REQUIRE_CLOSE(box1.getHigh(0),2.0,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1),2.3,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2),2.7,0.0001);

	//! [Enlarge the box with fixed P1]
	}

	{
	//! [Magnify the box]

	Box<3,float> box1({0.1,0.2,0.3},{1.0,1.1,1.3});

	box1.magnify(1.001);

	BOOST_REQUIRE_CLOSE(box1.getLow(0),0.1001,0.001);
	BOOST_REQUIRE_CLOSE(box1.getLow(1),0.2002,0.001);
	BOOST_REQUIRE_CLOSE(box1.getLow(2),0.3003,0.001);

	BOOST_REQUIRE_CLOSE(box1.getHigh(0),1.001,0.00001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1),1.1011,0.00001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2),1.3013,0.00001);

	// Check that the dimensions of the box are magnified of 0.1%
	BOOST_REQUIRE_CLOSE(box1.getHigh(0) - box1.getLow(0), 1.001 * 0.9 ,0.001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1) - box1.getLow(1), 1.001 * 0.9, 0.001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2) - box1.getLow(2), 1.001 * 1.0, 0.001);

	//! [Magnify the box]
	}

	{
	//! [Magnify the box with fixed P1]

	Box<3,float> box1({0.1,0.2,0.3},{1.0,1.1,1.3});

	box1.magnify_fix_P1(1.001);

	BOOST_REQUIRE_CLOSE(box1.getLow(0),0.1,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(1),0.2,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(2),0.3,0.0001);

	BOOST_REQUIRE_CLOSE(box1.getHigh(0),1.0009,0.00001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1),1.1009,0.00001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2),1.301,0.00001);

	// Check that the dimensions of the box are magnified of 0.1%
	BOOST_REQUIRE_CLOSE(box1.getHigh(0) - box1.getLow(0), 1.001 * 0.9 ,0.00001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1) - box1.getLow(1), 1.001 * 0.9, 0.00001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2) - box1.getLow(2), 1.001 * 1.0, 0.00001);

	//! [Magnify the box with fixed P1]
	}

	{
	//! [Translate a box]

	Box<3,float> box1({0.1,0.5,0.6},{1.0,1.2,1.4});
	Point<3,float> pnt({0.1,0.2,0.3});

	box1 -= pnt;

	BOOST_REQUIRE_CLOSE(box1.getLow(0),0.0,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(1),0.3,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(2),0.3,0.0001);

	BOOST_REQUIRE_CLOSE(box1.getHigh(0),0.9,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1),1.0,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2),1.1,0.0001);

	//! [Translate a box]

	box1 -= box1.getP2();

	BOOST_REQUIRE_CLOSE(box1.getLow(0),-0.9,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(1),-0.7,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(2),-0.8,0.0001);

	BOOST_REQUIRE_CLOSE(box1.getHigh(0),0.0,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1),0.0,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2),0.0,0.0001);

	box1 -= box1.getP1();

	BOOST_REQUIRE_CLOSE(box1.getLow(0),0.0,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(1),0.0,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getLow(2),0.0,0.0001);

	BOOST_REQUIRE_CLOSE(box1.getHigh(0),0.9,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(1),0.7,0.0001);
	BOOST_REQUIRE_CLOSE(box1.getHigh(2),0.8,0.0001);
	}

	{
	Box<2,size_t> invalid1({5,7},{3,9});
	Box<2,size_t> invalid2({5,11},{9,9});
	Box<2,size_t> invalid3({12,11},{9,9});

	Box<2,size_t> valid1({1,5},{6,9});
	Box<2,size_t> valid2({1,1},{1,1});

	bool val = invalid1.isValid();
	BOOST_REQUIRE_EQUAL(val,false);
	val = invalid2.isValid();
	BOOST_REQUIRE_EQUAL(invalid2.isValid(),false);
	val = invalid3.isValid();
	BOOST_REQUIRE_EQUAL(invalid3.isValid(),false);
	val = valid1.isValid();
	BOOST_REQUIRE_EQUAL(valid1.isValid(),true);
	val = valid2.isValid();
	BOOST_REQUIRE_EQUAL(valid2.isValid(),true);

	}

	{
	//! [Box operators]

	Box<2,float> box({0.5,0.6},{1.4,1.3});
	Point<2,float> p({3.0,4.0});
	Box<2,float> result;

	result = box * p;

	BOOST_REQUIRE_CLOSE(result.getHigh(0),4.2,0.0001);
	BOOST_REQUIRE_CLOSE(result.getHigh(1),5.2,0.0001);

	BOOST_REQUIRE_CLOSE(result.getLow(0),1.5,0.0001);
	BOOST_REQUIRE_CLOSE(result.getLow(1),2.4,0.0001);

	//! [Box operators]
	}

	{
	//! [Box is insideNP]

	Box<2,float> box({0.0,0.0},{1.0,1.0});
	Point<2,float> p1({0.0,0.0});
	Point<2,float> p2({0.5,0.5});
	Point<2,float> p3({1.0,0.5});

	BOOST_REQUIRE_EQUAL(box.isInsideNP(p1),true);
	BOOST_REQUIRE_EQUAL(box.isInsideNP(p2),true);
	BOOST_REQUIRE_EQUAL(box.isInsideNP(p3),false);

	//! [Box is insideNP]
	}

	{
	//! [Box is inside]

	Box<2,float> box({0.0,0.0},{1.0,1.0});
	Point<2,float> p1({0.0,0.0});
	Point<2,float> p2({0.5,0.5});
	Point<2,float> p3({1.0,0.5});

	BOOST_REQUIRE_EQUAL(box.isInside(p1),true);
	BOOST_REQUIRE_EQUAL(box.isInside(p2),true);
	BOOST_REQUIRE_EQUAL(box.isInside(p3),true);

	//! [Box is inside]
	}

	{
	//! [Box is insideNB]

	Box<2,float> box({0.0,0.0},{1.0,1.0});
	Point<2,float> p1({0.0,0.0});
	Point<2,float> p2({0.5,0.5});
	Point<2,float> p3({1.0,0.5});

	BOOST_REQUIRE_EQUAL(box.isInsideNB(p1),false);
	BOOST_REQUIRE_EQUAL(box.isInsideNB(p2),true);
	BOOST_REQUIRE_EQUAL(box.isInsideNB(p3),false);

	//! [Box is insideNB]
	}

	Box<3,float> b1({1.0,2.4,5.3},{6.0,7.9,9.9});
	Box<3,float> b2 = b1;

	bool ret = (b1 == b2);

	BOOST_REQUIRE_EQUAL(ret,true);

	b1.invalidate();

	BOOST_REQUIRE_EQUAL(b1.isValid(),false);

	std::cout << "Box unit test stop" << "\n";
}

BOOST_AUTO_TEST_CASE( box_min_distance_test )
{
	// 2D

	Box<2,float> b1({0.0,0.0},{1.0,1.0});

	// Tounching boxes
	Box<2,float> b2({-1.0,-1.0},{0.0,0.0});
	Box<2,float> b3({-1.0, 0.0},{0.0,1.0});
	Box<2,float> b4({-1.0, 1.0},{0.0,2.0});
	Box<2,float> b5({ 0.0,-1.0},{1.0,0.0});
	Box<2,float> b6({ 0.0, 0.0},{1.0,1.0});
	Box<2,float> b7({ 0.0, 1.0},{1.0,2.0});
	Box<2,float> b8({ 1.0,-1.0},{2.0,0.0});
	Box<2,float> b9({ 1.0, 0.0},{2.0,1.0});
	Box<2,float> b10({1.0, 1.0},{2.0,2.0});

	BOOST_REQUIRE_EQUAL(b1.min_distance(b2),0.0);
	BOOST_REQUIRE_EQUAL(b1.min_distance(b3),0.0);
	BOOST_REQUIRE_EQUAL(b1.min_distance(b4),0.0);
	BOOST_REQUIRE_EQUAL(b1.min_distance(b5),0.0);
	BOOST_REQUIRE_EQUAL(b1.min_distance(b6),0.0);
	BOOST_REQUIRE_EQUAL(b1.min_distance(b7),0.0);
	BOOST_REQUIRE_EQUAL(b1.min_distance(b8),0.0);
	BOOST_REQUIRE_EQUAL(b1.min_distance(b9),0.0);
	BOOST_REQUIRE_EQUAL(b1.min_distance(b10),0.0);

	// shift 0.1 on X

	Point<2,float> p({-0.1,0.0});

	b2 += p;
	b3 += p;
	b4 += p;
	b5 += p;
	b6 += p;
	b7 += p;
	b8 += p;
	b9 += p;
	b10 += p;

	BOOST_REQUIRE_CLOSE(b1.min_distance(b2),0.1,0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b3),0.1,0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b4),0.1,0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b5),0.1,0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b6),0.1,0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b7),0.1,0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b8),0.1,0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b9),0.1,0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b10),0.1,0.01);

	// shift out -0.1 on Y

	Point<2,float> p2({0.0,-0.1});

	// Tounching boxes
	b2 += p2;
	b3 += p2;
	b4 += p2;
	b5 += p2;
	b6 += p2;
	b7 += p2;
	b8 += p2;
	b9 += p2;
	b10 += p2;

	BOOST_REQUIRE_CLOSE(b1.min_distance(b2),sqrt(0.01 + 0.01),0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b3),sqrt(0.01 + 0.01),0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b4),sqrt(0.01 + 0.01),0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b5),sqrt(0.01 + 0.01),0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b6),sqrt(0.01 + 0.01),0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b7),sqrt(0.01 + 0.01),0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b8),sqrt(0.01 + 0.01),0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b9),sqrt(0.01 + 0.01),0.01);
	BOOST_REQUIRE_CLOSE(b1.min_distance(b10),sqrt(0.01 + 0.01),0.01);
}

BOOST_AUTO_TEST_CASE( box_is_inside_with_border )
{
	// 2D

	size_t bc_nn[] = {NON_PERIODIC,NON_PERIODIC};
	size_t bc_pn[] = {PERIODIC,NON_PERIODIC};
	size_t bc_np[] = {NON_PERIODIC,PERIODIC};
	size_t bc_pp[] = {PERIODIC,PERIODIC};

	Box<2,float> border({0.0,0.0},{1.0,1.0});

	// non-tounching box
	Box<2,float> b2({0.1,0.1},{0.2,0.2});
	Box<2,float> b3({0.1,0.1},{1.0,1.0});

	Point<2,float> p1({0.15,0.15});
	Point<2,float> p2({0.25,0.25});
	Point<2,float> p3({0.15,0.2});

	Point<2,float> p4({0.1,0.25});
	Point<2,float> p5({0.15,1.0});
	Point<2,float> p6({0.1,1.0});

	////// NON PERIODIC TEST

	bool result = b2.isInsideNP_with_border(p1,border,bc_nn);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b2.isInsideNP_with_border(p2,border,bc_nn);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b2.isInsideNP_with_border(p3,border,bc_nn);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b3.isInsideNP_with_border(p1,border,bc_nn);

	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p4,border,bc_nn);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p5,border,bc_nn);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p6,border,bc_nn);
	BOOST_REQUIRE_EQUAL(result,true);

	//////////////

	result = b2.isInsideNP_with_border(p1,border,bc_pn);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b2.isInsideNP_with_border(p2,border,bc_pn);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b2.isInsideNP_with_border(p3,border,bc_pn);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b3.isInsideNP_with_border(p1,border,bc_pn);

	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p4,border,bc_pn);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p5,border,bc_pn);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p6,border,bc_pn);
	BOOST_REQUIRE_EQUAL(result,true);

	////////////

	result = b2.isInsideNP_with_border(p1,border,bc_np);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b2.isInsideNP_with_border(p2,border,bc_np);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b2.isInsideNP_with_border(p3,border,bc_np);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b3.isInsideNP_with_border(p1,border,bc_np);

	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p4,border,bc_np);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p5,border,bc_np);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b3.isInsideNP_with_border(p6,border,bc_np);
	BOOST_REQUIRE_EQUAL(result,false);

	////////////

	result = b2.isInsideNP_with_border(p1,border,bc_pp);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b2.isInsideNP_with_border(p2,border,bc_pp);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b2.isInsideNP_with_border(p3,border,bc_pp);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b3.isInsideNP_with_border(p1,border,bc_pp);

	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p4,border,bc_pp);
	BOOST_REQUIRE_EQUAL(result,true);

	result = b3.isInsideNP_with_border(p5,border,bc_pp);
	BOOST_REQUIRE_EQUAL(result,false);

	result = b3.isInsideNP_with_border(p6,border,bc_pp);
	BOOST_REQUIRE_EQUAL(result,false);
}

BOOST_AUTO_TEST_SUITE_END()



#endif /* BOX_UNIT_TESTS_HPP_ */
