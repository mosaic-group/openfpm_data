/*
 * Point_unit_test.hpp
 *
 *  Created on: Aug 11, 2015
 *      Author: i-bird
 */

#ifndef SRC_SPACE_SHAPE_POINT_UNIT_TEST_HPP_
#define SRC_SPACE_SHAPE_POINT_UNIT_TEST_HPP_

#include <iostream>
#include "Space/Shape/Point.hpp"
#include "Grid/comb.hpp"
#include "util/convert.hpp"

BOOST_AUTO_TEST_SUITE( Point_test_suite )

BOOST_AUTO_TEST_CASE( Point_use )
{
	std::cout << "Point unit test start" << "\n";

	//! [assign and operation]

	Point<3,float> p;

	p.get(0) = 1.0;
	p.get(1) = 2.0;
	p.get(2) = 3.0;

	p = pmul(p,p) + p;

	BOOST_REQUIRE_EQUAL(p.get(0),2.0);
	BOOST_REQUIRE_EQUAL(p.get(1),6.0);
	BOOST_REQUIRE_EQUAL(p.get(2),12.0);

	// Combination 3D
	comb<3> c;
	c.c[0] = 1;
	c.c[1] = 1;
	c.c[2] = 1;
	Point<3,float> pc = toPoint<3,float>::convert(c);

	Point<3,float> one;
	// fill the point with one
	one.one();
	Point<3,float> middle = p / 2.0;

	// middle * (one + comb_p) = 0
	p = p + pmul(middle,(one + pc));

	BOOST_REQUIRE_EQUAL(p.get(0),4.0);
	BOOST_REQUIRE_EQUAL(p.get(1),12.0);
	BOOST_REQUIRE_EQUAL(p.get(2),24.0);

	//! [assign and operation]

	//! [norm]

	Point<3,float> pn({1.0,1.0,1.0});

	BOOST_REQUIRE_CLOSE(pn.norm(),1.732050,0.0001);

	//! [norm]

	// distance between two points

	Point<3,float> p1({1.0,1.0,1.0});
	Point<3,float> p2({2.0,3.0,4.0});

	BOOST_REQUIRE_EQUAL(p1.distance(p2),p2.distance(p1));
	BOOST_REQUIRE_CLOSE(p1.distance(p2),3.74165738677,0.0001);

	BOOST_REQUIRE_EQUAL(p1.distance2(p2),p2.distance2(p1));
	BOOST_REQUIRE_CLOSE(p1.distance2(p2),14.0,0.0001);

	std::cout << "Point unit test stop" << "\n";

	// check two point are equal

	BOOST_REQUIRE( p1 != p2 );
	BOOST_REQUIRE( p1 == p1 );

	// Check division and operator*

	p2 /= 2.0;
	p1 = p1 * 2.0;

	BOOST_REQUIRE_CLOSE(p2.get(0),1.0,0.001);
	BOOST_REQUIRE_CLOSE(p2.get(1),1.5,0.001);
	BOOST_REQUIRE_CLOSE(p2.get(2),2.0,0.001);

	BOOST_REQUIRE_CLOSE(p1.get(0),2.0,0.001);
	BOOST_REQUIRE_CLOSE(p1.get(1),2.0,0.001);
	BOOST_REQUIRE_CLOSE(p1.get(2),2.0,0.001);


}

BOOST_AUTO_TEST_CASE( Point_expression_usage )
{
	float scal = 0.0;

	Point<3,float> p1({1.0,1.0,1.0});
	Point<3,float> p2({2.0,3.0,4.0});
	Point<3,float> p3({6.0,7.0,9.0});

	p3 = p1 + 2.0;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p1.get(i) + 2.0);}
	p3 = 2.0 + p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),2.0 + p2.get(i));}
	p3 = p2 + p1;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p2.get(i) + p1.get(i));}

	p3 = p1 - 2.0;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p1.get(i) - 2.0);}
	p3 = 2.0 - p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),2.0 - p2.get(i));}
	p3 = p2 - p1;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p2.get(i) - p1.get(i));}

	p3 = p1 * 2.0;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p1.get(i) * 2.0);}
	p3 = 2.0f * p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),2.0 * p2.get(i));}
	p3 = p2 * p1;
	for (size_t i = 0 ; i < 3 ; i++)	{scal += p1.get(i) * p2.get(i);}
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),scal);}

	p3 = p1 / 2.0;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p1.get(i) / 2.0);}
	p3 = 2.0 / p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),(float)2.0 / p2.get(i));}
	p3 = p2 / p1;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p2.get(i) / p1.get(i));}

	// Variuos combination 3 operator

	p3 = p1 + (p2 + p1);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p1.get(i) + (p2.get(i) + p1.get(i))) ;}
	p3 = (p1 + p2) + p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),(p1.get(i) + p2.get(i)) + p2.get(i)) ;}
	p3 = (p1 + p2) + (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),(p1.get(i) + p2.get(i)) + (p2.get(i) + p1.get(i))) ;}

	p3 = p2 - (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p2.get(i) - (p1.get(i) + p2.get(i))) ;}
	p3 = (p1 + p2) - p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),(p1.get(i) + p2.get(i)) - p2.get(i)) ;}
	p3 = (p1 + p2) - (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),(p1.get(i) + p2.get(i)) - (p1.get(i) + p2.get(i))) ;}

	scal = 0;
	p3 = p2 * (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{scal += p2.get(i) * (p2.get(i) + p1.get(i));}
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),scal) ;}
	p3 = (p1 + p2) * p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),scal) ;}
	scal = 0;
	p3 = (p1 + p2) * (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{scal += (p2.get(i) + p1.get(i)) * (p2.get(i) + p1.get(i));}
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),scal) ;}

	// norm test
	p3 = p1 / norm(p1);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p1.get(i) / p1.norm()) ;}
	// distance calculation
	float dist = norm(p1 - p2);
	float dist2 = 0.0;
	for (size_t i = 0 ; i < 3 ; i++)	{dist2 += (p1.get(i) - p2.get(i)) * (p1.get(i) - p2.get(i));}
	dist2 = sqrt(dist2);
	BOOST_REQUIRE_EQUAL(dist,dist2);
	float dist3 = distance(p1,p2);
	BOOST_REQUIRE_EQUAL(dist,dist3);

	p3 = p2 / (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p2.get(i) / (p2.get(i) + p1.get(i))) ;}
	p3 = (p1 + p2) / p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),(p2.get(i) + p1.get(i)) / p2.get(i) ) ;}
	p3 = (p1 + p2) / (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),(p2.get(i) + p1.get(i)) / (p2.get(i) + p1.get(i)) ) ;}


	// Point function test

	p3 = abs(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::abs(p1.get(i) + p2.get(i))) ;}
	p3 = exp(p1 - p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::exp(p1.get(i) - p2.get(i))) ;}
	p3 = exp2(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::exp2(p1.get(i) + p2.get(i))) ;}
	p3 = expm1(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::expm1(p1.get(i) + p2.get(i))) ;}
	p3 = log(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::log(p1.get(i) + p2.get(i))) ;}
	p3 = log10(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::log10(p1.get(i) + p2.get(i))) ;}
	p3 = log2(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::log2(p1.get(i) + p2.get(i))) ;}
	p3 = log1p(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::log1p(p1.get(i) + p2.get(i))) ;}
	p3 = sqrt(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::sqrt(p1.get(i) + p2.get(i))) ;}
	p3 = cbrt(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::cbrt(p1.get(i) + p2.get(i))) ;}
	p3 = sin(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::sin(p1.get(i) + p2.get(i))) ;}
	p3 = cos(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::cos(p1.get(i) + p2.get(i))) ;}
	p3 = tan(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::tan(p1.get(i) + p2.get(i))) ;}
	p3 = atan(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::atan(p1.get(i) + p2.get(i))) ;}
	p3 = asin(p1/5.0f + p2/6.0f);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::asin(p1.get(i)/5.0f + p2.get(i)/6.0f)) ;}
	p3 = acos(p1/5.0f + p2/6.0f);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::acos(p1.get(i)/5.0f + p2.get(i)/6.0f)) ;}
	p3 = sinh(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::sinh(p1.get(i) + p2.get(i))) ;}
	p3 = cosh(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::cosh(p1.get(i) + p2.get(i))) ;}
	p3 = tanh(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::tanh(p1.get(i) + p2.get(i))) ;}
	p3 = atanh(p1/5.0f + p2/6.0f);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::atanh(p1.get(i)/5.0f + p2.get(i)/6.0f)) ;}
	p3 = asinh(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::asinh(p1.get(i) + p2.get(i))) ;}
	p3 = acosh(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::acosh(p1.get(i) + p2.get(i))) ;}
	p3 = erf(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::erf(p1.get(i) + p2.get(i))) ;}
	p3 = erfc(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::erfc(p1.get(i) + p2.get(i))) ;}
	p3 = tgamma(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::tgamma(p1.get(i) + p2.get(i))) ;}
	p3 = lgamma(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::lgamma(p1.get(i) + p2.get(i))) ;}
	p3 = ceil(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::ceil(p1.get(i) + p2.get(i))) ;}
	p3 = floor(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::floor(p1.get(i) + p2.get(i))) ;}
	p3 = trunc(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::trunc(p1.get(i) + p2.get(i))) ;}
	p3 = round(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::round(p1.get(i) + p2.get(i))) ;}
	p3 = nearbyint(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::nearbyint(p1.get(i) + p2.get(i))) ;}
	p3 = rint(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),std::rint(p1.get(i) + p2.get(i))) ;}


	double tmp = 5.0 + (p1 * p2);
	BOOST_REQUIRE_EQUAL(tmp,14.0);
	p3 = 5.0 + (p1 * p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),14.0) ;}
	tmp = 5.0 - (p1 * p2);
	BOOST_REQUIRE_EQUAL(tmp,-4.0);
	p3 = 5.0 - (p1 * p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),-4.0) ;}
	tmp = 5.0 * (p1 * p2);
	BOOST_REQUIRE_EQUAL(tmp,45.0);
	p3 = 5.0 * (p1 * p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),45.0) ;}
	tmp = 5.0f / (p1 * p2);
	BOOST_REQUIRE_EQUAL(tmp,5.0f/9.0f);
	p3 = 5.0f / (p1 * p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),5.0f/9.0f) ;}

	float p[3] = {1.0,2.0,3.0};
	auto p_e = getExpr(p);

	p3 = p_e - p1;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p[i]-p1.get(i)) ;}
	p3 = p_e + p1;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3.get(i),p[i]+p1.get(i)) ;}

	tmp = norm(p_e);
	BOOST_REQUIRE_EQUAL(tmp,sqrt(14.0f));

	Point<3,float> p4({1.0,2.0,3.0});

	p4 += p2;
	BOOST_REQUIRE_EQUAL(p4.get(0),3.0);
	BOOST_REQUIRE_EQUAL(p4.get(1),5.0);
	BOOST_REQUIRE_EQUAL(p4.get(2),7.0);


	p4 = Point<3,float>({1.0,2.0,3.0});

	p4 += p2 + p1;

	BOOST_REQUIRE_EQUAL(p4.get(0),4.0);
	BOOST_REQUIRE_EQUAL(p4.get(1),6.0);
	BOOST_REQUIRE_EQUAL(p4.get(2),8.0);

	double s = 0.0;
	s += p1+(p1+p2);
	BOOST_REQUIRE_EQUAL(s,4.0);


}


BOOST_AUTO_TEST_CASE( Point_expression_usage_with_array )
{
	float scal = 0.0;

	float p1_p[] = {1.0,1.0,1.0};
	float p2_p[] = {2.0,3.0,4.0};
	float p3_p[] = {6.0,7.0,9.0};

	Point<3,float> pp1({1.0,1.0,1.0});
	Point<3,float> pp2({2.0,3.0,4.0});
	Point<3,float> pp3({6.0,7.0,9.0});

	auto p1 = getExpr(p1_p);
	auto p2 = getExpr(p2_p);
	auto p3 = getExpr(p3_p);

	p3 = p1 + 2.0;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p1_p[i] + 2.0);}
	p3 = 2.0 + p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],2.0 + p2_p[i]);}
	p3 = p2 + p1;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p2_p[i] + p1_p[i]);}

	p3 = p1 - 2.0;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p1_p[i] - 2.0);}
	p3 = 2.0 - p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],2.0 - p2_p[i]);}
	p3 = p2 - p1;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p2_p[i] - p1_p[i]);}

	p3 = p1 * 2.0;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p1_p[i] * 2.0);}
	p3 = 2.0f * p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],2.0 * p2_p[i]);}
	p3 = p2 * p1;
	for (size_t i = 0 ; i < 3 ; i++)	{scal += p1_p[i] * p2_p[i];}
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],scal);}

	p3 = p1 / 2.0;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p1_p[i] / 2.0);}
	p3 = 2.0 / p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],(float)2.0 / p2_p[i]);}
	p3 = p2 / p1;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p2_p[i] / p1_p[i]);}

	// Some 2 operator combination

	p3 = p1 + p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p1_p[i] + p2_p[i]);}
	p3 = p1 + pp2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p1_p[i] + pp2.get(i));}
	p3 = pp2 + p1;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],pp2.get(i) + p1_p[i]);}

	// Variuos combination 3 operator

	p3 = p1 + (p2 + p1);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p1_p[i] + (p2_p[i] + p1_p[i])) ;}
	p3 = (p1 + p2) + p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],(p1_p[i] + p2_p[i]) + p2_p[i]) ;}
	p3 = (p1 + p2) + (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],(p1_p[i] + p2_p[i]) + (p2_p[i] + p1_p[i])) ;}

	p3 = p2 - (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p2_p[i] - (p1_p[i] + p2_p[i])) ;}
	p3 = (p1 + p2) - p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],(p1_p[i] + p2_p[i]) - p2_p[i]) ;}
	p3 = (p1 + p2) - (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],(p1_p[i] + p2_p[i]) - (p1_p[i] + p2_p[i])) ;}

	scal = 0;
	p3 = p2 * (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{scal += p2_p[i] * (p2_p[i] + p1_p[i]);}
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],scal) ;}
	p3 = (p1 + p2) * p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],scal) ;}
	scal = 0;
	p3 = (p1 + p2) * (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{scal += (p2_p[i] + p1_p[i]) * (p2_p[i] + p1_p[i]);}
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],scal) ;}

	// norm test
	p3 = p1 / norm(p1);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p1_p[i] / sqrt(3.0f)) ;}
	// distance calculation
	float dist = norm(p1 - p2);
	float dist2 = 0.0;
	for (size_t i = 0 ; i < 3 ; i++)	{dist2 += (p1_p[i] - p2_p[i]) * (p1_p[i] - p2_p[i]);}
	dist2 = sqrt(dist2);
	BOOST_REQUIRE_EQUAL(dist,dist2);
	float dist3 = distance(p1,p2);
	BOOST_REQUIRE_EQUAL(dist,dist3);

	p3 = p2 / (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],p2_p[i] / (p2_p[i] + p1_p[i])) ;}
	p3 = (p1 + p2) / p2;
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],(p2_p[i] + p1_p[i]) / p2_p[i] ) ;}
	p3 = (p1 + p2) / (p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],(p2_p[i] + p1_p[i]) / (p2_p[i] + p1_p[i]) ) ;}


	// Point function test

	p3 = abs(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::abs(p1_p[i] + p2_p[i])) ;}
	p3 = exp(p1 - p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::exp(p1_p[i] - p2_p[i])) ;}
	p3 = exp2(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::exp2(p1_p[i] + p2_p[i])) ;}
	p3 = expm1(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::expm1(p1_p[i] + p2_p[i])) ;}
	p3 = log(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::log(p1_p[i] + p2_p[i])) ;}
	p3 = log10(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::log10(p1_p[i] + p2_p[i])) ;}
	p3 = log2(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::log2(p1_p[i] + p2_p[i])) ;}
	p3 = log1p(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::log1p(p1_p[i] + p2_p[i])) ;}
	p3 = sqrt(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::sqrt(p1_p[i] + p2_p[i])) ;}
	p3 = cbrt(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::cbrt(p1_p[i] + p2_p[i])) ;}
	p3 = sin(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::sin(p1_p[i] + p2_p[i])) ;}
	p3 = cos(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::cos(p1_p[i] + p2_p[i])) ;}
	p3 = tan(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::tan(p1_p[i] + p2_p[i])) ;}
	p3 = atan(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::atan(p1_p[i] + p2_p[i])) ;}
	p3 = asin(p1/5.0f + p2/6.0f);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::asin(p1_p[i]/5.0f + p2_p[i]/6.0f)) ;}
	p3 = acos(p1/5.0f + p2/6.0f);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::acos(p1_p[i]/5.0f + p2_p[i]/6.0f)) ;}
	p3 = sinh(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::sinh(p1_p[i] + p2_p[i])) ;}
	p3 = cosh(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::cosh(p1_p[i] + p2_p[i])) ;}
	p3 = tanh(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::tanh(p1_p[i] + p2_p[i])) ;}
	p3 = atanh(p1/5.0f + p2/6.0f);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::atanh(p1_p[i]/5.0f + p2_p[i]/6.0f)) ;}
	p3 = acosh(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::acosh(p1_p[i] + p2_p[i])) ;}
	p3 = erf(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::erf(p1_p[i] + p2_p[i])) ;}
	p3 = erfc(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::erfc(p1_p[i] + p2_p[i])) ;}
	p3 = tgamma(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::tgamma(p1_p[i] + p2_p[i])) ;}
	p3 = lgamma(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::lgamma(p1_p[i] + p2_p[i])) ;}
	p3 = ceil(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::ceil(p1_p[i] + p2_p[i])) ;}
	p3 = floor(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::floor(p1_p[i] + p2_p[i])) ;}
	p3 = trunc(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::trunc(p1_p[i] + p2_p[i])) ;}
	p3 = round(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::round(p1_p[i] + p2_p[i])) ;}
	p3 = nearbyint(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::nearbyint(p1_p[i] + p2_p[i])) ;}
	p3 = rint(p1 + p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],std::rint(p1_p[i] + p2_p[i])) ;}


	double tmp = 5.0 + (p1 * p2);
	BOOST_REQUIRE_EQUAL(tmp,14.0);
	p3 = 5.0 + (p1 * p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],14.0) ;}
	tmp = 5.0 - (p1 * p2);
	BOOST_REQUIRE_EQUAL(tmp,-4.0);
	p3 = 5.0 - (p1 * p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],-4.0) ;}
	tmp = 5.0 * (p1 * p2);
	BOOST_REQUIRE_EQUAL(tmp,45.0);
	p3 = 5.0 * (p1 * p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],45.0) ;}
	tmp = 5.0f / (p1 * p2);
	BOOST_REQUIRE_EQUAL(tmp,5.0f/9.0f);
	p3 = 5.0f / (p1 * p2);
	for (size_t i = 0 ; i < 3 ; i++)	{BOOST_REQUIRE_EQUAL(p3_p[i],5.0f/9.0f) ;}
}


BOOST_AUTO_TEST_SUITE_END()




#endif /* SRC_SPACE_SHAPE_POINT_UNIT_TEST_HPP_ */
