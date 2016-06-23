/*
 * Packer_unit_tests.hpp
 *
 *  Created on: Jul 15, 2015
 *      Author: Pietro Incardona
 */

#ifndef SRC_PACKER_UNIT_TESTS_HPP_
#define SRC_PACKER_UNIT_TESTS_HPP_

#include "Pack_selector.hpp"
#include "Packer.hpp"
#include "Unpacker.hpp"
#include "Grid/grid_util_test.hpp"
#include <iostream>
#include "Vector/vector_test_util.hpp"
#include "data_type/aggregate.hpp"

BOOST_AUTO_TEST_SUITE( packer_unpacker )

BOOST_AUTO_TEST_CASE ( packer_unpacker_test )
{
	//! [Pack selector usage]

	int val = Pack_selector<unsigned char>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);
	val = Pack_selector<char>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);
	val = Pack_selector<short>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);
	val = Pack_selector<unsigned short>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);
	val = Pack_selector<int>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);
	val = Pack_selector<unsigned int>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);
	val = Pack_selector<long int>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);
	val = Pack_selector<unsigned long int>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);
	val = Pack_selector<float>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);
	val = Pack_selector<double>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_PRIMITIVE);

	val = Pack_selector<Point_test<float>>::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_OBJECTS_WITH_POINTER_CHECK);


	val = Pack_selector< openfpm::vector<Point_test<float>> >::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_VECTOR);
	val = Pack_selector< grid_cpu<3,Point_test<float>> >::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_GRID);
	val = Pack_selector< encapc<3,Point_test<float>, memory_traits_lin<Point_test<float>>::type > >::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_ENCAP_OBJECTS);

	struct test_s
	{
		float a;
		float b;

		static bool noPointers() {return true;}
	};

	val = Pack_selector< test_s >::value;
	BOOST_REQUIRE_EQUAL(val,PACKER_OBJECTS_WITH_POINTER_CHECK);

	//! [Pack selector usage]

	{

	//! [Pack into a message primitives objects vectors and grids]

	typedef Point_test<float> pt;

	// Create all the objects we want to pack
	unsigned char uc = 1;
	char c = 2;
	short s = 3;
	unsigned short us = 4;
	int i = 5;
	unsigned int ui = 6;
	long int li = 7;
	unsigned long int uli = 8;
	float f = 9;
	double d = 10;

	openfpm::vector<Point_test<float>> v = allocate_openfpm<openfpm::vector<Point_test<float>>>(1024);

	Point_test<float> p;
	p.fill();

	size_t sz[] = {16,16,16};
	grid_cpu<3,Point_test<float>> g(sz);
	g.setMemory();
	fill_grid<3>(g);
	grid_key_dx_iterator_sub<3> sub(g.getGrid(),{1,2,3},{5,6,7});

	// Here we start to push all the allocations required to pack all the data

	std::vector<size_t> pap_prp;

	Packer<unsigned char,HeapMemory>::packRequest(uc,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(unsigned char));
	Packer<char,HeapMemory>::packRequest(c,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(char));
	Packer<short,HeapMemory>::packRequest(s,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(short));
	Packer<unsigned short,HeapMemory>::packRequest(us,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(unsigned short));
	Packer<int,HeapMemory>::packRequest(i,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(int));
	Packer<unsigned int,HeapMemory>::packRequest(ui,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(unsigned int));
	Packer<long int,HeapMemory>::packRequest(li,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(long int));
	Packer<long unsigned int,HeapMemory>::packRequest(uli,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(long unsigned int));
	Packer<float,HeapMemory>::packRequest(f,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(float));
	Packer<double,HeapMemory>::packRequest(d,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(double));
	Packer<Point_test<float>,HeapMemory>::packRequest(p,pap_prp);
	Packer<openfpm::vector<Point_test<float>>,HeapMemory>::packRequest<pt::x,pt::v>(v,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],(sizeof(float) + sizeof(float[3])) * v.size());
	Packer<grid_cpu<3,Point_test<float>>,HeapMemory>::packRequest<pt::x,pt::v>(g,sub,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],(sizeof(float) + sizeof(float[3])) * sub.getVolume());

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	Pack_stat sts;

	// try to pack
	Packer<unsigned char,HeapMemory>::pack(mem,1,sts);
	Packer<char,HeapMemory>::pack(mem,2,sts);
	Packer<short,HeapMemory>::pack(mem,3,sts);
	Packer<unsigned short,HeapMemory>::pack(mem,4,sts);
	Packer<int,HeapMemory>::pack(mem,5,sts);
	Packer<unsigned int, HeapMemory>::pack(mem,6,sts);
	Packer<long int,HeapMemory>::pack(mem,7,sts);
	Packer<long unsigned int,HeapMemory>::pack(mem,8,sts);
	Packer<float,HeapMemory>::pack(mem,9,sts);
	Packer<double,HeapMemory>::pack(mem,10,sts);
	Packer<Point_test<float>,HeapMemory>::pack(mem,p,sts);
	Packer<openfpm::vector<Point_test<float>>,HeapMemory>::pack<pt::x,pt::v>(mem,v,sts);
	Packer<grid_cpu<3,Point_test<float>>,HeapMemory>::pack<pt::x,pt::v>(mem,g,sub,sts);

	//! [Pack into a message primitives objects vectors and grids]

	//! [Unpack a message into primitives objects vectors and grids]

	Unpack_stat ps;

	unsigned char uc2;
	Unpacker<unsigned char,HeapMemory>::unpack(mem,uc2,ps);
	char c2;
	Unpacker<char,HeapMemory>::unpack(mem,c2,ps);
	short s2;
	Unpacker<short,HeapMemory>::unpack(mem,s2,ps);
	unsigned short us2;
	Unpacker<unsigned short,HeapMemory>::unpack(mem,us2,ps);
	int i2;
	Unpacker<int,HeapMemory>::unpack(mem,i2,ps);
	unsigned int ui2;
	Unpacker<unsigned int,HeapMemory>::unpack(mem,ui2,ps);
	long int li2;
	Unpacker<long int,HeapMemory>::unpack(mem,li2,ps);
	unsigned long int uli2;
	Unpacker<unsigned long int,HeapMemory>::unpack(mem,uli2,ps);
	float f2;
	Unpacker<float,HeapMemory>::unpack(mem,f2,ps);
	double d2;
	Unpacker<double,HeapMemory>::unpack(mem,d2,ps);

	// Unpack the point and check
	Point_test<float> p_test;
	Unpacker<Point_test<float>,HeapMemory>::unpack(mem,p_test,ps);

	// Unpack the vector and check
	openfpm::vector<Point_test<float>> v_test;
	v_test.resize(v.size());
	Unpacker<openfpm::vector<Point_test<float>>,HeapMemory>::unpack<pt::x,pt::v>(mem,v_test,ps);

	//! [Unpack a message into primitives objects vectors and grids]

	BOOST_REQUIRE_EQUAL(uc2,uc);
	BOOST_REQUIRE_EQUAL(c2,c);
	BOOST_REQUIRE_EQUAL(s2,s);
	BOOST_REQUIRE_EQUAL(us2,us);
	BOOST_REQUIRE_EQUAL(i2,i);
	BOOST_REQUIRE_EQUAL(ui2,ui);
	BOOST_REQUIRE_EQUAL(li2,li);
	BOOST_REQUIRE_EQUAL(uli2,uli);
	BOOST_REQUIRE_EQUAL(f2,f);
	BOOST_REQUIRE_EQUAL(d2,d);

	bool val = (p_test == p);
	BOOST_REQUIRE_EQUAL(true,val);

	auto it = v_test.getIterator();

	while (it.isNext())
	{
		float f1 = v_test.template get<pt::x>(it.get());
		float f2 = v.template get<pt::x>(it.get());

		BOOST_REQUIRE_EQUAL(f1,f2);

		for (size_t i = 0 ; i < 3 ; i++)
		{
			f1 = v_test.template get<pt::v>(it.get())[i];
			f2 = v.template get<pt::v>(it.get())[i];

			BOOST_REQUIRE_EQUAL(f1,f2);
		}

		++it;
	}

	// Unpack the grid and check

	size_t sz2[] = {16,16,16};
	grid_cpu<3,Point_test<float>> g_test(sz2);
	g_test.setMemory();
	grid_key_dx_iterator_sub<3> sub2(g_test.getGrid(),{1,2,3},{5,6,7});

	Unpacker<grid_cpu<3,Point_test<float>>,HeapMemory>::unpack<pt::x,pt::v>(mem,sub2,g_test,ps);

	// Check the unpacked grid
	sub2.reset();

	while (sub2.isNext())
	{
		float f1 = g_test.template get<pt::x>(sub2.get());
		float f2 = g.template get<pt::x>(sub2.get());

		BOOST_REQUIRE_EQUAL(f1,f2);

		for (size_t i = 0 ; i < 3 ; i++)
		{
			f1 = g_test.template get<pt::v>(sub2.get())[i];
			f2 = g.template get<pt::v>(sub2.get())[i];

			BOOST_REQUIRE_EQUAL(f1,f2);
		}

		++sub2;
	}

	// destroy the packed memory
	mem.decRef();
	delete &mem;

	//! [Unpack the object]

	}
}

BOOST_AUTO_TEST_SUITE_END()



#endif /* SRC_PACKER_UNIT_TESTS_HPP_ */
