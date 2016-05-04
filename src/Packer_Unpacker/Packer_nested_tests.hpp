/*!
 * This file contains the tests for vector/vector_std nested packer and unpacker
 * Created on: Oct 11, 2015
 *     Author: Yaroslav Zaluzhnyi
 */


#ifndef SRC_PACKER_NESTED_TESTS_HPP_
#define SRC_PACKER_NESTED_TESTS_HPP_

#include "Pack_selector.hpp"
#include "Packer.hpp"
#include "Unpacker.hpp"
#include "Grid/grid_util_test.hpp"
#include <iostream>
#include "data_type/aggregate.hpp"

struct test_box_vpack
{
	size_t i;
	size_t j;
	Box<3,float> bx;
};

//Testing packing and unpacking for different vectors

BOOST_AUTO_TEST_SUITE( packer_unpacker )

BOOST_AUTO_TEST_CASE ( vector_ptst_packer_unpacker )
{
	std::cout << "Vector pack/unpack test start" << "\n";

	openfpm::vector<openfpm::vector<openfpm::vector<Point_test<float>>>> v;
	for (size_t i = 0; i < 5; i++) {
		openfpm::vector<openfpm::vector<Point_test<float>>> v4;
		for (size_t j = 0; j < 6; j++) {
			v4.add(allocate_openfpm<openfpm::vector<Point_test<float>>>(7));
		}
		v.add(v4);
	}

	typedef Point_test<float> pt;

	//Pack request vector
	std::vector<size_t> pap_prp;
	
	//Pack requesting
	
	Packer<decltype(v),HeapMemory>::packRequest<pt::x, pt::v>(v,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],(sizeof(float) + sizeof(float[3])) * 7);

	//Just to see the elements of pack request vector
#ifdef DEBUG
	for (size_t i = 0; i < pap_prp.size(); i++)
		std::cout << pap_prp[i] << std::endl;
#endif

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	//Packing
	Pack_stat sts;

	Packer<decltype(v),HeapMemory>::pack<pt::x,pt::v>(mem,v,sts);

	//Unpacking
	Unpack_stat ps;

	openfpm::vector<openfpm::vector<openfpm::vector<Point_test<float>>>> v_unp;

 	Unpacker<decltype(v_unp),HeapMemory>::unpack<pt::x,pt::v>(mem,v_unp,ps);

 	//Check equality vector v_unp
	for (size_t k = 0; k < v_unp.size(); k++)
	{
		for (size_t i = 0; i < v_unp.get(k).size(); i++)
		{
			auto it = v_unp.get(k).get(i).getIterator();

			while (it.isNext())
			{
				float f1 = v_unp.get(k).get(i).template get<pt::x>(it.get());
				float f2 = v.get(k).get(i).template get<pt::x>(it.get());

				BOOST_REQUIRE_EQUAL(f1,f2);

				for (size_t j = 0 ; j < 3 ; j++)
				{
					f1 = v_unp.get(k).get(i).template get<pt::v>(it.get())[j];
					f2 = v.get(k).get(i).template get<pt::v>(it.get())[j];

					BOOST_REQUIRE_EQUAL(f1,f2);
				}
				++it;
			}
		}
	}

	mem.decRef();
	delete &mem;
}

BOOST_AUTO_TEST_CASE ( vector_std_packer_unpacker )
{
	openfpm::vector<openfpm::vector<openfpm::vector<float>>> v2;
	for (size_t i = 0; i < 5; i++) {
		openfpm::vector<openfpm::vector<float>> v6;
		for (size_t j = 0; j < 6; j++) {
			openfpm::vector<float> v7;
			for (size_t k = 0; k < 7; k++) {
				v7.add(1);
			}
			v6.add(v7);
		}
		v2.add(v6);
	}

	//Pack request vector
	std::vector<size_t> pap_prp;

	//Pack requesting

	Packer<decltype(v2),HeapMemory>::packRequest<>(v2,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(float) * 7);

	//Just to see the elements of pack request vector
#ifdef DEBUG
	for (size_t i = 0; i < pap_prp.size(); i++)
		std::cout << pap_prp[i] << std::endl;
#endif

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	//Packing

	Pack_stat sts;

	Packer<decltype(v2),HeapMemory>::pack<>(mem,v2,sts);

	//Unpacking

	Unpack_stat ps;

	openfpm::vector<openfpm::vector<openfpm::vector<float>>> v2_unp;

 	Unpacker<decltype(v2_unp),HeapMemory>::unpack<>(mem,v2_unp,ps);

	//Check equality vector v2_unp
	for (size_t k = 0; k < v2_unp.size(); k++)
	{
		for (size_t i = 0; i < v2_unp.get(k).size(); i++)
		{
			for (size_t j = 0; j < v2_unp.get(k).get(i).size(); j++)
			{
				float f1 = v2_unp.get(k).get(i).get(j);
				float f2 = v2.get(k).get(i).get(j);

				BOOST_REQUIRE_EQUAL(f1,f2);
			}
		}
	}

	mem.decRef();
	delete &mem;
}

BOOST_AUTO_TEST_CASE ( vector_zerosize_packer_unpacker )
{
	openfpm::vector<openfpm::vector<openfpm::vector<Point_test<float>>>> v5;
	for (size_t i = 0; i < 5; i++) {
		openfpm::vector<openfpm::vector<Point_test<float>>> v4;
		for (size_t j = 0; j < 6; j++) {
			v4.add(allocate_openfpm<openfpm::vector<Point_test<float>>>(7));
		}
		v5.add(v4);
	}
	openfpm::vector<openfpm::vector<Point_test<float>>> v51;
	for (size_t j = 0; j < 5; j++) {
		v51.add(allocate_openfpm<openfpm::vector<Point_test<float>>>(7));
	}
	v51.add(allocate_openfpm<openfpm::vector<Point_test<float>>>(0));

	v5.add(v51);

	typedef Point_test<float> pt;

	//Pack request vector
	std::vector<size_t> pap_prp;

	//Pack requesting
	Packer<decltype(v5),HeapMemory>::packRequest<pt::x, pt::v>(v5,pap_prp);
	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],0ul);

	//Just to see the elements of pack request vector
#ifdef DEBUG
	for (size_t i = 0; i < pap_prp.size(); i++)
		std::cout << pap_prp[i] << std::endl;
#endif

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	//Packing
	Pack_stat sts;

	Packer<decltype(v5),HeapMemory>::pack<pt::x,pt::v>(mem,v5,sts);

	//Unpacking
	Unpack_stat ps;

	openfpm::vector<openfpm::vector<openfpm::vector<Point_test<float>>>> v5_unp;

 	Unpacker<decltype(v5_unp),HeapMemory>::unpack<pt::x,pt::v>(mem,v5_unp,ps);

 	//Check equality vector v5_unp
	for (size_t k = 0; k < v5_unp.size(); k++)
	{
		for (size_t i = 0; i < v5_unp.get(k).size(); i++)
		{
			auto it = v5_unp.get(k).get(i).getIterator();

			while (it.isNext())
			{
				float f1 = v5_unp.get(k).get(i).template get<pt::x>(it.get());
				float f2 = v5.get(k).get(i).template get<pt::x>(it.get());

				BOOST_REQUIRE_EQUAL(f1,f2);

				for (size_t j = 0 ; j < 3 ; j++)
				{
					f1 = v5_unp.get(k).get(i).template get<pt::v>(it.get())[j];
					f2 = v5.get(k).get(i).template get<pt::v>(it.get())[j];

					BOOST_REQUIRE_EQUAL(f1,f2);
				}
				++it;
			}
		}
	}

	mem.decRef();
	delete &mem;
}

BOOST_AUTO_TEST_CASE ( vector_zerosize__lvl_2_packer_unpacker )
{
	openfpm::vector<openfpm::vector<openfpm::vector<Point_test<float>>>> v5;
	for (size_t i = 0; i < 5; i++) {
		openfpm::vector<openfpm::vector<Point_test<float>>> v4;
		for (size_t j = 0; j < 6; j++) {
			v4.add(allocate_openfpm<openfpm::vector<Point_test<float>>>(7));
		}
		v5.add(v4);
	}
	openfpm::vector<openfpm::vector<Point_test<float>>> v51;
	v51.clear();
	v5.add(v51);

	v5.get(2).get(1).clear();

	typedef Point_test<float> pt;

	//Pack request vector
	std::vector<size_t> pap_prp;

	//Pack requesting
	Packer<decltype(v5),HeapMemory>::packRequest<pt::x, pt::v>(v5,pap_prp);

	//Just to see the elements of pack request vector
#ifdef DEBUG
	for (size_t i = 0; i < pap_prp.size(); i++)
		std::cout << pap_prp[i] << std::endl;
#endif

	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(size_t));

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	//Packing
	Pack_stat sts;

	Packer<decltype(v5),HeapMemory>::pack<pt::x,pt::v>(mem,v5,sts);

	//Unpacking
	Unpack_stat ps;

	openfpm::vector<openfpm::vector<openfpm::vector<Point_test<float>>>> v5_unp;

 	Unpacker<decltype(v5_unp),HeapMemory>::unpack<pt::x,pt::v>(mem,v5_unp,ps);

 	//Check equality vector v5_unp
 	BOOST_REQUIRE_EQUAL(v5_unp.size(),v5.size());
	for (size_t k = 0; k < v5_unp.size(); k++)
	{
		BOOST_REQUIRE_EQUAL(v5_unp.get(k).size(),v5.get(k).size());
		for (size_t i = 0; i < v5_unp.get(k).size(); i++)
		{
			BOOST_REQUIRE_EQUAL(v5_unp.get(k).get(i).size(),v5.get(k).get(i).size());
			auto it = v5_unp.get(k).get(i).getIterator();

			while (it.isNext())
			{
				float f1 = v5_unp.get(k).get(i).template get<pt::x>(it.get());
				float f2 = v5.get(k).get(i).template get<pt::x>(it.get());

				BOOST_REQUIRE_EQUAL(f1,f2);

				for (size_t j = 0 ; j < 3 ; j++)
				{
					f1 = v5_unp.get(k).get(i).template get<pt::v>(it.get())[j];
					f2 = v5.get(k).get(i).template get<pt::v>(it.get())[j];

					BOOST_REQUIRE_EQUAL(f1,f2);
				}
				++it;
			}
		}
	}

	mem.decRef();
	delete &mem;
}


BOOST_AUTO_TEST_CASE ( vector_zerosize__lvl_2_packer_unpacker_float )
{
	openfpm::vector<openfpm::vector<openfpm::vector<float>>> v2;
	for (size_t i = 0; i < 5; i++) {
		openfpm::vector<openfpm::vector<float>> v6;
		for (size_t j = 0; j < 6; j++) {
			openfpm::vector<float> v7;
			for (size_t k = 0; k < 7; k++) {
				v7.add(1);
			}
			v6.add(v7);
		}
		v2.add(v6);
	}

	v2.get(0).clear();

	//Pack request vector
	std::vector<size_t> pap_prp;

	//Pack requesting

	Packer<decltype(v2),HeapMemory>::packRequest<>(v2,pap_prp);
//	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],sizeof(float) * 7);

	//Just to see the elements of pack request vector
#ifdef DEBUG
	for (size_t i = 0; i < pap_prp.size(); i++)
		std::cout << pap_prp[i] << std::endl;
#endif

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	//Packing

	Pack_stat sts;

	Packer<decltype(v2),HeapMemory>::pack<>(mem,v2,sts);

	//Unpacking

	Unpack_stat ps;

	openfpm::vector<openfpm::vector<openfpm::vector<float>>> v2_unp;

 	Unpacker<decltype(v2_unp),HeapMemory>::unpack<>(mem,v2_unp,ps);

	//Check equality vector v2_unp
	for (size_t k = 0; k < v2_unp.size(); k++)
	{
		for (size_t i = 0; i < v2_unp.get(k).size(); i++)
		{
			for (size_t j = 0; j < v2_unp.get(k).get(i).size(); j++)
			{
				float f1 = v2_unp.get(k).get(i).get(j);
				float f2 = v2.get(k).get(i).get(j);

				BOOST_REQUIRE_EQUAL(f1,f2);
			}
		}
	}

	mem.decRef();
	delete &mem;
}

BOOST_AUTO_TEST_CASE ( vector_std_smarter_packer_unpacker )
{
	//Create an object

	openfpm::vector<openfpm::vector<test_box_vpack>> v;

	//Fill it with data

	v.resize(4);
	for (size_t i = 0; i < v.size(); i++)
	{
		v.get(i).resize(5);
		for (size_t j = 0; j < v.get(i).size(); j++)
		{
			v.get(i).get(j).i = 1;
			v.get(i).get(j).j = 2;
			v.get(i).get(j).bx.setHigh(0, 3.0);
			v.get(i).get(j).bx.setHigh(1, 4.0);
			v.get(i).get(j).bx.setHigh(2, 5.0);
			v.get(i).get(j).bx.setLow(0, 6.0);
			v.get(i).get(j).bx.setLow(1, 7.0);
			v.get(i).get(j).bx.setLow(2, 8.0);
		}
	}

	//Pack_request vector
	std::vector<size_t> pap_prp;

	typedef Point_test<float> pt;
	//Pack request
	Packer<decltype(v),HeapMemory>::packRequest<>(v,pap_prp);

	//Just to see the elements of pack request vector
#ifdef DEBUG
	for (size_t i = 0; i < pap_prp.size(); i++)
		std::cout << pap_prp[i] << std::endl;
#endif

	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],(sizeof(float)*6 + sizeof(size_t)*2) * 5);

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	//Packing

	Pack_stat sts;

	Packer<decltype(v),HeapMemory>::pack<>(mem,v,sts);

	//Unpacking

	Unpack_stat ps;

	openfpm::vector<openfpm::vector<test_box_vpack>> v_unp;

	Unpacker<decltype(v_unp),HeapMemory>::unpack<>(mem,v_unp,ps);

	for (size_t i = 0; i < v.size(); i++)
	{
		for (size_t j = 0; j < v.get(i).size(); j++)
		{
			float s1 = v_unp.get(i).get(j).i;
			float s2 = v.get(i).get(j).i;
			float s3 = v_unp.get(i).get(j).j;
			float s4 = v.get(i).get(j).j;

			BOOST_REQUIRE_EQUAL(s1,s2);
			BOOST_REQUIRE_EQUAL(s3,s4);

			Box<3,float> b1 = v_unp.get(i).get(j).bx.getBox();
			Box<3,float> b2 = v.get(i).get(j).bx.getBox();
			BOOST_REQUIRE(b1 == b2);
		}
	}

	mem.decRef();
	delete &mem;
}

BOOST_AUTO_TEST_CASE ( vector_smarter_packer_unpacker )
{
	//Create an object
	openfpm::vector<openfpm::vector<aggregate<float,float,openfpm::vector<Point_test<float>>>>> v4;

	//Fill it with data

	v4.resize(1);
	Point_test<float> p;
	p.fill();
	for (size_t i = 0 ; i < v4.size() ; i++)
	{
		v4.get(i).resize(1);
		for (size_t j = 0 ; j < v4.get(i).size() ; j++)
		{
			v4.get(i).template get<0>(j) = 1.0;
			v4.get(i).template get<1>(j) = 2.0;
			v4.get(i).template get<2>(j).resize(2);

			for (size_t k = 0 ; k < v4.get(i).template get<2>(j).size() ; k++)
				v4.get(i).template get<2>(j).get(k) = p;
		}
	}

	//Pack_request vector
	std::vector<size_t> pap_prp;

	//Pack request
	Packer<decltype(v4),HeapMemory>::packRequest<1,2>(v4,pap_prp);

	//Just to see the elements of pack request vector
#ifdef DEBUG
	for (size_t i = 0; i < pap_prp.size(); i++)
		std::cout << pap_prp[i] << std::endl;
#endif

	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],((sizeof(float)*4 + sizeof(float[3])) + sizeof(float[3][3]))*2);

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	//Packing

	Pack_stat sts;

	Packer<decltype(v4),HeapMemory>::pack<1,2>(mem,v4,sts);

	//Unpacking

	Unpack_stat ps;

	openfpm::vector<openfpm::vector<aggregate<float,float,openfpm::vector<Point_test<float>>>>> v4_unp;

	Unpacker<decltype(v4_unp),HeapMemory>::unpack<1,2>(mem,v4_unp,ps);

	//Check the data

	for (size_t i = 0 ; i < v4.size() ; i++)
	{
		for (size_t j = 0 ; j < v4.get(i).size() ; j++)
		{
			//float f1 = v4_unp.get(i).template get<0>(j);
			//float f2 = v4.get(i).template get<0>(j);
			float f3 = v4_unp.get(i).template get<1>(j);
			float f4 = v4.get(i).template get<1>(j);

			//BOOST_REQUIRE_EQUAL(f1,f2);
			BOOST_REQUIRE_EQUAL(f3,f4);

			for (size_t k = 0 ; k < v4.get(i).template get<2>(j).size() ; k++)
			{
				Point_test<float> p1 = v4_unp.get(i).template get<2>(j).get(k);
				Point_test<float> p2 = v4.get(i).template get<2>(j).get(k);

				BOOST_REQUIRE(p1 == p2);
			}
		}
	}

	mem.decRef();
	delete &mem;
}

BOOST_AUTO_TEST_CASE ( vector_smarter_packer_unpacker_2 )
{
	//Create an object
	openfpm::vector<openfpm::vector<aggregate<float,float,Point_test<float>,Point_test<float>,Point_test<float>,Point_test<float>,Point_test<float>,Point_test<float>,Point_test<float>>>> v4;

	//Fill it with data
	v4.resize(1);
	Point_test<float> p;
	p.fill();
	for (size_t i = 0 ; i < v4.size() ; i++)
	{
		v4.get(i).resize(5000);
		for (size_t j = 0 ; j < v4.get(i).size() ; j++)
		{
			v4.get(i).template get<0>(j) = 1.0;
			v4.get(i).template get<1>(j) = 2.0;
			v4.get(i).template get<2>(j) = p;
			v4.get(i).template get<3>(j) = p;
			v4.get(i).template get<4>(j) = p;
			v4.get(i).template get<5>(j) = p;
			v4.get(i).template get<6>(j) = p;
			v4.get(i).template get<7>(j) = p;
			v4.get(i).template get<8>(j) = p;
		}
	}

	//Pack_request vector
	std::vector<size_t> pap_prp;

	//Pack request
	Packer<decltype(v4),HeapMemory>::packRequest<0,1,2,3,4,5,6,7,8>(v4,pap_prp);

	//Just to see the elements of pack request vector
#ifdef DEBUG
	for (size_t i = 0; i < pap_prp.size(); i++)
		std::cout << pap_prp[i] << std::endl;
#endif

	BOOST_REQUIRE_EQUAL(pap_prp[pap_prp.size()-1],((sizeof(float)*30 + sizeof(float[3])*7) + sizeof(float[3][3])*7)*5000);

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	//Packing

	Pack_stat sts;

	Packer<decltype(v4),HeapMemory>::pack<0,1,2,3,4,5,6,7,8>(mem,v4,sts);

	//Unpacking

	Unpack_stat ps;

	openfpm::vector<openfpm::vector<aggregate<float,float,Point_test<float>,Point_test<float>,Point_test<float>,Point_test<float>,Point_test<float>,Point_test<float>,Point_test<float>>>> v4_unp;

	Unpacker<decltype(v4_unp),HeapMemory>::unpack<0,1,2,3,4,5,6,7,8>(mem,v4_unp,ps);

	//Checking the data
	for (size_t i = 0; i < v4.size(); i++)
	{
		for (size_t j = 0; j < v4.get(i).size(); j++)
		{
			float f1 = v4_unp.get(i).template get<0>(j);
			float f2 = v4.get(i).template get<0>(j);
			float f3 = v4_unp.get(i).template get<1>(j);
			float f4 = v4.get(i).template get<1>(j);

			Point_test<float> p1 = v4_unp.get(i).template get<2>(j);
			Point_test<float> p2 = v4.get(i).template get<2>(j);

			BOOST_REQUIRE_EQUAL(f1,f2);
			BOOST_REQUIRE_EQUAL(f3,f4);
			BOOST_REQUIRE(p1 == p2);
		}
	}

	std::cout << "Vector pack/unpack test stop" << "\n";

	mem.decRef();
	delete &mem;
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_PACKER_NESTED_TESTS_HPP_ */
