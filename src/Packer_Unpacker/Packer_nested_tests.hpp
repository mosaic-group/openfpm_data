#ifndef SRC_PACKER_NESTED_TESTS_HPP_
#define SRC_PACKER_NESTED_TESTS_HPP_

#include "Pack_selector.hpp"
#include "Packer.hpp"
#include "Unpacker.hpp"
#include "Grid/grid_util_test.hpp"
#include <iostream>

//Testing packing and unpacking for different vectors

BOOST_AUTO_TEST_SUITE( packer_unpacker )

BOOST_AUTO_TEST_CASE ( vector_ptst_packer_unpacker )
{
	std::cout << "Vector pack/unpack test start" << "\n";

	openfpm::vector<openfpm::vector<openfpm::vector<Point_test<float>>>> v;
	for (size_t i = 0; i < 5; i++) {
		openfpm::vector<openfpm::vector<Point_test<float>>> v4;
		for (size_t j = 0; j < 6; j++) {
			v4.add(allocate_openfpm(7));
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

	std::cout << "Vector pack/unpack test stop" << "\n";

}
BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_PACKER_NESTED_TESTS_HPP_ */
