/*
 * zmorton_unit_tests.hpp
 *
 *  Created on: Aug 1, 2019
 *      Author: i-bird
 */

#ifndef ZMORTON_UNIT_TESTS_HPP_
#define ZMORTON_UNIT_TESTS_HPP_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "util/zmorton.hpp"

template<typename T>
bool check(size_t res, grid_key_dx<2,T> & k)
{
	bool check = true;

	check &= (k.get(0) & 0x1) == (res & 0x1);
	check &= (k.get(0) & (0x1ul << 1)) == ((res & (0x1ul << 2)) >> 1);
	check &= (k.get(0) & (0x1ul << 2)) == ((res & (0x1ul << 4)) >> 2);
	check &= (k.get(0) & (0x1ul << 3)) == ((res & (0x1ul << 6)) >> 3);
	check &= (k.get(0) & (0x1ul << 4)) == ((res & (0x1ul << 8)) >> 4);
	check &= (k.get(0) & (0x1ul << 5)) == ((res & (0x1ul << 10)) >> 5);
	check &= (k.get(0) & (0x1ul << 6)) == ((res & (0x1ul << 12)) >> 6);
	check &= (k.get(0) & (0x1ul << 7)) == ((res & (0x1ul << 14)) >> 7);
	check &= (k.get(0) & (0x1ul << 8)) == ((res & (0x1ul << 16)) >> 8);
	check &= (k.get(0) & (0x1ul << 9)) == ((res & (0x1ul << 19)) >> 9);
	check &= (k.get(0) & (0x1ul << 10)) == ((res & (0x1ul << 20)) >> 10);
	check &= (k.get(0) & (0x1ul << 11)) == ((res & (0x1ul << 22)) >> 11);
	check &= (k.get(0) & (0x1ul << 12)) == ((res & (0x1ul << 24)) >> 12);
	check &= (k.get(0) & (0x1ul << 13)) == ((res & (0x1ul << 26)) >> 13);
	check &= (k.get(0) & (0x1ul << 14)) == ((res & (0x1ul << 28)) >> 14);
	check &= (k.get(0) & (0x1ul << 15)) == ((res & (0x1ul << 30)) >> 15);
	check &= (k.get(0) & (0x1ul << 16)) == ((res & (0x1ul << 32)) >> 16);
	check &= (k.get(0) & (0x1ul << 17)) == ((res & (0x1ul << 34)) >> 17);
	check &= (k.get(0) & (0x1ul << 18)) == ((res & (0x1ul << 36)) >> 19);
	check &= (k.get(0) & (0x1ul << 19)) == ((res & (0x1ul << 38)) >> 19);
	check &= (k.get(0) & (0x1ul << 20)) == ((res & (0x1ul << 40)) >> 20);
	check &= (k.get(0) & (0x1ul << 21)) == ((res & (0x1ul << 42)) >> 21);
	check &= (k.get(0) & (0x1ul << 22)) == ((res & (0x1ul << 44)) >> 22);
	check &= (k.get(0) & (0x1ul << 23)) == ((res & (0x1ul << 46)) >> 23);
	check &= (k.get(0) & (0x1ul << 24)) == ((res & (0x1ul << 48)) >> 24);
	check &= (k.get(0) & (0x1ul << 25)) == ((res & (0x1ul << 50)) >> 25);
	check &= (k.get(0) & (0x1ul << 26)) == ((res & (0x1ul << 52)) >> 26);
	check &= (k.get(0) & (0x1ul << 27)) == ((res & (0x1ul << 54)) >> 27);
	check &= (k.get(0) & (0x1ul << 28)) == ((res & (0x1ul << 56)) >> 28);
	check &= (k.get(0) & (0x1ul << 29)) == ((res & (0x1ul << 58)) >> 29);
	check &= (k.get(0) & (0x1ul << 30)) == ((res & (0x1ul << 60)) >> 30);
	check &= (k.get(0) & (0x1ul << 31)) == ((res & (0x1ul << 62)) >> 31);

	res = res >> 1;

	check &= (k.get(1) & 0x1) == (res & 0x1);
	check &= (k.get(1) & (0x1ul << 1)) == ((res & (0x1ul << 2)) >> 1);
	check &= (k.get(1) & (0x1ul << 2)) == ((res & (0x1ul << 4)) >> 2);
	check &= (k.get(1) & (0x1ul << 3)) == ((res & (0x1ul << 6)) >> 3);
	check &= (k.get(1) & (0x1ul << 4)) == ((res & (0x1ul << 8)) >> 4);
	check &= (k.get(1) & (0x1ul << 5)) == ((res & (0x1ul << 10)) >> 5);
	check &= (k.get(1) & (0x1ul << 6)) == ((res & (0x1ul << 12)) >> 6);
	check &= (k.get(1) & (0x1ul << 7)) == ((res & (0x1ul << 14)) >> 7);
	check &= (k.get(1) & (0x1ul << 8)) == ((res & (0x1ul << 16)) >> 8);
	check &= (k.get(1) & (0x1ul << 9)) == ((res & (0x1ul << 19)) >> 9);
	check &= (k.get(1) & (0x1ul << 10)) == ((res & (0x1ul << 20)) >> 10);
	check &= (k.get(1) & (0x1ul << 11)) == ((res & (0x1ul << 22)) >> 11);
	check &= (k.get(1) & (0x1ul << 12)) == ((res & (0x1ul << 24)) >> 12);
	check &= (k.get(1) & (0x1ul << 13)) == ((res & (0x1ul << 26)) >> 13);
	check &= (k.get(1) & (0x1ul << 14)) == ((res & (0x1ul << 28)) >> 14);
	check &= (k.get(1) & (0x1ul << 15)) == ((res & (0x1ul << 30)) >> 15);
	check &= (k.get(1) & (0x1ul << 16)) == ((res & (0x1ul << 32)) >> 16);
	check &= (k.get(1) & (0x1ul << 17)) == ((res & (0x1ul << 34)) >> 17);
	check &= (k.get(1) & (0x1ul << 18)) == ((res & (0x1ul << 36)) >> 19);
	check &= (k.get(1) & (0x1ul << 19)) == ((res & (0x1ul << 38)) >> 19);
	check &= (k.get(1) & (0x1ul << 20)) == ((res & (0x1ul << 40)) >> 20);
	check &= (k.get(1) & (0x1ul << 21)) == ((res & (0x1ul << 42)) >> 21);
	check &= (k.get(1) & (0x1ul << 22)) == ((res & (0x1ul << 44)) >> 22);
	check &= (k.get(1) & (0x1ul << 23)) == ((res & (0x1ul << 46)) >> 23);
	check &= (k.get(1) & (0x1ul << 24)) == ((res & (0x1ul << 48)) >> 24);
	check &= (k.get(1) & (0x1ul << 25)) == ((res & (0x1ul << 50)) >> 25);
	check &= (k.get(1) & (0x1ul << 26)) == ((res & (0x1ul << 52)) >> 26);
	check &= (k.get(1) & (0x1ul << 27)) == ((res & (0x1ul << 54)) >> 27);
	check &= (k.get(1) & (0x1ul << 28)) == ((res & (0x1ul << 56)) >> 28);
	check &= (k.get(1) & (0x1ul << 29)) == ((res & (0x1ul << 58)) >> 29);
	check &= (k.get(1) & (0x1ul << 30)) == ((res & (0x1ul << 60)) >> 30);
	check &= (k.get(1) & (0x1ul << 31)) == ((res & (0x1ul << 62)) >> 31);

	return check;
}

BOOST_AUTO_TEST_SUITE( zmorton_suite_test )

BOOST_AUTO_TEST_CASE( zmorton_linearization_test )
{
	{
		grid_key_dx<1> key(0);

		BOOST_REQUIRE_EQUAL(lin_zid(key),0);

		grid_key_dx<1> key2(100);

		BOOST_REQUIRE_EQUAL(lin_zid(key2),100);
	}

	{
		grid_key_dx<2> key({0,0});

		BOOST_REQUIRE_EQUAL(lin_zid(key),0);

		grid_key_dx<2> key2({2,2});

		BOOST_REQUIRE_EQUAL(lin_zid(key2),12);

		grid_key_dx<2> key3({3,2});

		BOOST_REQUIRE_EQUAL(lin_zid(key3),13);

		grid_key_dx<2> key4({0,3});

		BOOST_REQUIRE_EQUAL(lin_zid(key4),10);

		grid_key_dx<2> key5({3,0});

		BOOST_REQUIRE_EQUAL(lin_zid(key5),5);

		grid_key_dx<2> key6({165,347});

		size_t res = lin_zid(key6);

		BOOST_REQUIRE_EQUAL(check(res,key6),true);

		grid_key_dx<2> key7({0xF,0XF});

		res = lin_zid(key7);

		BOOST_REQUIRE_EQUAL(res,0xFF);

		grid_key_dx<2> key8({0xFF,0XFF});

		res = lin_zid(key8);

		BOOST_REQUIRE_EQUAL(res,0xFFFF);

		grid_key_dx<2> key9({0xFFFF,0XFFFF});

		res = lin_zid(key9);

		BOOST_REQUIRE_EQUAL(res,0xFFFFFFFF);

		grid_key_dx<2> key10({0xFFFFFFFF,0XFFFFFFFF});

		res = lin_zid(key10);

		BOOST_REQUIRE_EQUAL(res,0xFFFFFFFFFFFFFFFF);
	}

	{
		grid_key_dx<3> key({0,0,0});

		BOOST_REQUIRE_EQUAL(lin_zid(key),0);

		grid_key_dx<3> key2({2,2,2});

		BOOST_REQUIRE_EQUAL(lin_zid(key2),56);

		grid_key_dx<3> key3({4,4,4});

		BOOST_REQUIRE_EQUAL(lin_zid(key3),448);

		grid_key_dx<3> key4({8,8,8});

		BOOST_REQUIRE_EQUAL(lin_zid(key4),3584);

		grid_key_dx<3> key5({16,16,16});

		BOOST_REQUIRE_EQUAL(lin_zid(key5),(1 << 12) + (1 << 13) + (1 << 14));

		grid_key_dx<3> key6({32,32,32});

		BOOST_REQUIRE_EQUAL(lin_zid(key6),(1 << 15) + (1 << 16) + (1 << 17));

		grid_key_dx<3> key7({64,64,64});

		BOOST_REQUIRE_EQUAL(lin_zid(key7),(1 << 18) + (1 << 19) + (1 << 20));

		grid_key_dx<3> key8({128,128,128});

		BOOST_REQUIRE_EQUAL(lin_zid(key8),(1 << 21) + (1 << 22) + (1 << 23));

		grid_key_dx<3> key9({256,256,256});

		BOOST_REQUIRE_EQUAL(lin_zid(key9),(1 << 24) + (1 << 25) + (1 << 26));

		grid_key_dx<3> key10({512,512,512});

		BOOST_REQUIRE_EQUAL(lin_zid(key10),(1 << 27) + (1 << 28) + (1 << 29));

		grid_key_dx<3> key11({1024,1024,1024});

		BOOST_REQUIRE_EQUAL(lin_zid(key11),(1ul << 30) + (1ul << 31) + (1ul << 32));

		grid_key_dx<3> key12({2048,2048,2048});

		BOOST_REQUIRE_EQUAL(lin_zid(key12),(1ul << 33) + (1ul << 34) + (1ul << 35));

		grid_key_dx<3> key13({4096,4096,4096});

		BOOST_REQUIRE_EQUAL(lin_zid(key13),(1ul << 36) + (1ul << 37) + (1ul << 38));

		grid_key_dx<3> key14({8192,8192,8192});

		BOOST_REQUIRE_EQUAL(lin_zid(key14),(1ul << 39) + (1ul << 40) + (1ul << 41));

		grid_key_dx<3> key15({16384,16384,16384});

		BOOST_REQUIRE_EQUAL(lin_zid(key15),(1ul << 42) + (1ul << 43) + (1ul << 44));

		grid_key_dx<3> key16({32768,32768,32768});

		BOOST_REQUIRE_EQUAL(lin_zid(key16),(1ul << 45) + (1ul << 46) + (1ul << 47));

		grid_key_dx<3> key17({65536,65536,65536});

		BOOST_REQUIRE_EQUAL(lin_zid(key17),(1ul << 48) + (1ul << 49) + (1ul << 50));

		grid_key_dx<3> key18({131072,131072,131072});

		BOOST_REQUIRE_EQUAL(lin_zid(key18),(1ul << 51) + (1ul << 52) + (1ul << 53));

		grid_key_dx<3> key19({262144,262144,262144});

		BOOST_REQUIRE_EQUAL(lin_zid(key19),(1ul << 54) + (1ul << 55) + (1ul << 56));

		grid_key_dx<3> key20({524288,524288,524288});

		BOOST_REQUIRE_EQUAL(lin_zid(key20),(1ul << 57) + (1ul << 58) + (1ul << 59));

		grid_key_dx<3> key21({1048576,1048576,1048576});

		BOOST_REQUIRE_EQUAL(lin_zid(key21),(1ul << 60) + (1ul << 61) + (1ul << 62));
	}
}


BOOST_AUTO_TEST_CASE( zmorton_invlinearization_test )
{
	{
		grid_key_dx<1> key(0);
		grid_key_dx<1> ikey;

		size_t lin = lin_zid(key);
		invlin_zid(lin,ikey);

		BOOST_REQUIRE(key == ikey);

		grid_key_dx<1> key2(100);
		grid_key_dx<1> ikey2;

		lin = lin_zid(key2);
		invlin_zid(lin,ikey2);

		BOOST_REQUIRE(key2 == ikey2);
	}

	{
		grid_key_dx<2> key({0,0});
		grid_key_dx<2> ikey;

		size_t lin = lin_zid(key);
		invlin_zid(lin,ikey);

		BOOST_REQUIRE(key == ikey);

		grid_key_dx<2> key2({2,2});
		grid_key_dx<2> ikey2;

		lin = lin_zid(key2);
		invlin_zid(lin,ikey2);

		BOOST_REQUIRE(key2 == ikey2);

		grid_key_dx<2> key3({3,2});
		grid_key_dx<2> ikey3;

		lin = lin_zid(key3);
		invlin_zid(lin,ikey3);

		BOOST_REQUIRE(key3 == ikey3);

		grid_key_dx<2> key4({0,3});
		grid_key_dx<2> ikey4;

		lin = lin_zid(key4);
		invlin_zid(lin,ikey4);

		BOOST_REQUIRE(key4 == ikey4);

		grid_key_dx<2> key5({3,0});
		grid_key_dx<2> ikey5({3,0});

		lin = lin_zid(key5);
		invlin_zid(lin,ikey5);

		BOOST_REQUIRE(key5 == ikey5);

		grid_key_dx<2> key6({165,347});
		grid_key_dx<2> ikey6;

		lin = lin_zid(key6);
		invlin_zid(lin,ikey6);

		BOOST_REQUIRE(key6 == ikey6);

		grid_key_dx<2> key7({0xF,0XF});
		grid_key_dx<2> ikey7;

		lin = lin_zid(key7);
		invlin_zid(lin,ikey7);

		BOOST_REQUIRE(key7 == ikey7);

		grid_key_dx<2> key8({0xFF,0XFF});
		grid_key_dx<2> ikey8;

		lin = lin_zid(key8);
		invlin_zid(lin,ikey8);

		BOOST_REQUIRE(key8 == ikey8);

		grid_key_dx<2> key9({0xFFFF,0XFFFF});
		grid_key_dx<2> ikey9;

		lin = lin_zid(key9);
		invlin_zid(lin,ikey9);

		BOOST_REQUIRE(key9 == ikey9);

		grid_key_dx<2> key10({0xFFFFFFFF,0XFFFFFFFF});
		grid_key_dx<2> ikey10;

		lin = lin_zid(key10);
		invlin_zid(lin,ikey10);

		BOOST_REQUIRE(key10 == ikey10);
	}

	{
		grid_key_dx<3> key({0,0,0});
		grid_key_dx<3> ikey;

		size_t lin = lin_zid(key);
		invlin_zid(lin,ikey);

		BOOST_REQUIRE(key == ikey);

		grid_key_dx<3> key2({2,2,2});
		grid_key_dx<3> ikey2;

		lin = lin_zid(key2);
		invlin_zid(lin,ikey2);

		BOOST_REQUIRE(key2 == ikey2);

		grid_key_dx<3> key3({4,4,4});
		grid_key_dx<3> ikey3;

		lin = lin_zid(key3);
		invlin_zid(lin,ikey3);

		BOOST_REQUIRE(key3 != ikey);

		grid_key_dx<3> key4({8,8,8});
		grid_key_dx<3> ikey4;

		lin = lin_zid(key4);
		invlin_zid(lin,ikey4);

		BOOST_REQUIRE(key4 == ikey4);

		grid_key_dx<3> key5({16,16,16});
		grid_key_dx<3> ikey5;

		lin = lin_zid(key5);
		invlin_zid(lin,ikey5);

		BOOST_REQUIRE(key5 == ikey5);

		grid_key_dx<3> key6({32,32,32});
		grid_key_dx<3> ikey6;

		lin = lin_zid(key6);
		invlin_zid(lin,ikey6);

		BOOST_REQUIRE(key6 == ikey6);

		grid_key_dx<3> key7({64,64,64});
		grid_key_dx<3> ikey7;

		lin = lin_zid(key7);
		invlin_zid(lin,ikey7);

		BOOST_REQUIRE(key7 == ikey7);

		grid_key_dx<3> key8({128,128,128});
		grid_key_dx<3> ikey8;

		lin = lin_zid(key8);
		invlin_zid(lin,ikey8);

		BOOST_REQUIRE(key8 == ikey8);

		grid_key_dx<3> key9({256,256,256});
		grid_key_dx<3> ikey9;

		lin = lin_zid(key9);
		invlin_zid(lin,ikey9);

		BOOST_REQUIRE(key9 == ikey9);

		grid_key_dx<3> key10({512,512,512});
		grid_key_dx<3> ikey10;

		lin = lin_zid(key10);
		invlin_zid(lin,ikey10);

		BOOST_REQUIRE(key10 == ikey10);

		grid_key_dx<3> key11({1024,1024,1024});
		grid_key_dx<3> ikey11;

		lin = lin_zid(key11);
		invlin_zid(lin,ikey11);

		BOOST_REQUIRE(key11 == ikey11);

		grid_key_dx<3> key12({2048,2048,2048});
		grid_key_dx<3> ikey12;

		lin = lin_zid(key12);
		invlin_zid(lin,ikey12);

		BOOST_REQUIRE(key12 == ikey12);

		grid_key_dx<3> key13({4096,4096,4096});
		grid_key_dx<3> ikey13;

		lin = lin_zid(key13);
		invlin_zid(lin,ikey13);

		BOOST_REQUIRE(key13 == ikey13);

		grid_key_dx<3> key14({8192,8192,8192});
		grid_key_dx<3> ikey14;

		lin = lin_zid(key14);
		invlin_zid(lin,ikey14);

		BOOST_REQUIRE(key14 == ikey14);

		grid_key_dx<3> key15({16384,16384,16384});
		grid_key_dx<3> ikey15;

		lin = lin_zid(key15);
		invlin_zid(lin,ikey15);

		BOOST_REQUIRE(key15 == ikey15);

		grid_key_dx<3> key16({32768,32768,32768});
		grid_key_dx<3> ikey16;

		lin = lin_zid(key16);
		invlin_zid(lin,ikey16);

		BOOST_REQUIRE(key16 == ikey16);

		grid_key_dx<3> key17({65536,65536,65536});
		grid_key_dx<3> ikey17;

		lin = lin_zid(key17);
		invlin_zid(lin,ikey17);

		BOOST_REQUIRE(key17 == ikey17);

		grid_key_dx<3> key18({131072,131072,131072});
		grid_key_dx<3> ikey18;

		lin = lin_zid(key18);
		invlin_zid(lin,ikey18);

		BOOST_REQUIRE(key18 == ikey18);

		grid_key_dx<3> key19({262144,262144,262144});
		grid_key_dx<3> ikey19;

		lin = lin_zid(key19);
		invlin_zid(lin,ikey19);

		BOOST_REQUIRE(key19 == ikey19);

		grid_key_dx<3> key20({524288,524288,524288});
		grid_key_dx<3> ikey20;

		lin = lin_zid(key20);
		invlin_zid(lin,ikey20);

		BOOST_REQUIRE(key20 == ikey20);

		grid_key_dx<3> key21({1048576,1048576,1048576});
		grid_key_dx<3> ikey21;

		lin = lin_zid(key21);
		invlin_zid(lin,ikey21);

		BOOST_REQUIRE(key21 == ikey21);

		bool match = true;

		for (unsigned int i = 0 ; i < 100 ; i++)
		{
			for (unsigned int j = 0 ; j < 100 ; j++)
			{
				for (unsigned int k = 0 ; k < 100 ; k++)
				{
					grid_key_dx<3> key23({i,j,k});
					grid_key_dx<3> ikey23;

					lin = lin_zid(key23);
					invlin_zid(lin,ikey23);

					match &= key23 == ikey23;
				}
			}
		}

		BOOST_REQUIRE(match == true);
	}
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* ZMORTON_UNIT_TESTS_HPP_ */
