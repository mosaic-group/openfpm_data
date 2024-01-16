#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "multi_array_ref_openfpm.hpp"
#include <boost/mpl/vector.hpp>
#include <iostream>

BOOST_AUTO_TEST_SUITE( multi_array_ref_openfpm_suite )

BOOST_AUTO_TEST_CASE( multi_array_ref_openfpm_use )
{
	std::cout << "Test multi array start" << "\n";

	{

	float test_mem[10][3];
	float * p_test = &test_mem[0][0];

	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			p_test[j*10+i] = i*10 + j;  // <--------- Note I am inverting from the natural i*10+j
		}
	}

	openfpm::multi_array_ref_openfpm<float,2,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>>> ar((float *)test_mem,10);

	BOOST_REQUIRE_EQUAL(ar[0][0],0);
	BOOST_REQUIRE_EQUAL(ar[0][1],1);
	BOOST_REQUIRE_EQUAL(ar[1][0],10);
	BOOST_REQUIRE_EQUAL(ar[2][2],22);
	BOOST_REQUIRE_EQUAL(ar[1][2],12);
	BOOST_REQUIRE_EQUAL(ar[9][2],92);

	}

	{

	float test_mem[10][3][7];
	float * p_test = &test_mem[0][0][0];

	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			for (size_t k = 0 ; k < 7 ; k++)
			{
				p_test[(j*7+k)*10+i] = i*100 + j*10 + k;  // <--------- Note I am inverting from the natural i*10+j
			}
		}
	}

	openfpm::multi_array_ref_openfpm<float,3,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>,boost::mpl::int_<7>>> ar2((float *)test_mem,10);

	BOOST_REQUIRE_EQUAL(ar2[0][0][0],0);
	BOOST_REQUIRE_EQUAL(ar2[0][0][1],1);
	BOOST_REQUIRE_EQUAL(ar2[0][0][2],2);
	BOOST_REQUIRE_EQUAL(ar2[0][0][3],3);
	BOOST_REQUIRE_EQUAL(ar2[0][0][4],4);
	BOOST_REQUIRE_EQUAL(ar2[0][0][6],6);
	BOOST_REQUIRE_EQUAL(ar2[0][1][0],10);
	BOOST_REQUIRE_EQUAL(ar2[1][0][4],104);
	BOOST_REQUIRE_EQUAL(ar2[2][2][6],226);
	BOOST_REQUIRE_EQUAL(ar2[1][2][6],126);
	BOOST_REQUIRE_EQUAL(ar2[9][2][3],923);
	}


	{

	float test_mem[10][3][7][2];
	float * p_test = &test_mem[0][0][0][0];

	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			for (size_t k = 0 ; k < 7 ; k++)
			{
				for (size_t s = 0 ; s < 2 ; s++)
				{
					p_test[(j*7*2+k*2+s)*10+i] = i*1000 + j*100 + k*10 + s;  // <--------- Note I am inverting from the natural
				}
			}
		}
	}

	openfpm::multi_array_ref_openfpm<float,4,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>,boost::mpl::int_<7>,boost::mpl::int_<2>>> ar2((float *)test_mem,10);

	BOOST_REQUIRE_EQUAL(ar2[0][0][0][0],0);
	BOOST_REQUIRE_EQUAL(ar2[0][0][0][1],1);
	BOOST_REQUIRE_EQUAL(ar2[0][0][2][0],20);
	BOOST_REQUIRE_EQUAL(ar2[0][0][2][1],21);
	BOOST_REQUIRE_EQUAL(ar2[0][0][3][0],30);
	BOOST_REQUIRE_EQUAL(ar2[0][0][4][0],40);
	BOOST_REQUIRE_EQUAL(ar2[0][2][0][1],201);
	BOOST_REQUIRE_EQUAL(ar2[0][1][0][0],100);
	BOOST_REQUIRE_EQUAL(ar2[0][1][6][1],161);
	BOOST_REQUIRE_EQUAL(ar2[2][2][6][1],2261);
	BOOST_REQUIRE_EQUAL(ar2[1][1][5][1],1151);
	BOOST_REQUIRE_EQUAL(ar2[9][2][4][1],9241);

	// write_test

	ar2[1][1][5][1] = 1111;
	BOOST_REQUIRE_EQUAL(ar2[1][1][5][1],1111);

	}

	std::cout << "End multi array stop" << "\n";
}

BOOST_AUTO_TEST_CASE( multi_array_ref_openfpm_copy )
{
	std::cout << "Test multi array copy start" << "\n";

	{

	float test_mem[10][3][7][2];
	float * p_test = &test_mem[0][0][0][0];

	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			for (size_t k = 0 ; k < 7 ; k++)
			{
				for (size_t s = 0 ; s < 2 ; s++)
				{
					p_test[(j*7*2+k*2+s)*10+i] = i*1000 + j*100 + k*10 + s;  // <--------- Note I am inverting from the natural
				}
			}
		}
	}

	openfpm::multi_array_ref_openfpm<float,4,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>,boost::mpl::int_<7>,boost::mpl::int_<2>>> ar2((float *)test_mem,10);

	float test_mem2[10][3][7][2];
	openfpm::multi_array_ref_openfpm<float,4,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>,boost::mpl::int_<7>,boost::mpl::int_<2>>> ar3((float *)test_mem2,10);

	// Copy and check

	ar3 = ar2;

	bool check = true;
	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			for (size_t k = 0 ; k < 7 ; k++)
			{
				for (size_t s = 0 ; s < 2 ; s++)
				{
					check &= ar3[i][j][k][s] == ar2[i][j][k][s];
				}
			}
		}
	}

	BOOST_REQUIRE_EQUAL(check,true);

	}

	std::cout << "End multi array copy stop" << "\n";
}

BOOST_AUTO_TEST_CASE( multi_array_ref_openfpm_bind_ref )
{
	std::cout << "Test multi array bind_ref start" << "\n";

	{

	float test_mem[10][3][7][2];
	float * p_test = &test_mem[0][0][0][0];

	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			for (size_t k = 0 ; k < 7 ; k++)
			{
				for (size_t s = 0 ; s < 2 ; s++)
				{
					p_test[(j*7*2+k*2+s)*10+i] = i*1000 + j*100 + k*10 + s;  // <--------- Note I am inverting from the natural
				}
			}
		}
	}

	openfpm::multi_array_ref_openfpm<float,4,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>,boost::mpl::int_<7>,boost::mpl::int_<2>>> ar2((float *)test_mem,10);

	openfpm::multi_array_ref_openfpm<float,4,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>,boost::mpl::int_<7>,boost::mpl::int_<2>>> ar3((float *)NULL,10);
	ar3.bind_ref(ar2);

	// Copy and check

	bool check = true;
	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			for (size_t k = 0 ; k < 7 ; k++)
			{
				for (size_t s = 0 ; s < 2 ; s++)
				{
					check &= ar3[i][j][k][s] == ar2[i][j][k][s];
				}
			}
		}
	}

	BOOST_REQUIRE_EQUAL(check,true);

	}

	std::cout << "End multi array copy bind_ref" << "\n";
}

BOOST_AUTO_TEST_CASE( multi_array_ref_openfpm_swap )
{
	std::cout << "Test multi array swap start" << "\n";

	{

	float test_mem[10][3][7][2];
	float * p_test = &test_mem[0][0][0][0];

	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			for (size_t k = 0 ; k < 7 ; k++)
			{
				for (size_t s = 0 ; s < 2 ; s++)
				{
					p_test[(j*7*2+k*2+s)*10+i] = i*1000 + j*100 + k*10 + s;  // <--------- Note I am inverting from the natural
				}
			}
		}
	}

	openfpm::multi_array_ref_openfpm<float,4,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>,boost::mpl::int_<7>,boost::mpl::int_<2>>> ar2((float *)test_mem,10);

	float test_mem2[10][3][7][2];
	float * p_test2 = &test_mem2[0][0][0][0];

	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			for (size_t k = 0 ; k < 7 ; k++)
			{
				for (size_t s = 0 ; s < 2 ; s++)
				{
					p_test2[(j*7*2+k*2+s)*10+i] = i*1000 + j*100 + k*10 + s + 100000;  // <--------- Note I am inverting from the natural
				}
			}
		}
	}

	openfpm::multi_array_ref_openfpm<float,4,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>,boost::mpl::int_<7>,boost::mpl::int_<2>>> ar3((float *)test_mem2,10);
	ar3.swap(ar2);

	// Copy and check

	bool check = true;
	for (size_t i = 0 ; i < 10 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{
			for (size_t k = 0 ; k < 7 ; k++)
			{
				for (size_t s = 0 ; s < 2 ; s++)
				{
					check &= ar3[i][j][k][s] == i*1000 + j*100 + k*10 + s;
					check &= ar2[i][j][k][s] == i*1000 + j*100 + k*10 + s + 100000;
				}
			}
		}
	}

	BOOST_REQUIRE_EQUAL(check,true);

	}

	std::cout << "End multi array swap stop" << "\n";
}

BOOST_AUTO_TEST_CASE( test_is_multi_array )
{
	bool test = openfpm::is_multi_array<int>::value;

	BOOST_REQUIRE_EQUAL(test,false);

	test = openfpm::is_multi_array<openfpm::multi_array_ref_openfpm<float,4,boost::mpl::vector<boost::mpl::int_<-1>,boost::mpl::int_<3>,boost::mpl::int_<7>,boost::mpl::int_<2>>>>::value;

	BOOST_REQUIRE_EQUAL(test,true);
}

BOOST_AUTO_TEST_SUITE_END()


