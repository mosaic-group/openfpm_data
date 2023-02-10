/*
 * unit_test_init_cleanup.hpp
 *
 *  Created on: May 1, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_IO_SRC_UNIT_TEST_INIT_CLEANUP_IO_HPP_
#define OPENFPM_IO_SRC_UNIT_TEST_INIT_CLEANUP_IO_HPP_


struct ut_start
{
	//!
    ut_start()
    {
    	BOOST_TEST_MESSAGE("Initialize global VCluster");

    	openfpm_init(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);
    }

    ~ut_start()
    {
    	BOOST_TEST_MESSAGE("Delete global VClster");
    	openfpm_finalize();
    }
};

//____________________________________________________________________________//

BOOST_GLOBAL_FIXTURE( ut_start );


#endif /* OPENFPM_IO_SRC_UNIT_TEST_INIT_CLEANUP_IO_HPP_ */
