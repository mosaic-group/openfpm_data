/*
 * unit_test_init_cleanup.hpp
 *
 *  Created on: Apr 17, 2015
 *      Author: Pietro Incardona
 */

#ifndef UNIT_TEST_INIT_CLEANUP_HPP_
#define UNIT_TEST_INIT_CLEANUP_HPP_

#include "util/cudify/cudify.hpp"

//! boost unit test fixation (start procedure to call before testing)
struct ut_start
{

	//! start procedure before running the test
    ut_start()
    {
#ifdef CUDA_ON_CPU
        init_wrappers();
#endif
#ifdef PERFORMANCE_TEST
    	test_dir = getenv("OPENFPM_PERFORMANCE_TEST_DIR");

    	if (test_dir == NULL)
    	{
    		std::cerr << "Error: " __FILE__ << ":" << __LINE__ << " in order to run the performance test you must set the environment variable $OPENFPM_PERFORMANCE_TEST_DIR to the test or an empty directory";
    		exit(1);
    	}
#endif

    }

    //! post procedure to call after the test
    ~ut_start()
    {
    }
};

//____________________________________________________________________________//

BOOST_GLOBAL_FIXTURE( ut_start );



#endif /* UNIT_TEST_INIT_CLEANUP_HPP_ */
