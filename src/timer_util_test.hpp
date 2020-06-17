/*
 * timer_util_test.hpp
 *
 *  Created on: Jul 11, 2015
 *      Author: i-bird
 */

#ifndef TIMER_UTIL_TEST_HPP_
#define TIMER_UTIL_TEST_HPP_

#include "timer.hpp"

BOOST_AUTO_TEST_SUITE( timer_test )

BOOST_AUTO_TEST_CASE( timer_use )
{
	//! [timer usage and behavior]
	timer t;

	// start the timer
	t.start();

	sleep(1);

	// get the elapsed real time and cpu time without stop
	BOOST_REQUIRE_CLOSE(t.getwct(),1.0,20.0);
	BOOST_REQUIRE_SMALL(t.getcputime(),10.0);

	sleep(1);

	// stop the timer
	t.stop();

	sleep(1);

	// unusefull
	t.stop();

	// get the cpu time and real time
	t.getcputime();
	t.getwct();

	// get the elapsed real time and cpu time without stop
	BOOST_REQUIRE_CLOSE(t.getwct(),2.0,20.0);
	BOOST_REQUIRE_SMALL(t.getcputime(),10.0);

	t.reset();

	BOOST_REQUIRE_CLOSE(t.getwct(),0.0,20.0);
	BOOST_REQUIRE_SMALL(t.getcputime(),10.0);

	//! [timer usage and behavior]
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* TIMER_UTIL_TEST_HPP_ */
