//
//  Timer.h
//
//

#ifndef TIMER_HPP
#define TIMER_HPP

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

/*! \brief Class for cpu time benchmarking
 *
 * Usage:
 *
 * timet t;
 * t.start();
 *
 * ... Do something
 *
 * t.stop();
 * t.getwct();
 *
 */

class timer
{
    struct timespec tsstart;
    clock_t cstart;

    struct timespec tsstop;
    clock_t cstop;

public:
    timer()
    {
    }

    void start()
    {
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
        clock_serv_t cclock;
        mach_timespec_t mts;
        host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
        clock_get_time(cclock, &mts);
        mach_port_deallocate(mach_task_self(), cclock);
        tsstart.tv_sec = mts.tv_sec;
        tsstart.tv_nsec = mts.tv_nsec;
#else
        clock_gettime(CLOCK_REALTIME, &tsstart);
#endif
        cstart = clock();

    }

    void stop()
    {
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
        clock_serv_t cclock;
        mach_timespec_t mts;
        host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
        clock_get_time(cclock, &mts);
        mach_port_deallocate(mach_task_self(), cclock);
        tsstop.tv_sec = mts.tv_sec;
        tsstop.tv_nsec = mts.tv_nsec;
#else
        clock_gettime(CLOCK_REALTIME, &tsstop);
#endif
        cstop = clock();
    }

    double getwct()
    {
        return ((double)(tsstop.tv_sec - tsstart.tv_sec)) + ((1e-9) * ((double)(tsstop.tv_nsec - tsstart.tv_nsec)));
    }

    double getcputime()
    {
        return (((double)(cstop - cstart)) / CLOCKS_PER_SEC);
    }

};

#endif


