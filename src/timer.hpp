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

#include "util/cuda_util.hpp"

/*! \brief Class for cpu time benchmarking
 *
 * Usage:
 *
 * \snippet timer_util_test.hpp timer usage and behavior
 *
 */

class timer
{
	//! Flag that indicate if the timer is running or not
	bool running;

	//! starting time
    struct timespec tsstart;

    //! start time from epoch
    clock_t cstart;

    // stop time
    struct timespec tsstop;

    //! stop time from epoch
    clock_t cstop;

#if defined(__NVCC__) && !defined(CUDA_ON_CPU)
    #ifdef __HIP__
        hipEvent_t start_g, stop_g;
    #else
        cudaEvent_t start_g, stop_g;
    #endif
#endif

    // Fill the stop point
    void check()
    {
#if defined(SYNC_BEFORE_TAKE_TIME) && defined(__NVCC__)
        #ifdef __HIP__
        hipDeviceSynchronize();
        #else
        cudaDeviceSynchronize();
        #endif
#endif

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

public:

    //! Default constructor
    timer()
	:running(false),cstart(0),cstop(0)
    {
    	tsstart = timespec();
    	tsstop = timespec();
    }

    /*! \brief Start the timer
     *
     */
    void start()
    {
    	// time is running
    	running = true;

#if defined(SYNC_BEFORE_TAKE_TIME) && defined(__NVCC__) 
        cudaDeviceSynchronize();
#endif

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

    /*! \brief Stop the timer
     *
     *
     */
    void stop()
    {
    	if (running == false)	return;
    	running = false;
    	check();
    }

    /*! \brief Return the elapsed real time
     *
     *
     */
    double getwct()
    {
    	if (running == true)
    		check();

    	return ((double)(tsstop.tv_sec - tsstart.tv_sec)) + ((1e-9) * ((double)(tsstop.tv_nsec - tsstart.tv_nsec)));
    }

    /*! \brief Return the cpu time
     *
     *
     */
    double getcputime()
    {
    	if (running == true)
    		check();

        return (((double)(cstop - cstart)) / CLOCKS_PER_SEC);
    }

    /*! \brief Reset the timer
     *
     *
     */
    void reset()
    {
    	tsstart = tsstop;
    	cstart = cstop;
    }

    #if defined(__NVCC__) && !defined(CUDA_ON_CPU)

    void startGPU()
    {
        #ifdef __HIP__
        hipEventCreate(&start_g);
        hipEventRecord(start_g,0);
        #else
        cudaEventCreate(&start_g);
        cudaEventRecord(start_g,0);
        #endif
    }

    void stopGPU()
    {
        #ifdef __HIP__
        hipEventCreate(&stop_g);
        hipEventRecord(stop_g,0);
        hipEventSynchronize(stop_g);
        #else
        cudaEventCreate(&stop_g);
        cudaEventRecord(stop_g,0);
        cudaEventSynchronize(stop_g);
        #endif
    }

    double getwctGPU()
    {
        float elapsedTime;

        #ifdef __HIP__
        hipEventElapsedTime(&elapsedTime, start_g,stop_g);
        #else
        cudaEventElapsedTime(&elapsedTime, start_g,stop_g);
        #endif

        return elapsedTime;
    }

    #endif
};

#endif


