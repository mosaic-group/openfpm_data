/*
 * test_5.hpp
 *
 *  Created on: Feb 4, 2014
 *      Author: i-bird
 */

void test5(layout_gpu< grid<Point<float>>, memory_gpu<memory_gpu_type<Point<float>>::type> > & c3)
{
	   std::cout << "3D Array with grid_key runtime dimensions " << "\n";

	   grid_key<Point<float>::x> kx(3);
	   grid_key<Point<float>::y> ky(kx);
	   grid_key<Point<float>::z> kz(kx);
	   grid_key<Point<float>::s> ks(kx);
	   grid_key<Point<float>::v> kv(kx);
	   grid_key<Point<float>::t> kt(kx);

	   timespec ts_start;
	   // clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
	   clock_gettime(CLOCK_REALTIME, &ts_start); // Works on Linux

	   for (int i = 0 ; i < GS_SIZE ; i++)
	   {
	    for (int j = 0 ; j < GS_SIZE ; j++)
	    {
	      for (int k = 0 ; k < GS_SIZE ; k++)
	      {
	    	  kx.set(i,j,k);

	    	  c3.get(kx) = 1.1f;
	    	  c3.get(ky) = 1.2f;
	    	  c3.get(kz) = 1.3f;
	    	  c3.get(ks) = 1.0f;

	    	  c3.get(kv)[0] = 1.0f;
	    	  c3.get(kv)[0] = 2.0f;
	    	  c3.get(kv)[0] = 3.0f;

	    	  c3.get(kt)[0][0] = 1.0f;
	    	  c3.get(kt)[0][1] = 2.0f;
	    	  c3.get(kt)[0][2] = 3.0f;
	    	  c3.get(kt)[1][0] = 4.0f;
	    	  c3.get(kt)[1][1] = 5.0f;
	    	  c3.get(kt)[1][2] = 6.0f;
	    	  c3.get(kt)[2][0] = 7.0f;
	    	  c3.get(kt)[2][1] = 8.0f;
	    	  c3.get(kt)[2][2] = 9.0f;
	      }
	    }
	   }

	   timespec end_time;
	   clock_gettime(CLOCK_REALTIME, &end_time); // Works on Linux
	    float time_dif =(float)( end_time.tv_sec - ts_start.tv_sec  + (double)(end_time.tv_nsec - ts_start.tv_nsec)/1000000000.0 );

	   std::cout << "End : " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024 << " MB " << "  Bandwidth: " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024/time_dif << " MB/s  \n";
}

void test5(layout_cpu< grid<Point<float>>, memory_cpu<memory_cpu_type<Point<float>>::type> > & c3)
{
   std::cout << "3D Array with grid_key runtime dimensions " << "\n";

   grid_key<Point<float>::x> kx(3);
   grid_key<Point<float>::y> ky(kx);
   grid_key<Point<float>::z> kz(kx);
   grid_key<Point<float>::s> ks(kx);
   grid_key<Point<float>::v> kv(kx);
   grid_key<Point<float>::t> kt(kx);

   timespec ts_start;
   // clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
   clock_gettime(CLOCK_REALTIME, &ts_start); // Works on Linux

   for (int i = 0 ; i < GS_SIZE ; i++)
   {
    for (int j = 0 ; j < GS_SIZE ; j++)
    {
      for (int k = 0 ; k < GS_SIZE ; k++)
      {
	kx.set(i,j,k);

	c3.get(kx) = 1.1f;
	c3.get(ky) = 1.2f;
	c3.get(kz) = 1.3f;
	c3.get(ks) = 1.0f;

	c3.get(kv)[0] = 1.0f;
	c3.get(kv)[0] = 2.0f;
	c3.get(kv)[0] = 3.0f;

	c3.get(kt)[0][0] = 1.0f;
	c3.get(kt)[0][1] = 2.0f;
	c3.get(kt)[0][2] = 3.0f;
	c3.get(kt)[1][0] = 4.0f;
	c3.get(kt)[1][1] = 5.0f;
	c3.get(kt)[1][2] = 6.0f;
	c3.get(kt)[2][0] = 7.0f;
	c3.get(kt)[2][1] = 8.0f;
	c3.get(kt)[2][2] = 9.0f;
      }
    }
   }

   timespec end_time;
   clock_gettime(CLOCK_REALTIME, &end_time); // Works on Linux
    float time_dif =(float)( end_time.tv_sec - ts_start.tv_sec  + (double)(end_time.tv_nsec - ts_start.tv_nsec)/1000000000.0 );

   std::cout << "End : " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024 << " MB " << "  Bandwidth: " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024/time_dif << " MB/s  \n";
}
