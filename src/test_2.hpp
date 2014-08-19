#include <time.h>
#include <typeinfo>


template<typename g> void test_layout_grid3d(g & c3)
{
	  std::cout << "3D Array with grid_key (without redundant dimension): " << "\n";

	  typedef Point<float> P;

	   timespec ts_start;
	   // clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
	   clock_gettime(CLOCK_REALTIME, &ts_start); // Works on Linux

	  grid_key_dx<3> kk;

	  for (int i = 0 ; i < GS_SIZE ; i++)
	  {
	    for (int j = 0 ; j < GS_SIZE ; j++)
	    {
	      for (int k = 0 ; k < GS_SIZE ; k++)
	      {

		kk.set(i,j,k);

		c3.template get<P::x>(kk) = 1.1f;
		c3.template get<P::y>(kk) = 1.2f;
		c3.template get<P::z>(kk) = 1.3f;
		c3.template get<P::s>(kk) = 1.0f;

		c3.template get<P::v>(kk)[0] = 1.0f;
		c3.template get<P::v>(kk)[1] = 2.0f;
		c3.template get<P::v>(kk)[2] = 3.0f;

		c3.template get<P::t>(kk)[0][0] = 1.0f;
		c3.template get<P::t>(kk)[0][1] = 2.0f;
		c3.template get<P::t>(kk)[0][2] = 3.0f;
		c3.template get<P::t>(kk)[1][0] = 4.0f;
		c3.template get<P::t>(kk)[1][1] = 5.0f;
		c3.template get<P::t>(kk)[1][2] = 6.0f;
		c3.template get<P::t>(kk)[2][0] = 7.0f;
		c3.template get<P::t>(kk)[2][1] = 8.0f;
		c3.template get<P::t>(kk)[2][2] = 9.0f;

	      }
	    }
	   }

	  timespec end_time;
	  clock_gettime(CLOCK_REALTIME, &end_time); // Works on Linux
	   float time_dif =(float)( end_time.tv_sec - ts_start.tv_sec  + (double)(end_time.tv_nsec - ts_start.tv_nsec)/1000000000.0 );

	   std::cout << "End : " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024 << " MB " << "  Bandwidth: " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024/time_dif << " MB/s  ";

	   /////////////////////////////////// MEM CHECK ////////////////////////////////////////////////////////

	   bool passed = true;

	   for (int i = 0 ; i < GS_SIZE ; i++)
	   {
	     for (int j = 0 ; j < GS_SIZE ; j++)
	     {
	       for (int k = 0 ; k < GS_SIZE ; k++)
	       {
	    	   kk.set(i,j,k);

	    	   c3.template get<P::x>(kk) = i;
	    	   c3.template get<P::y>(kk) = j;
	    	   c3.template get<P::z>(kk) = k;
	    	   c3.template get<P::s>(kk) = i+j+k;

	    	   c3.template get<P::v>(kk)[0] = i;
	    	   c3.template get<P::v>(kk)[1] = j;
	    	   c3.template get<P::v>(kk)[2] = k;

	    	   c3.template get<P::t>(kk)[0][0] = i+i;
	    	   c3.template get<P::t>(kk)[0][1] = i+j;
	    	   c3.template get<P::t>(kk)[0][2] = i+k;
	    	   c3.template get<P::t>(kk)[1][0] = j+i;
	    	   c3.template get<P::t>(kk)[1][1] = j+j;
	    	   c3.template get<P::t>(kk)[1][2] = j+k;
	    	   c3.template get<P::t>(kk)[2][0] = k+i;
	    	   c3.template get<P::t>(kk)[2][1] = k+j;
	    	   c3.template get<P::t>(kk)[2][2] = k+k;
	       }
	     }
	   }

	   for (int i = 0 ; i < GS_SIZE ; i++)
	   {
	     for (int j = 0 ; j < GS_SIZE ; j++)
	     {
	       for (int k = 0 ; k < GS_SIZE ; k++)
	       {
	    	   kk.set(i,j,k);

	    	   if (c3.template get<P::x>(kk) != i) passed = false;
	    	   if (c3.template get<P::y>(kk) != j) passed = false;
	    	   if (c3.template get<P::z>(kk) != k) passed = false;
	    	   if (c3.template get<P::s>(kk) != i+j+k) passed = false;

	    	   if (c3.template get<P::v>(kk)[0] != i) passed = false;
	    	   if (c3.template get<P::v>(kk)[1] != j) passed = false;
	    	   if (c3.template get<P::v>(kk)[2] != k) passed = false;

	    	   if (c3.template get<P::t>(kk)[0][0] != i+i) passed = false;
	    	   if (c3.template get<P::t>(kk)[0][1] != i+j) passed = false;
	    	   if (c3.template get<P::t>(kk)[0][2] != i+k) passed = false;
	    	   if (c3.template get<P::t>(kk)[1][0] != j+i) passed = false;
	    	   if (c3.template get<P::t>(kk)[1][1] != j+j) passed = false;
	    	   if (c3.template get<P::t>(kk)[1][2] != j+k) passed = false;
	    	   if (c3.template get<P::t>(kk)[2][0] != k+i) passed = false;
	    	   if (c3.template get<P::t>(kk)[2][1] != k+j) passed = false;
	    	   if (c3.template get<P::t>(kk)[2][2] != k+k) passed = false;
	       }
	     }
	   }

	   if (passed == true)
		   std::cout << "PASSED"  << "\n";
	   else
		   std::cout << "FAILED"  << "\n";
}

template<unsigned int dim, typename g> void test_layout_gridNd(g & c3)
{
	  std::cout << dim << "D Array with grid_key (without redundant dimension): " << "\n";

	  typedef Point<float> P;

	   timespec ts_start;
	   // clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
	   clock_gettime(CLOCK_REALTIME, &ts_start); // Works on Linux

	  grid_key_dx<dim> kk;

	  grid_key_dx_iterator<dim> key_it;

	  while (key_it.hasNext())
	  {
		grid_key_dx<dim> & kk = key_it.next();

		c3.template get<P::x>(kk) = 1.1f;
		c3.template get<P::y>(kk) = 1.2f;
		c3.template get<P::z>(kk) = 1.3f;
		c3.template get<P::s>(kk) = 1.0f;

		c3.template get<P::v>(kk)[0] = 1.0f;
		c3.template get<P::v>(kk)[1] = 2.0f;
		c3.template get<P::v>(kk)[2] = 3.0f;

		c3.template get<P::t>(kk)[0][0] = 1.0f;
		c3.template get<P::t>(kk)[0][1] = 2.0f;
		c3.template get<P::t>(kk)[0][2] = 3.0f;
		c3.template get<P::t>(kk)[1][0] = 4.0f;
		c3.template get<P::t>(kk)[1][1] = 5.0f;
		c3.template get<P::t>(kk)[1][2] = 6.0f;
		c3.template get<P::t>(kk)[2][0] = 7.0f;
		c3.template get<P::t>(kk)[2][1] = 8.0f;
		c3.template get<P::t>(kk)[2][2] = 9.0f;

	  }

	  timespec end_time;
	  clock_gettime(CLOCK_REALTIME, &end_time); // Works on Linux
	   float time_dif =(float)( end_time.tv_sec - ts_start.tv_sec  + (double)(end_time.tv_nsec - ts_start.tv_nsec)/1000000000.0 );

	   std::cout << "End : " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024 << " MB " << "  Bandwidth: " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024/time_dif << " MB/s  ";

	   /////////////////////////////////// MEM CHECK ////////////////////////////////////////////////////////

	   bool passed = true;

	   for (int i = 0 ; i < GS_SIZE ; i++)
	   {
	     for (int j = 0 ; j < GS_SIZE ; j++)
	     {
	       for (int k = 0 ; k < GS_SIZE ; k++)
	       {
	    	   kk.set(i,j,k);

	    	   c3.template get<P::x>(kk) = i;
	    	   c3.template get<P::y>(kk) = j;
	    	   c3.template get<P::z>(kk) = k;
	    	   c3.template get<P::s>(kk) = i+j+k;

	    	   c3.template get<P::v>(kk)[0] = i;
	    	   c3.template get<P::v>(kk)[1] = j;
	    	   c3.template get<P::v>(kk)[2] = k;

	    	   c3.template get<P::t>(kk)[0][0] = i+i;
	    	   c3.template get<P::t>(kk)[0][1] = i+j;
	    	   c3.template get<P::t>(kk)[0][2] = i+k;
	    	   c3.template get<P::t>(kk)[1][0] = j+i;
	    	   c3.template get<P::t>(kk)[1][1] = j+j;
	    	   c3.template get<P::t>(kk)[1][2] = j+k;
	    	   c3.template get<P::t>(kk)[2][0] = k+i;
	    	   c3.template get<P::t>(kk)[2][1] = k+j;
	    	   c3.template get<P::t>(kk)[2][2] = k+k;
	       }
	     }
	   }

	   for (int i = 0 ; i < GS_SIZE ; i++)
	   {
	     for (int j = 0 ; j < GS_SIZE ; j++)
	     {
	       for (int k = 0 ; k < GS_SIZE ; k++)
	       {
	    	   kk.set(i,j,k);

	    	   if (c3.template get<P::x>(kk) != i) passed = false;
	    	   if (c3.template get<P::y>(kk) != j) passed = false;
	    	   if (c3.template get<P::z>(kk) != k) passed = false;
	    	   if (c3.template get<P::s>(kk) != i+j+k) passed = false;

	    	   if (c3.template get<P::v>(kk)[0] != i) passed = false;
	    	   if (c3.template get<P::v>(kk)[1] != j) passed = false;
	    	   if (c3.template get<P::v>(kk)[2] != k) passed = false;

	    	   if (c3.template get<P::t>(kk)[0][0] != i+i) passed = false;
	    	   if (c3.template get<P::t>(kk)[0][1] != i+j) passed = false;
	    	   if (c3.template get<P::t>(kk)[0][2] != i+k) passed = false;
	    	   if (c3.template get<P::t>(kk)[1][0] != j+i) passed = false;
	    	   if (c3.template get<P::t>(kk)[1][1] != j+j) passed = false;
	    	   if (c3.template get<P::t>(kk)[1][2] != j+k) passed = false;
	    	   if (c3.template get<P::t>(kk)[2][0] != k+i) passed = false;
	    	   if (c3.template get<P::t>(kk)[2][1] != k+j) passed = false;
	    	   if (c3.template get<P::t>(kk)[2][2] != k+k) passed = false;
	       }
	     }
	   }

	   if (passed == true)
		   std::cout << "PASSED"  << "\n";
	   else
		   std::cout << "FAILED"  << "\n";
}
