#define GS_SIZE 256


void test1()
{
	// Testing the grids

	std::cout << "Grid size known at runtime" << "\n";
	std::cout << "1D Array with index calculation: " << "\n";

	Point_orig<float> * pA = new Point_orig<float>[GS_SIZE*GS_SIZE*GS_SIZE];

	int gs_sq = GS_SIZE*GS_SIZE;
	int gs = GS_SIZE;

	clock_t begin_time = clock();

	for (int i = 0 ; i < GS_SIZE ; i++)
	{
		for (int j = 0 ; j < GS_SIZE ; j++)
		{
	      for (int k = 0 ; k < GS_SIZE ; k++)
	      {
	    	  pA[i*gs_sq+j*gs+k].x = 1.1f;
	    	  pA[i*gs_sq+j*gs+k].y = 1.2f;
	    	  pA[i*gs_sq+j*gs+k].z = 1.3f;
	    	  pA[i*gs_sq+j*gs+k].s = 1.0f;

	    	  pA[i*gs_sq+j*gs+k].v[0] = 1.0f;
	    	  pA[i*gs_sq+j*gs+k].v[1] = 2.0f;
	    	  pA[i*gs_sq+j*gs+k].v[2] = 3.0f;

	    	  pA[i*gs_sq+j*gs+k].t[0][0] = 1.0f;
	    	  pA[i*gs_sq+j*gs+k].t[0][1] = 2.0f;
	    	  pA[i*gs_sq+j*gs+k].t[0][2] = 3.0f;
	    	  pA[i*gs_sq+j*gs+k].t[1][0] = 4.0f;
	    	  pA[i*gs_sq+j*gs+k].t[1][1] = 5.0f;
	    	  pA[i*gs_sq+j*gs+k].t[1][2] = 6.0f;
	    	  pA[i*gs_sq+j*gs+k].t[2][0] = 7.0f;
	    	  pA[i*gs_sq+j*gs+k].t[2][1] = 8.0f;
	    	  pA[i*gs_sq+j*gs+k].t[2][2] = 9.0f;
	      }
	    }
	   }

	   clock_t end_time = clock ();
	   float time_dif =(float)( end_time - begin_time ) / (float) CLOCKS_PER_SEC;

	   std::cout << "End : " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024 << " MB " << "  Bandwidth: " << GS_SIZE*GS_SIZE*GS_SIZE*16*4/1024/1024/time_dif << " MB/s  \n";
}


