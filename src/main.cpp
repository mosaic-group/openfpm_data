#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <boost/mpl/int.hpp>
#include <typeinfo>
#include <ct_array.hpp>
#include "memory/CudaMemory.cuh"
#include "memory/HeapMemory.hpp"
#include "memory_conf.hpp"
#include "map_grid.hpp"
#include "map_vector.hpp"

// Include tests

#include <test_1.hpp>
#include <test_2.hpp>

/*! \brief Test all grid with dimensionality dim and size sz on all dimensions
 *
 * Test all grid with dimensionality dim and size sz on all dimensions
 *
 */

template<unsigned int dim> void test_all_grid(size_t sz)
{
	std::vector<size_t> szz;
	szz.clear();

	for (int i = 0 ; i < dim ; i++)
	{szz.push_back(sz);}

	{grid_cpu<dim, Point<float> > c3(szz);
	c3.template setMemory<CudaMemory>();
	test_layout_gridNd<dim>(c3,sz);}

	{grid_cpu<dim, Point<float> > c3(szz);
	c3.template setMemory<HeapMemory>();
	test_layout_gridNd<dim>(c3,sz);}

	{grid_gpu<dim, Point<float> > c3(szz);
	c3.template setMemory<CudaMemory>();
	test_layout_gridNd<dim>(c3,sz);}

	{grid_gpu<dim, Point<float> > c3(szz);
	c3.template setMemory<HeapMemory>();
	test_layout_gridNd<dim>(c3,sz);}
}


int main()
{
/*  tensor<int,3,3,3> c;
  tensor<tensor<int,3,3,3>,3,3,3> c2;*/
  
  std::vector<size_t> sz;
  sz.push_back(GS_SIZE);
  sz.push_back(GS_SIZE);
  sz.push_back(GS_SIZE);

  // This is an ordinary test simple 3D with plain C array
  
  test1();

  // test the grid from dimensionality 1 to 8 with several size non multiple of two

  // Dimension 8-1

  test_all_grid<8>(4);
  test_all_grid<7>(8);
  test_all_grid<6>(9);
  test_all_grid<5>(18);
  test_all_grid<4>(37);
  test_all_grid<3>(126);
  test_all_grid<2>(1414);
  test_all_grid<1>(2000000);

   // Test the 3d gpu grid with Cudamemory and HeapMemory with different size

  for (int i = 2 ; i <= GS_SIZE ; i++)
  {
	  sz.clear();
	  sz.push_back(i);
	  sz.push_back(i);
	  sz.push_back(i);

	  {grid_gpu<3, Point<float> > c3(sz);
	  c3.setMemory<CudaMemory>();
	  test_layout_grid3d(c3,i);}

	  {grid_gpu<3, Point<float> > c3(sz);
	  c3.setMemory<HeapMemory>();
	  test_layout_grid3d(c3,i);}

	  // Test the 3d cpu grid with Cudamemory and HeapMemory

	  {grid_cpu<3, Point<float> > c3(sz);
	  c3.setMemory<CudaMemory>();
	  test_layout_grid3d(c3,i);}

	  {grid_cpu<3, Point<float> > c3(sz);
	  c3.setMemory<HeapMemory>();
	  test_layout_grid3d(c3,i);}

  }

  // Test openfpm vector

  openfpm::vector<Point<float>> ofv;

  // Test another grid

//  test3(c3);
   
//  test4(c3);

//  test5(c3);


//   k.set(2,2,2);
//   c3.get(k).x = 6.0;
   
//   std::cout << c3.get(k) << "\n";
    
 //  c3.output("vtk");
   
   // iterator test
   
//   auto it1 = c3.iterator();
   
/*   while (it1.next())
   {
     it1.get().x = it1.getId(0);
     it1.get().y = it1.getId(1);
     it1.get().z = it1.getId(2);
   }*/
   
   // iterator test
   
/*   auto pit = p1.iterator();
   
   while (pit.next())
   {
     pit.get().x = rand() % 10;
     pit.get().y = rand() % 10;
     pit.get().z = rand() % 10;
   }
   
   c3.create_neighborhood(p1);
   
   // interpolation particle-cell
   
   auto git = c3.iterator();
   while (git.next())
   {
     n_git = git.get().neighboorhood().iterator();
     while(n_git.next())
     {
       float dist = abs(n_git.get().x - git.get().x);
       float Wi_t *= Wi(dist);
       
       dist = abs(n_git.get().y - git.get().y);
       Wi_t *= Wi(dist);
              
       float dist = abs(n_git.get().z - git.get().z);
       Wi_t *= Wi(dist);
       
       git.get().strenght += Wi_t * n_git.get().strenght;
     }
   }*/
   
   // interpolation cell-particle
   
   
   
   // GPU test
   
   
   
   // Xeon Phi test
   
   
   
//  layout_cpu<grid< tensor<int,3,3,3>, memory_cpu<float> >, grid_key > c4;
  
//  grid< tensor< tensor<int,3,3,3> ,3,3,3> > c5(sz);  

//  std::cout << c2.size() << "\n";
//  std::cout << c5.size() << "\n";

  

}
