#include "map.hpp"
#include "memory_cpu.hpp"
#include "memory_gpu.hpp"
#include "Particle.hpp"
#include <boost/mpl/int.hpp>
#include <typeinfo>
#include <test_1.hpp>
#include <test_2.hpp>
//#include <test_3.hpp>
//#include <test_4.hpp>
//#include <test_5.hpp>
#include "memory_gpu_thrust.hpp"

/*float Wi(float dist_x)
{
  if (dist_x <= 1)
    return 1.5f*dist_x * dist_x * dist_x - 2.5f * dist_x * dist_x + 1.0f;
  else if (dist_x <= 2)
    return -0.5f*dist_x * dist_x * dist_x * dist_x + 2.5f * dist_x * dist_x - 4.0f * dist_x + 2.0f;
  else
    return 0.0;
}*/

int main()
{
/*  tensor<int,3,3,3> c;
  tensor<tensor<int,3,3,3>,3,3,3> c2;*/
  
  std::vector<size_t> sz;
  sz.push_back(GS_SIZE);
  sz.push_back(GS_SIZE);
  sz.push_back(GS_SIZE);
  
  layout_gpu< grid<Point<float>>, memory_gpu<memory_gpu_type<Point<float>>::type> > c3(sz);

  // cpu test
  
  test1();
  
//  layout_cpu< Particles<Point<float>, memory_cpu<float> >, particle_key > p1;
   
  test2(c3);

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
