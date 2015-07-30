#define IS_GRID 2
#define IS_VECTOR 1
#define GRAPH
// dim dimensionalility
// T grid or vector
// St type of space float double ....

template <unsigned int dim, typename T, typename St, int impl = 2*is_grid<T>::value + is_vector<T>::value >
class weight_counter
{
    size_t weight(T & t, const Box<dim,St> & domain, Box<dim,St> & b)
    {     
      //error if T is not a grid or a vector
          std::cerr << "Not implemented";
        return 0;
    }
}

/////// When T = openfpm::grid_cpu

template <unsigned int dim, typename T, typename St>
class weight_counter<dim,T,St,IS_GRID>
{
    size_t weight(T & t, Box<dim,T> & b)
    {
       // Grid implementation
      // get a size vector
      size_t a[dim] = t.getGrid().getSize();
      size_t res = 1;
      //for each dimension
      for (int i = 0; i < dim; i++)
          size_t count = 0;
          //rounding coordinates of Low and High up and down
          int a1 = (int)(b.getLow(i));
          int a2 = (int)(b.getHigh(i) + 0.9999999);
          //counting points those get into a box, excluding borders
          for (int j = a1+1; j < a2; j++)
            count++;
          //adding negative borders, if the number was integer
          if (a1 == ceil(a1))
            count++;
          //multypling by all dimensions
          res *= count;
      return res;
    }
}

////// When T = openfpm::vector

template <unsigned int dim, typename T, typename St>
class weight_counter<dim,T,St,IS_VECTOR>
{
    size_t weight(openfpm::vector<Point<dim>> & t, const Box<dim,St> & domain, Box<dim,St> & b)
    {
       // Vector implementation
      //get size of a vector = number of points
      size_t size = t.size();
      int count = 0;
      //count points those are inside of a box, excluding positive borders
      for (int i = 0; i < size, i++) 
      {
        if (b.isInsidePE(t.get(i)))
          count++;
      }
      return count;
    }
}
