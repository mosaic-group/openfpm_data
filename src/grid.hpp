#include <iostream>
#include <boost/shared_array.hpp>
#include <vector>
#include <initializer_list>
#include <array>
#include "memory.hpp"

#define HARDWARE 1

/*! \brief grid_key_dx is the key to access any element in the grid
 *
 * grid_key_dx is the key to access any element in the grid
 *
 * \param dim dimensionality of the grid
 *
 */

template<unsigned int dim>
class grid_key_dx
{
public:
  
  //! Constructor
  grid_key_dx()
  {}
  
  //! Construct a grid key from a list of numbers
  template<typename a, typename ...T>grid_key_dx(a v,T...t)
  {
    k[dim-1] = v;
    invert_assign(t...);
  }
  
  //! set the grid key from a list of numbers
  template<typename a, typename ...T>void set(a v, T...t)
  {
    k[dim-1] = v;
    invert_assign(t...);
  }
  
  /*! \brief get the i index
   *
   * Get the i index
   *
   * \param i index to get
   *
   * \return the index value
   *
   */
  mem_id get(size_t i)
  {
	  return k[i];
  }
  
  /*! \brief Set the i index
   *
   * Set the i index
   *
   * \param i index to set
   * \param id value to set
   *
   */
  void set_d(size_t i, mem_id id)
  {
	  k[i] = id;
  }

  //! structure that store all the index
  mem_id k[dim];

private:

  template<typename a, typename ...T>void invert_assign(a v,T...t)
  {
    k[sizeof...(T)] = v;
    invert_assign(t...);
  }

  template<typename a, typename ...T>void invert_assign(a v)
  {
    k[0] = v;
  }

};


/*! \brief grid_key_d is the key to access any element in the grid
 *
 * grid_key_d is the key to access any element in the grid
 *
 * \param dim dimensionality of the grid
 * \param p object property to get from the element of the grid
 *
 */

template<unsigned int dim, unsigned int p>
class grid_key_d
{
public:
  
  
  template<typename a, typename ...T>grid_key_d(a v,T...t)
  {
    k[dim-1] = v;
    invert_assign(t...);
  }
  
  template<typename a, typename ...T>void invert_assign(a v,T...t)
  {
    k[sizeof...(T)] = v;
    invert_assign(t...);
  }

  template<typename a, typename ...T>void invert_assign(a v)
  {
    k[0] = v;
  }
  
  mem_id k[dim];
};


/*! \brief grid_key is the key to access any element in the grid
 *
 * grid_key is the key to access any element in the grid
 *
 * \param p dimensionality of the grid
 *
 */

template<unsigned int p>
class grid_key
{
public:
  
  template<unsigned int s>grid_key(grid_key<s> & key)
  {
    id = key.id;
  }
  
  grid_key(size_t sz)
  :id(new mem_id[sz])
  {
  }
  
  void set_d(size_t id_, mem_id val)
  {
    id[id_] = val;
  }
  
  void set(mem_id val1)
  {
    id[0] = val1;
  }
  
  void set(mem_id val1, mem_id val2)
  {
    id[0] = val2;
    id[1] = val1;
  }
  
  void set(mem_id val1, mem_id val2, mem_id val3)
  {
    id[0] = val3;
    id[1] = val2;
    id[2] = val1;
  }
  
  void set(mem_id val1, mem_id val2, mem_id val3, mem_id val4)
  {
    id[0] = val4;
    id[1] = val3;
    id[2] = val2;
    id[3] = val1;
  }
  
  mem_id * getId()
  {
    return id.get();
  }
  
  boost::shared_array<mem_id> id;
};

//#pragma openfpm create(layout)

/*! \brief class that store the information at runtime of the grid plus define the linearization
 *
 * class that store the information at runtime of the grid plus define the linearization
 *
 * \param N dimensionality
 * \param T type of object is going to store the grid
 *
 */

template<unsigned int N, typename T>
class grid
{
	//! total number of the elements in the grid
  size_t size_tot;

  //! size of the grid
  size_t sz[N];

  //! size of the grid on each stride (used for linearization)
  size_t sz_s[N];
  
public:
  
  /*! \brief construct a grid from another grid
   *
   * construct a grid from another grid, type can be different
   *
   */

  template<typename S> grid(grid<N,S> g)
  {
	  // copy all the members

	  size_tot = g.size_tot;

	  for (int i = 0 ; i < N ; i++)
	  {sz[i] = g.sz[i]; sz_s[i] = g.sz_s[i];}
  }

  // Static element to calculate total size

  size_t totalSize(std::vector<size_t> & sz)
  {
    size_t tSz = 1;
    
    for (size_t i = 0 ;  i < sz.size() ; i++)
    {
      tSz *= sz[i];
    }
    
    return tSz;
  }

  /*! \brief Construct a grid of a specified size
   *
   * Construct a grid of a specified size
   *
   * \param sz is an std::vector that contain the size of the grid on each dimension
   *
   */

  grid(std::vector<size_t> & sz)
  : size_tot(totalSize(sz))
  {
    sz_s[0] = sz[0];
    this->sz[0] = sz[0];
    for (size_t i = 1 ;  i < sz.size() ; i++)
    {
      sz_s[i] = sz[i]*sz_s[i-1];
      this->sz[i] = sz[i];
    }
  }
  
  /*! \brief Linearization of the grid_key_dx
   *
   * Linearization of the grid_key_dx given a key, it spit out a number that is just the 1D linearization
   * of the key. In this case is the linearization of N index
   *
   * \param grid_key_dx<dim> grid key to access the element on a key
   *
   */

  #pragma openfpm layout(get)
  template<unsigned int dim> mem_id LinId(grid_key_dx<dim> & gk)
  {
    mem_id lid = gk.k[0];
    for (mem_id i = 1 ; i < dim ; i++)
    {
      lid += gk.k[i] * sz_s[i-1];
    }
    
    return lid;
  }
  
  /*! \brief Linearization of the grid_key_d
   *
   * Linearization of the grid_key_d given a key, it spit out a number that is just the 1D linearization
   * of the key. In this case is the linearization of N index
   *
   * \param grid_key_d<dim,p> grid key to access the element on a key
   *
   */

  #pragma openfpm layout(get)
  template<unsigned int dim, unsigned int p> mem_id LinId(grid_key_d<dim,p> & gk)
  {
    mem_id lid = gk.k[0];
    for (mem_id i = 1 ; i < dim ; i++)
    {
      lid += gk.k[i] * sz_s[i-1];
    }
    
    return lid;
  }
  
  /*! \brief Linearization of an array of mem_id (long int)
   *
   * Linearization of an array of mem_id, it spit out a number that is just the 1D linearization
   * of the key. In this case is the linearization of N index
   *
   * \param id an array of mem_id index
   *
   */

  #pragma openfpm layout(get)
  mem_id LinId(mem_id * id)
  {
    mem_id lid = 0;
    lid += id[0];
    for (mem_id i = 1 ; i < N ; i++)
    {
      lid += id[i] * sz_s[i-1];
    }
    
    return lid;
  }
  
  //! Destructor
  ~grid() {};
  
  /*! \brief Return the size of the grid
   *
   * Return the size of the grid
   *
   * \return the size of the grid
   *
   */

  #pragma openfpm layout(size)
  size_t size()
  {
    return size_tot;
  };

  /**
   *
   * Get the size of the grid on the direction i
   *
   * \param i direction
   * \return the size on the direction i
   *
   */

  size_t size(unsigned int i)
  {
	  return sz[i];
  }

  //!  It simply mean that all the classes grid are friend of all its specialization
  template <unsigned int,typename> friend class grid;
};

/**
 *
 * Grid key class iterator, iterate through the grid element
 *
 * \param dim dimensionality of the grid
 *
 * Usage: In general you never create object directly, but you get it from a grid_cpu or grid_gpu with
 *        getIterator()
 *
 */

template<unsigned int dim>
class grid_key_dx_iterator
{
	grid<dim,void> grid_base;

	grid_key_dx<dim> gk;

public:

	/*! \brief Constructor require a grid
	 *
	 * Constructor require a grid<dim,T>
	 *
	 * \param T type of object that the grid store
	 *
	 */
	template<typename T> grid_key_dx_iterator(grid<dim,T> & g)
	: grid_base(g)
	{
		//! Initialize to 0 the index

		for (int i = 0 ; i < dim ; i++)
		{gk.set_d(i,0);}
	}

	/*! \brief Get the next element
	 *
	 * Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	grid_key_dx_iterator<dim> operator++()
	{
		//! increment the first index

		size_t id = gk.get(0);
		gk.set_d(0,id+1);

		//! check the overflow of all the index with exception of the last dimensionality

		int i = 0;
		for ( ; i < dim-1 ; i++)
		{
			size_t id = gk.get(i);
			if (id >= grid_base.size(i))
			{
				// ! overflow, increment the next index

				gk.set_d(i,0);
				id = gk.get(i+1);
				gk.set_d(i+1,id+1);
			}
			else
			{
				break;
			}
		}

		return *this;
	}

	/*! \brief Check if there is the next element
	 *
	 * Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */

	bool isEnd()
	{
		if (gk.get(dim-1) < grid_base.size(dim-1))
		{
			//! we did not reach the end of the grid

			return true;
		}

		//! we reach the end of the grid
		return false;
	}

	/*! \brief Get the actual key
	 *
	 * Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	grid_key_dx<dim> get()
	{
		return gk;
	}
};



