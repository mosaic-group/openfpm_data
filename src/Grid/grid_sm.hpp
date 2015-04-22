#ifndef GRID_HPP
#define GRID_HPP

#include "config.h"
#include <boost/shared_array.hpp>
#include <vector>
#include <initializer_list>
#include <array>
#include "memory.hpp"
#include "Space/Shape/Box.hpp"
#include "grid_key.hpp"

// Box need the definition of grid_key_dx_r


#define HARDWARE 1

/*! \brief Class to check if the edge can be created or not
 *
 * Class to check if the edge can be created or not, in this case no check is implemented
 *
 */

class NoCheck
{
public:
	/*! \brief No check is performed
	 *
	 * No check is performed
	 *
	 * \param v_id Vertex id
	 * \param sz Size limit for the vertex id
	 *
	 */
	static bool valid(size_t v_id, size_t sz)
	{
		return true;
	}
};

/*! \brief Class to check if the edge can be created or not
 *
 * Class to check if the edge can be created or not, in this case no check is implemented
 *
 */

class CheckExistence
{
public:
	/*! \brief Check is performed
	 *
	 * Check is performed
	 *
	 * \param v_id Vertex id
	 * \param sz Size limit for the vertex id
	 *
	 */
	static bool valid(size_t v_id, size_t sz)
	{
		return v_id < sz;
	}
};



//#pragma openfpm create(layout)

/*! \brief class that store the information at runtime of the grid and define the index linearization
 * by stride
 *
 * class that store the information at runtime of the grid plus define the linearization
 *
 * \param N dimensionality
 * \param T type of object is going to store the grid
 *
 */

template<unsigned int N, typename T>
class grid_sm
{
  //! Box enclosing the grid
  Box<N,size_t> box;

  //! total number of the elements in the grid
  size_t size_tot;

  //! size of the grid
  size_t sz[N];

  //! size of the grid on each stride (used for linearization)
  size_t sz_s[N];
  
  //! ghost margin, how far from the margin is the ghost layer bound (High bound)
  size_t mrgsH[N];

  //! ghost margin, how far from the margin is the ghost layer bound (Low bound)
  size_t mrgsL[N];

  /*! \brief It multiplicate two number and return the result
   *
   * It multiplicate two number and return the result, mainly used for LinId
   *
   * \param a operand 1
   * \param b operand 2
   *
   */

  inline size_t mulLin(size_t a, size_t b)
  {
	  return a*b;
  }

  /*! \brief Initialize the basic structure
   *
   * Initialize the basic structure
   *
   * \param sz vector that store the size of the grid on each
   *           dimensions
   *
   */

  void Initialize(std::vector<size_t> & sz)
  {
	  // Convert the vector to an array
	  size_t sz_a[N];

	  // Copy
	  for(int i = 0 ; i < N ; i++)
	  {
		  sz_a[i] = sz[i];
	  }

	  // Initialize
	  Initialize(sz_a);
  }

  /*! \brief Initialize the basic structure
   *
   * Initialize the basic structure
   *
   * \param sz vector that store the size of the grid on each
   *           dimensions
   *
   */

  void Initialize(const size_t (& sz)[N])
  {
	  //! Initialize the basic structure for each dimension
	  sz_s[0] = sz[0];
	  this->sz[0] = sz[0];

	  // set the box
      box.setHigh(0,sz[0]);
      box.setLow(0,0);
      mrgsH[0] = 0;
      mrgsL[0] = 0;

	  for (size_t i = 1 ;  i < N ; i++)
	  {
		  sz_s[i] = sz[i]*sz_s[i-1];
	      this->sz[i] = sz[i];

	      // set the box
	      box.setHigh(i,sz[i]);
	      box.setLow(i,0);

	      // High margin, Low margin
	      mrgsH[i] = 0;
	      mrgsL[i] = 0;
	  }
  }

  /*! \brief Initialize the basic structure
   *
   * Produce a grid of size 0 on each dimension
   *
   */

  void Initialize()
  {
	  //! Initialize the basic structure for each dimension
	  sz_s[0] = 0;
	  this->sz[0] = 0;

	  // set the box
      box.setHigh(0,0);
      box.setLow(0,0);
      mrgsH[0] = 0;
      mrgsL[0] = 0;

	  for (size_t i = 1 ;  i < N ; i++)
	  {
		  sz_s[i] = sz[i]*sz_s[i-1];
	      this->sz[i] = sz[i];

	      // set the box
	      box.setHigh(i,sz[i]);
	      box.setLow(i,0);

	      // High margin, Low margin
	      mrgsH[i] = 0;
	      mrgsL[i] = 0;
	  }
  }

public:
  
  /*! \brief Return the box enclosing the grid
   *
   */
  const Box<N,size_t> & getBox()
  {
	  return box;
  }

  /*! \brief Reset the dimension of the grid
   *
   * \param std::vector that store on each dimension the size of the grid
   *
   */

  void setDimensions(std::vector<size_t> & dims)
  {
	  Initialize(dims);
	  size_tot = totalSize(dims);
  }

  /*! \brief Reset the dimension of the grid
   *
   * \param std::vector that store on each dimension the size of the grid
   *
   */

  void setDimensions(const size_t  (& dims)[N])
  {
	  Initialize(dims);
	  size_tot = totalSize(dims);
  }

  /*! \brief Set the ghost layer margins High bound
   *
   * \param margin border
   *
   */

  void setGhostH(size_t margin[])
  {
	  for (size_t s = 0; s < N ; s++)
	  {
		  mrgsH[s] = margin[s];
	  }
  }

  /*! \brief Set the ghost layer margins Low bound
   *
   * \param margin border
   *
   */

  void setGhostL(size_t margin[])
  {
	  for (size_t s = 0; s < N ; s++)
	  {
		  mrgsL[s] = margin[s];
	  }
  }

  /*! \brief Return the point where the domain start
   *
   * Return the point where the domain start
   *
   */
  grid_key_dx<N> getDomainStart()
  {
	  //! Start key

	  grid_key_dx<N> key_start;

	  // Calculate the starting point of the domain

	  for (unsigned int i = 0 ; i < N ; i++)
	  {
		  key_start.set_d(i,mrgsL[i]);
	  }

	  return key_start;
  }

  /*! \brief Return the point where the domain stop
   *
   * Return the point where the domain stop
   *
   */

  grid_key_dx<N> getDomainStop()
  {
	  //! Stop key

	  grid_key_dx<N> key_stop;

	  for (unsigned int i = 0 ; i < N ; i++)
	  {
		  // Calculate the ending point
		  key_stop.set_d(i,sz[i]-mrgsH[i]);
	  }

	  return key_stop;
  }

  /*! \brief Return the point where the domain start and stop
   *
   * Return the point where the domain start and stop
   *
   * \param start point to set
   * \param stop point to set
   *
   */

  void getDomainStartStop(grid_key_dx<N>& start, grid_key_dx<N> & stop)
  {
	  // Iterate on all dimension and calculate the starting point and
	  // the ending point of the hyper-cube

	  for (unsigned int i = 0 ; i < N ; i++)
	  {
		  // Calculate the starting point
		  start.set_d(i,mrgsL[i]);

		  // Calculate the ending point
		  stop.set_d(i,sz[i]-mrgsH[i]);
	  }
  }

  /*! \brief Is linearize additive
   *
   * Is linearize a linear function, in this case for stride return true
   * because linearize respect the property
   *
   * Linearize(key1 + key2) = Linearize(key1) + Linearize(key2)
   *
   */

  bool isLinearizeLinear()
  {
	  return true;
  }

  /*! \brief Default constructor
   *
   * It produce a grid of size 0 on each dimension
   *
   */

  grid_sm()
  {
	  Initialize();
  }

  /*! \brief construct a grid from another grid
   *
   * construct a grid from another grid, type can be different
   *
   */

  template<typename S> grid_sm(const grid_sm<N,S> & g)
  {
	  // copy all the members

	  size_tot = g.size_tot;

	  for (int i = 0 ; i < N ; i++)
	  {sz[i] = g.sz[i]; sz_s[i] = g.sz_s[i]; mrgsL[i] = g.mrgsL[i] ; mrgsH[i] = g.mrgsH[i];}
  }

  // Static element to calculate total size

  size_t totalSize(const size_t (& sz)[N])
  {
    size_t tSz = 1;

    for (size_t i = 0 ;  i < N ; i++)
    {
      tSz *= sz[i];
    }

    return tSz;
  }

  // Static element to calculate total size

  size_t totalSize(const std::vector<size_t> & sz)
  {
	  // Convert the vector to an array
	  size_t sz_a[N];

	  // Copy
	  for(int i = 0 ; i < N ; i++)
	  {
		  sz_a[i] = sz[i];
	  }
    
	  return totalSize(sz_a);
  }

  /*! \brief Construct a grid of a specified size
   *
   * Construct a grid of a specified size
   *
   * \param sz is an std::vector that contain the size of the grid on each dimension
   *
   */

  grid_sm(std::vector<size_t> & sz)
  : size_tot(totalSize(sz))
  {
	  Initialize(sz);
  }
  
  /*! \brief Construct a grid of a specified size
   *
   * Construct a grid of a specified size
   *
   * \param sz is an array that contain the size of the grid on each dimension
   *
   */

  grid_sm(size_t (& sz)[N])
  : size_tot(totalSize(sz))
  {
	  Initialize(sz);
  }

  /*! \brief Construct a grid of a specified size
   *
   * Construct a grid of a specified size
   *
   * \param sz is an std::vector that contain the size of the grid on each dimension
   *
   */

  grid_sm(std::vector<size_t> && sz)
  : size_tot(totalSize(sz))
  {
    sz_s[0] = sz[0];
    this->sz[0] = sz[0];
    for (size_t i = 1 ;  i < sz.size() ; i++)
    {
      sz_s[i] = sz[i]*sz_s[i-1];
      this->sz[i] = sz[i];

      mrgsL[i] = 0;
      mrgsH[i] = 0;
    }
  }

  /*! \brief Linearization of the grid_key_dx
   *
   * Linearization of a shifted grid_key_dx
   *
   * \tparam class that check the linearization, if this check the function return -1
   * \param grid_key_dx to linearize
   * \param shift of the grid key
   *
   * \return The linearization of the gk key shifted by c, or -1 if the check fail
   */

  template<typename check=NoCheck> mem_id LinId(const grid_key_dx<N> & gk, char sum_id[N]) const
  {
	  // Check the sum produce a valid key

	  if (check::valid(gk.k[0] + sum_id[0],sz[0]) == false)
		  return -1;

	  mem_id lid = gk.k[0] + sum_id[0];
	  for (mem_id i = 1 ; i < N ; i++)
	  {
		  // Check the sum produce a valid key

		  if (check::valid(gk.k[i] + sum_id[i],sz[i]) == false)
			  return -1;

		  lid += (gk.k[i] + sum_id[i]) * sz_s[i-1];
	  }

	  return lid;
  }

  /*! \brief Linearization of the set of indexes
   *
   * Linearization of the set of indexes, it spit out a number that is just the 1D linearization.
   * In this case is the linearization of N index
   *
   * \param k set of indexes to linearize
   *
   */

  mem_id LinId(size_t * k) const
  {
    mem_id lid = k[0];
    for (mem_id i = 1 ; i < N ; i++)
    {
      lid += k[i] * sz_s[i-1];
    }

    return lid;
  }

  /*! \brief Linearization of the grid_key_dx
   *
   * Linearization of the grid_key_dx given a key, it spit out a number that is just the 1D linearization
   * of the key. In this case is the linearization of N index
   *
   * \param grid_key_dx<dim> grid key to access the element on a key
   *
   */

  mem_id LinId(const size_t (& k)[N]) const
  {
    mem_id lid = k[0];
    for (mem_id i = 1 ; i < N ; i++)
    {
      lid += k[i] * sz_s[i-1];
    }

    return lid;
  }

  /*! \brief Linearization of the grid_key_dx
   *
   * Linearization of the grid_key_dx given a key, it spit out a number that is just the 1D linearization
   * of the key. In this case is the linearization of N index
   *
   * \param grid_key_dx<dim> grid key to access the element on a key
   *
   */

  mem_id LinId(const grid_key_dx<N> & gk) const
  {
    mem_id lid = gk.k[0];
    for (mem_id i = 1 ; i < N ; i++)
    {
      lid += gk.k[i] * sz_s[i-1];
    }
    
    return lid;
  }
  
  /*! \brief linearize an arbitrary set of index
   *
   * linearize an arbitrary set of index
   *
   */
  template<typename a, typename ...lT>mem_id LinId(a v,lT...t) const
  {
#ifdef DEBUG
	  if (sizeof...(t)+1 > N)
	  {
		  std::cerr << "Error incorrect grid cannot linearize more index than its dimensionality" << "\n";
	  }
#endif

	  return v*sz_s[sizeof...(t)-1] + LinId(t...);
  }

  //! Linearize a set of index
  template<typename a>mem_id LinId(a v) const
  {
	  return v;
  }

  //! Construct

  /*! \brief inversion of the linearization of the grid_key_dx
   *
   * \param mem_id id of the object
   * \param grid_key, key of the grid that id identify
   *
   */

  //#pragma openfpm layout(get)
  grid_key_dx<N> InvLinId(mem_id id) const
  {
    // Inversion of linearize

	grid_key_dx<N> gk;

	for (mem_id i = 0 ; i < N ; i++)
	{
      gk.set_d(i,id % sz[i]);
      id /= sz[i];
	}

    return gk;
  }

  /*! \brief Linearization of the grid_key_d
   *
   * Linearization of the grid_key_d given a key, it spit out a number that is just the 1D linearization
   * of the key. In this case is the linearization of N index
   *
   * \param grid_key_d<dim,p> grid key to access the element on a key
   *
   */

  //#pragma openfpm layout(get)
  template<unsigned int dim, unsigned int p> mem_id LinId(const grid_key_d<dim,p> & gk) const
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

  //#pragma openfpm layout(get)
  mem_id LinId(mem_id * id) const
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
  ~grid_sm() {};
  
  /*! \brief Return the size of the grid
   *
   * Return the size of the grid
   *
   * \return the size of the grid
   *
   */

  //#pragma openfpm layout(size)
  size_t size()
  {
    return size_tot;
  };

  /*! \brief Copy the grid from another grid
   *
   * \param g grid from witch to copy
   *
   */

  inline grid_sm<N,T> & operator=(const grid_sm<N,T> & g)
  {
	  box = g.box;
	  size_tot = g.size_tot;

	  for (int i = 0 ; i < N ; i++)
	  {
		  sz[i] = g.sz[i];
		  sz_s[i] = g.sz_s[i];
		  mrgsH[i] = g.mrgsH[i];
		  mrgsL[i] = g.mrgsL[i];
	  }

	  return *this;
  }

  /**
   *
   * Get the stride-size of the grid on the direction i
   *
   * [Example] on a grid 16*16*16 it return 16,256
   *
   * \param i direction
   * \return the size on the direction i
   *
   */

  size_t size_s(unsigned int i) const
  {
	  return sz_s[i];
  }

  /**
   *
   * Get the size of the grid on the direction i
   *
   * \param i direction
   * \return the size on the direction i
   *
   */

  size_t size(unsigned int i) const
  {
	  return sz[i];
  }

  /*! \brief Return the size of the grid as an std::vector
   *
   * \return get the size of the grid as an std::vector
   *
   */
  std::vector<size_t> getVectorSize()
  {
	  std::vector<size_t> vect_sz;

	  for (int i = 0 ; i < N ; i++)
	  {
		  vect_sz.push_back(sz[i]);
	  }

	  return vect_sz;
  }

  /*! \brief Return the size of the grid as an array
   *
   * \return get the size of the grid as an array
   *
   */
  size_t (& getSize())[N]
  {
	  return sz;
  }

  //!  It simply mean that all the classes grid are friend of all its specialization
  template <unsigned int,typename> friend class grid_sm;
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
#ifdef DEBUG
	// Actual status of the iterator, when the iterator is not initialized cannot be used
	// and reinitialize must be called
	bool initialized = false;
#endif

	grid_sm<dim,void> grid_base;

	/*! \brief return the index i of the gk key
	 *
	 * \param i index to get
	 *
	 * \return index value
	 *
	 */

	size_t get_gk(size_t i) const
	{
		return gk.get(i);
	}

protected:

	grid_key_dx<dim> gk;

public:

	/*! \brief Default constructor
	 *
	 * \WARNING entremly unsafe
	 * Before use the iterator you have call reinitialize
	 *
	 */
	grid_key_dx_iterator()
	{
#ifdef DEBUG
		initialized = false;
#endif
	}

	/*! \brief Constructor require a grid
	 *
	 * Constructor require a grid<dim,T>
	 *
	 * \param T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 */
	grid_key_dx_iterator(const grid_key_dx_iterator<dim> & g_it)
	: grid_base(g_it.grid_base)
	{
		//! Initialize to 0 the index

		for (int i = 0 ; i < dim ; i++)
		{
			gk.set_d(i,g_it.get_gk(i));
		}

#ifdef DEBUG
		initialized = true;
#endif
	}

	/*! \brief Constructor require a grid
	 *
	 * Constructor require a grid<dim,T>
	 *
	 * \param T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 */
	template<typename T> grid_key_dx_iterator(const grid_sm<dim,T> & g)
	: grid_base(g)
	{
		//! Initialize to 0 the index

		for (int i = 0 ; i < dim ; i++)
		{gk.set_d(i,0);}

#ifdef DEBUG
		initialized = true;
#endif
	}

	/*! \brief Constructor require a grid
	 *
	 * Constructor require a grid<dim,T>
	 *
	 * \param T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 */
	grid_key_dx_iterator<dim> operator=(const grid_key_dx_iterator<dim> & key_it)
	{
		grid_base = key_it.grid_base;

		//! Initialize the index using key_it

		for (int i = 0 ; i < dim ; i++)
		{gk.set_d(i,key_it.get_gk(i));}

		return *this;
	}

	/*! \brief Get the next element
	 *
	 * Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	grid_key_dx_iterator<dim> & operator++()
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

	/*! \brief Set the dimension
	 *
	 * Set the dimension
	 *
	 * \param dim is the dimension
	 * \param sz set the counter to sz
	 *
	 */
	void set(int d, size_t sz)
	{
		// set the counter dim to sz

		gk.set_d(d,sz);
	}

	/*! \brief Check if there is the next element
	 *
	 * Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */

	bool isNext()
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
	const grid_key_dx<dim> & get()
	{
		return gk;
	}

	/*! \brief Reinitialize the grid_key_dx_iterator
	 *
	 * \param key form
	 *
	 */
	void reinitialize(const grid_key_dx_iterator<dim> & key)
	{

	}

	/*! \brief Reset the counter
	 *
	 */
	void reset()
	{
		//! Initialize to 0 the index

		for (int i = 0 ; i < dim ; i++)
		{gk.set_d(i,0);}
	}
};


/**
 *
 * Grid key class iterator, iterate through a starting grid element
 * to a stop grid element
 *
 * \param dim dimensionality of the grid
 *
 * Usage: In general you never create object directly, but you get it from a grid_cpu or grid_gpu with
 *        getIteratorLinStartStop()
 *
 */

template<unsigned int dim>
class grid_key_dx_iterator_sp : public grid_key_dx_iterator<dim>
{
	//! grid base where we are iterating
	grid_sm<dim,void> grid_base;

	//! stop point
	grid_key_dx<dim> gk_stop;

public:

	/*! \brief Constructor require a grid
	 *
	 * Constructor require a grid<dim,T>
	 *
	 * It construct an iterator from one index to another, in particular
	 * if linearize is the function that linearize all the grid_key, it
	 * create an iterator that pass through Linearize^(-1)(start)
	 * Linearize^(-1)(start+1) ....... Linearize^(-1)(stop)
	 *
	 * \param T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 * \param from starting point
	 * \param to end point
	 *
	 */
	template<typename T> grid_key_dx_iterator_sp(grid_sm<dim,T> & g, mem_id from, mem_id to)
	:grid_base(g)
	{
		//! Convert to a grid_key
		this->gk = g.InvLinId(from);

		//! Convert to a grid_key
		gk_stop = g.InfLinId(to);
	}

	/*! \brief Check if there is the next element
	 *
	 * Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */

	bool isNext()
	{
		//! for all dimensions
		for (int i = dim-1 ; i >= 0 ; i++ )
		{
			//! if the index overflow the stop point
			if (this->gk.get(i) > gk_stop.get(i))
			{
				return true;
			}
			else if (this->gk.get(i) < gk_stop.get(i))
			{
				// we have still point
				return false;
			}
		}

		//! we reach the end of the grid
		return false;
	}
};

/**
 *
 * Grid key class iterator, iterate through a subgrid defined by an hyper-cube
 *
 * \param dim dimensionality of the grid
 *
 * Usage: In general you never create object directly, but you get it from a grid_cpu or grid_gpu with
 *        getIteratorSub()
 *
 */

template<unsigned int dim>
class grid_key_dx_iterator_sub : public grid_key_dx_iterator<dim>
{
#ifdef DEBUG
	bool initialized = false;
#endif

	//! grid base where we are iterating
	grid_sm<dim,void> grid_base;

	//! start point
	grid_key_dx<dim> gk_start;

	//! stop point
	grid_key_dx<dim> gk_stop;

	/*! \brief Initialize gk
	 *
	 */
	void Initialize()
	{
		// Check that start and stop are inside the domain otherwise crop them

		for (int i = 0 ; i < dim ; i++)
		{
			// if start smaller than 0
			if (gk_start.get(i) < 0)
			{
				if (gk_start.get(i) < gk_stop.get(i))
					gk_start.set_d(i,0);
				else
				{
					// No points are available
					gk_start.set_d(dim-1,gk_stop.get(dim-1)+1);
					break;
				}
			}

			// if stop bigger than the domain
			if (gk_stop.get(i) >= grid_base.size(i))
			{
				if (gk_start.get(i) < grid_base.size(i))
					gk_stop.set_d(i,grid_base.size(i)-1);
				else
				{
					// No point are available
					gk_start.set_d(dim-1,gk_stop.get(dim-1)+1);
					break;
				}
			}
		}

		//! Initialize gk
		for (unsigned int i = 0 ; i < dim ; i++)
		{
			this->gk.set_d(i,gk_start.get(i));
		}

#ifdef DEBUG
		initialized = true;
#endif
	}

public:

	/*! \brief Default constructor
	 *
	 * WARNING: extremly unsafe
	 * If you use this constructor before use the iterator you should call reinitialize first
	 *
	 */
	grid_key_dx_iterator_sub()
	{}

	/*! \brief Constructor require a grid
	 *
	 * Constructor require a grid<dim,T>
	 *
	 * It construct an iterator over an hyper-cube defined by start and stop,
	 * \WARNING if start and stop are outside the domain defined by g the intersection
	 * will be considered
	 *
	 * \param T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	grid_key_dx_iterator_sub(const grid_key_dx_iterator_sub<dim> & g_s_it)
	: grid_key_dx_iterator_sub(g_s_it.grid_base,g_s_it.gk_start,g_s_it.gk_stop)
	{}


	/*! \brief Constructor require a grid
	 *
	 * Constructor require a grid<dim,T>
	 *
	 * It construct an iterator over an hyper-cube defined by start and stop,
	 * \WARNING if start and stop are outside the domain defined by g the intersection
	 * will be considered
	 *
	 * \param T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	template<typename T> grid_key_dx_iterator_sub(const grid_sm<dim,T> & g, const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop)
	: grid_key_dx_iterator<dim>(g),grid_base(g),gk_start(start), gk_stop(stop)
	{
#ifdef DEBUG
		//! If we are on debug check that the stop grid_key id bigger than the start
		//! grid_key

		for (unsigned int i = 0 ; i < dim ; i++)
		{
			if (gk_start.get(i) > gk_stop.get(i))
			{
				std::cerr << "Error grid_key_dx_iterator : the starting point of the grid cannot be bigger than the stop point at any coordinate" << "\n";
			}
		}
#endif

		Initialize();
	}


	/*! \brief Constructor require a grid
	 *
	 * Constructor require a grid<dim,T>
	 *
	 * It construct an iterator over an hyper-cube defined by start and stop,
	 * \WARNING if start and stop are outside the domain defined by g the intersection
	 * will be considered
	 *
	 * \param T type of object that the grid store
	 *
	 * \param g grid info we are iterating
	 * \param m Margin of the domain
	 *
	 */
	template<typename T> grid_key_dx_iterator_sub(const grid_sm<dim,T> & g, const size_t m)
	: grid_key_dx_iterator<dim>(g),grid_base(g)
	{
		// Initialize the start and stop point
		for (unsigned int i = 0 ; i < dim ; i++)
		{
			gk_start.set(i,m);
			gk_stop.set(i,g.size(i)-m);
		}

		//
		Initialize();
	}

	/*! \brief Constructor require a grid
	 *
	 * Constructor require a grid<dim,T>
	 *
	 * It construct an iterator over an hyper-cube defined by start and stop,
	 *
	 * \param T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	template<typename T> grid_key_dx_iterator_sub(const grid_sm<dim,T> & g, size_t (& start)[dim], size_t (& stop)[dim])
	: grid_key_dx_iterator<dim>(g),grid_base(g),gk_start(start), gk_stop(stop)
	{
#ifndef DEBUG
		//! If we are on debug check that the stop grid_key id bigger than the start
		//! grid_key

		for (unsigned int i = 0 ; i < dim ; i++)
		{
			if (start[i] > stop[i])
			{
				std::cerr << "Error grid_key_dx_iterator : the starting point of the grid cannot be bigger than the stop point at any coordinate" << "\n";
			}
		}
#endif

		Initialize();
	}

	/*! \brief Get the next element
	 *
	 * Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	grid_key_dx_iterator<dim> & operator++()
	{
#ifdef DEBUG
		if (initialized == false)
		{std::cerr << "Error: " << __FILE__ << __LINE__ << " using unitialized iterator" << "\n";}
#endif

		//! increment the first index

		size_t id = this->gk.get(0);
		this->gk.set_d(0,id+1);

		//! check the overflow of all the index with exception of the last dimensionality

		int i = 0;
		for ( ; i < dim-1 ; i++)
		{
			size_t id = this->gk.get(i);
			if (id > gk_stop.get(i))
			{
				// ! overflow, increment the next index

				this->gk.set_d(i,gk_start.get(i));
				id = this->gk.get(i+1);
				this->gk.set_d(i+1,id+1);
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

	bool isNext()
	{
#ifdef DEBUG
		if (initialized == false)
		{std::cerr << "Error: " << __FILE__ << __LINE__ << " using unitialized iterator" << "\n";}
#endif

		if (this->gk.get(dim-1) <= gk_stop.get(dim-1))
		{
			//! we did not reach the end of the grid

			return true;
		}

		//! we reach the end of the grid
		return false;
	}

	/*! \brief Return the actual grid key iterator
	 *
	 */
	grid_key_dx<dim> get()
	{
#ifdef DEBUG
		if (initialized == false)
		{std::cerr << "Error: " << __FILE__ << __LINE__ << " using unitialized iterator" << "\n";}
#endif

		return grid_key_dx_iterator<dim>::get();
	}

	/*! \brief Reinitialize the iterator
	 *
	 * it reinitialize the iterator with the passed grid_key_dx_iterator_sub, it became like a clone
	 *
	 * \param grid_key_dx_iterator_sub
	 *
	 */

	void reinitialize(const grid_key_dx_iterator_sub<dim> & g_s_it)
	{
		// Reinitialize the iterator

		grid_key_dx_iterator<dim>::reinitialize(g_s_it);
		grid_base = g_s_it.grid_base;
		gk_start = g_s_it.gk_start;
		gk_stop = g_s_it.gk_stop;


		#ifdef DEBUG
			//! If we are on debug check that the stop grid_key id bigger than the start
			//! grid_key

			for (unsigned int i = 0 ; i < dim ; i++)
			{
				if (gk_start.get(i) > gk_stop.get(i))
				{
					std::cerr << "Error grid_key_dx_iterator : the starting point of the grid cannot be bigger than the stop point at any coordinate" << "\n";
				}
			}

			initialized = true;
		#endif

		Initialize();
	}
};


/*! \brief Emulate grid_key_dx with runtime dimensionality
 *
 * Emulate grid_key_dx with runtime dimensionality
 *
 */

class grid_key_dx_r
{
	size_t dim;

public:

	/*! \brief Get the dimensionality of the key
	 *
	 * Get the dimensionality of the key
	 *
	 */

	size_t getDim()
	{
		return dim;
	}

	  /*! \brief constructor
	   *
	   * constructor
	   *
	   * \param dim Dimensionality
	   *
	   */
	  grid_key_dx_r(grid_key_dx_r & key)
	  :dim(key.dim)
	  {
		  // Allocate the key
		  k = new mem_id[dim];

		  // Copy the key
		  for(unsigned int i = 0 ; i < dim ; i++)
		  {
			  k[i] = key.k[i];
		  }
	  }

  /*! \brief constructor
   *
   * constructor
   *
   * \param dim Dimensionality
   *
   */
  grid_key_dx_r(size_t dim)
  :dim(dim)
  {
	  // Allocate the key
	  k = new mem_id[dim];
  }

  ~grid_key_dx_r()
  {
	  delete [] k;
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
  mem_id * k;

private:

  /*! \brief Recursively invert the assignment
   *
   * Recursively invert the assignment at compile-time
   *
   */
  template<typename a, typename ...T>void invert_assign(a v,T...t)
  {
    k[sizeof...(T)] = v;
    invert_assign(t...);
  }

  template<typename a, typename ...T>void invert_assign(a v)
  {
    k[0] = v;
  }

  void invert_assign()
  {
  }

};


/**
 *
 * Iterate through the elements (i1,i2,....,in) with i1 ... in unsigned integers
 * with the following constrain (i1>i2>......>in)
 *
 */

class Iterator_g_const
{
	size_t dim;

	//! size of the grid (the grid is assumed a square so equal on each dimension)
	size_t sz;

	// Actual grid_key position
	grid_key_dx_r gk;

public:

	/*! \brief Get the dimensionality of the iterator
	 *
	 * Get the dimensionality of the iterator
	 *
	 */

	size_t getDim()
	{
		return dim;
	}

	/*! \brief Constructor
	 *
	 * \param n Dimensionality (how many i1 ... in you have)
	 * \param sz Size of the grid on all dimensions range of the value i1 ... in can assume
	 *
	 */

	Iterator_g_const(size_t n, size_t sz)
	:dim(n),sz(sz),gk(n)
	{
		// fill gk with the first grid element that satisfied the constrain: 0,1,2,3... dim

		for (size_t i = 0 ; i < dim ; i++)
		{
			gk.set_d(i,dim-i-1);
		}
	}

	/*! \brief Get the next element
	 *
	 * Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	Iterator_g_const & operator++()
	{
		//! increment the first index

		gk.set_d(0,gk.get(0)+1);

		//! check the overflow of all the index with exception of the last dimensionality

		unsigned int i = 0;
		for ( ; i < dim-1 ; i++)
		{
			size_t id = gk.get(i);
			if (id >= sz)
			{
				// ! overflow, increment the next index

				id = gk.get(i+1);
				if (id+i+2 >= sz)
				{
					// there is no-way to produce a valid key
					// there is not need to check the previous index
					// overflow i+1

					gk.set_d(i+1,sz);
				}
				else
				{

					// reinitialize the previous index

					for (unsigned int s = 0 ; s <= i+1 ; s++)
					{
						gk.set_d(i+1-s,id+1+s);
					}
				}
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

	bool isNext()
	{
		// If dimensionless return immediately
		if (dim == 0)
			return false;

		if (gk.get(dim-1) < static_cast<mem_id>(sz-dim+1))
		{
			//! we did not reach the end of the grid

			return true;
		}

		//! we reach the end of the grid
		return false;
	}

	/*! brief Return the actual key
	 *
	 * Return the actual key
	 *
	 * \return The actual key that identify with the set of index
	 *
	 */

	grid_key_dx_r & get()
	{
		return gk;
	}

};

#endif
