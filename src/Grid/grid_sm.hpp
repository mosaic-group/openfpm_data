#ifndef GRID_HPP
#define GRID_HPP

#include "config.h"
#include <boost/shared_array.hpp>
#include <vector>
#include <initializer_list>
#include <array>
#include "memory/memory.hpp"
#include "Space/Shape/Box.hpp"
#include "grid_key.hpp"
#include <iostream>

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

// Declarations;

template<unsigned int N, typename T> class grid_sm;
template <unsigned int dim> class print_warning_on_adjustment;
template<unsigned int dim,typename warn=print_warning_on_adjustment<dim>> class grid_key_dx_iterator_sub;

/*! \brief class that store the information of the grid like number of point on each direction and
 *  define the index linearization by stride
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
		for(size_t i = 0 ; i < N ; i++)
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

		for (size_t i = 1 ;  i < N ; i++)
		{
			sz_s[i] = sz[i]*sz_s[i-1];
			this->sz[i] = sz[i];

			// set the box
			box.setHigh(i,sz[i]);
			box.setLow(i,0);
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

		for (size_t i = 1 ;  i < N ; i++)
		{
			sz_s[i] = sz[i]*sz_s[i-1];
			this->sz[i] = sz[i];

			// set the box
			box.setHigh(i,sz[i]);
			box.setLow(i,0);
		}
	}

public:

	/*! \brief Return the box enclosing the grid
	 *
	 */
	const Box<N,size_t> & getBox() const
	{
		return box;
	}

	/*! \brief Reset the dimension of the grid
	 *
	 * \param dims std::vector that store on each dimension the size of the grid
	 *
	 */

	void setDimensions(std::vector<size_t> & dims)
	{
		Initialize(dims);
		size_tot = totalSize(dims);
	}

	/*! \brief Reset the dimension of the grid
	 *
	 * \param dims store on each dimension the size of the grid
	 *
	 */
	void setDimensions(const size_t  (& dims)[N])
	{
		Initialize(dims);
		size_tot = totalSize(dims);
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

		for (size_t i = 0 ; i < N ; i++)
		{sz[i] = g.sz[i]; sz_s[i] = g.sz_s[i];}
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
		for(size_t i = 0 ; i < N ; i++)
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

	grid_sm(const size_t (& sz)[N])
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
		for (size_t i = 1 ;  i < sz.size() && N > 1 ; i++)
		{
			sz_s[i] = sz[i]*sz_s[i-1];
			this->sz[i] = sz[i];
		}
	}

	/*! \brief Linearization of the grid_key_dx with a specified shift
	 *
	 * \tparam check class that check the linearization, if this check fail the function return -1
	 * \param gk grid_key_dx to linearize
	 * \param sum_id shift on each dimension
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

	mem_id LinIdPtr(size_t * k) const
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
	 * \param k grid key to access the element on the grid
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
	 * \param gk grid key to access the element of the grid
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
	template<typename a, typename ...lT>mem_id Lin(a v,lT...t) const
	{
#ifdef DEBUG
		if (sizeof...(t)+1 > N)
		{
			std::cerr << "Error incorrect grid cannot linearize more index than its dimensionality" << "\n";
		}
#endif

		return v*sz_s[sizeof...(t)-1] + Lin(t...);
	}

	//! Linearize a set of index
	template<typename a>mem_id Lin(a v) const
	{
		return v;
	}

	//! Construct

	/*! \brief inversion of the linearization of the grid_key_dx
	 *
	 * \param id of the object
	 * \return key of the grid that id identify
	 *
	 */
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
	 * \param gk grid key to access the element on a key
	 * \return index of the memory
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
	size_t size() const
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

		for (size_t i = 0 ; i < N ; i++)
		{
			sz[i] = g.sz[i];
			sz_s[i] = g.sz_s[i];
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
	const size_t (& getSize() const)[N]
	{
		return sz;
	}

	/*! \brief Return a sub-grid iterator
	 *
	 * Return a sub-grid iterator, to iterate through the grid
	 *
	 * \param start start point
	 * \param stop stop point
	 *
	 */
	inline grid_key_dx_iterator_sub<N> getSubIterator(grid_key_dx<N> & start, grid_key_dx<N> & stop) const
	{
		return grid_key_dx_iterator_sub<N>(*this,start,stop);
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
 * \note if you have a grid you can get this object from getIterator()
 *
 * ### Grid iterator declaration and usage
 * \snippet grid_unit_tests.hpp Grid iterator test usage
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
	 * \warning entremly unsafe
	 * Before use the iterator you have call reinitialize
	 *
	 */
	grid_key_dx_iterator()
	{
#ifdef DEBUG
		initialized = false;
#endif
	}

	/*! \brief Constructor from a grid_key_dx_iterator<dim>
	 *
	 * \param g_it grid_key_dx_iterator<dim>
	 */
	grid_key_dx_iterator(const grid_key_dx_iterator<dim> & g_it)
	: grid_base(g_it.grid_base)
	{
		//! Initialize to 0 the index

		for (size_t i = 0 ; i < dim ; i++)
		{
			gk.set_d(i,g_it.get_gk(i));
		}

#ifdef DEBUG
		initialized = true;
#endif
	}

	/*! \brief Constructor require a grid_sm<dim,T>
	 *
	 * \param g info of the grid on which iterate
	 */
	template<typename T> grid_key_dx_iterator(const grid_sm<dim,T> & g)
	: grid_base(g)
	{
		reset();

#ifdef DEBUG
		initialized = true;
#endif
	}

	/*! \brief Constructor from another grid_key_dx_iterator
	 *
	 * \param key_it grid_key_dx_iterator
	 */
	grid_key_dx_iterator<dim> operator=(const grid_key_dx_iterator<dim> & key_it)
	{
		grid_base = key_it.grid_base;

		//! Initialize the index using key_it

		for (size_t i = 0 ; i < dim ; i++)
		{gk.set_d(i,key_it.get_gk(i));}

		return *this;
	}

	/*! \brief Get the next element
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

		size_t i = 0;
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
	 * \param d is the dimension
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
		if (gk.get(dim-1) < (long int)grid_base.size(dim-1))
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
		grid_base = key.grid_base;
		reset();
	}

	/*! \brief Reset the iterator (it restart from the beginning)
	 *
	 */
	void reset()
	{
		//! Initialize to 0 the index

		for (size_t i = 0 ; i < dim ; i++)
		{gk.set_d(i,0);}
	}
};


/**
 *
 * Grid key class iterator, iterate through a starting linearized grid element
 * to a stop linearized grid element in particular if linearize is the function
 *  that linearize all the grid_key_dx, it create an iterator that pass through
 *  Linearize^(-1)(start) Linearize^(-1)(start+1) ....... Linearize^(-1)(stop)
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
	//! stop point
	grid_key_dx<dim> gk_stop;

public:

	/*! \brief Constructor require a grid grid<dim,T>
	 *
	 * It construct an iterator from one index to another, in particular
	 * if linearize is the function that linearize all the grid_key, it
	 * create an iterator that pass through Linearize^(-1)(start)
	 * Linearize^(-1)(start+1) ....... Linearize^(-1)(stop)
	 *
	 * For example for start (1,1) and stop (3,3) the point indicated with # are
	 * explored by the iterator
	 *
	 * \verbatim
	 *
                      +-----+-----+-----+-----+-----+-----+ (6,5)
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      #-----#-----#-----#-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      #-----#-----#-----#-----#-----#-----#
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----#-----#-----#-----#-----#-----#
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                    (0,0)
	 *
	 *
	 * \endverbatim
	 *
	 *
	 * \tparam T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 * \param from starting point
	 * \param to end point
	 *
	 */
	template<typename T> grid_key_dx_iterator_sp(grid_sm<dim,T> & g, mem_id from, mem_id to)
	:grid_key_dx_iterator<dim>(g)
	 {
		//! Convert to a grid_key
		this->gk = g.InvLinId(from);

		//! Convert to a grid_key
		gk_stop = g.InvLinId(to);
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
		for (int i = dim-1 ; i >= 0 ; i-- )
		{
			//! check if we still have points
			if (this->gk.get(i) < gk_stop.get(i))
				return true;
			else if (this->gk.get(i) > gk_stop.get(i))
				return false;
		}

		//! (Final point) we we still have one point
		return true;
	}
};

/* \brief grid_key_dx_iterator_sub can adjust the domain if the border go out-of-side
 *        in this case a warning is produced
 *
 * \tparam dim dimensionality
 *
 */
template <unsigned int dim>
class print_warning_on_adjustment
{
public:
	inline static void pw1(size_t i, const grid_key_dx<dim> & gk_start) {std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " start with index smaller than zero x[" << i << "]=" << gk_start.get(i) << "\n";}
	inline static void pw2(size_t i, const grid_key_dx<dim> & gk_start) {std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " Cropping start point x[" << i << "]=" << gk_start.get(i) << " to x[" << i << "]=0 \n"  ;}
	inline static void pw3() {std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " stop with smaller index than start \n";}
	inline static void pw4(size_t i, const grid_key_dx<dim> & gk_stop, const grid_sm<dim,void> & grid_base) {std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " stop index bigger than cell domain x[" << i << "]=" << gk_stop.get(i) << " >= " << grid_base.size(i) << "\n";}
	inline static void pw5() {std::cerr << "Warning grid_key_dx_iterator_sub: " << __FILE__ << ":" << __LINE__ << " the starting point is not smaller or equal than than the stop point in all the coordinates" << "\n";}
};

/* \brief grid_key_dx_iterator_sub can adjust the domain if the border go out-of-side
 *        in this case a warning is produced
 *
 * \tparam dim dimensionality
 * \tparam gb type of grid
 *
 */
template <unsigned int dim>
class do_not_print_warning_on_adjustment
{
public:
	inline static void pw1(size_t i, const grid_key_dx<dim> & gk_start) {}
	inline static void pw2(size_t i, const grid_key_dx<dim> & gk_start) {}
	inline static void pw3() {}
	inline static void pw4(size_t i, const grid_key_dx<dim> & gk_stop, const grid_sm<dim,void> & grid_base) {}
	inline static void pw5() {}
};

/**
 *
 * Grid key class iterator, iterate through a sub-grid defined by an hyper-cube
 *
 * \param dim dimensionality of the grid
 *
 * \note if you have a grid you can get this object from getIteratorSub()
 *
 * ### Sub grid iterator declaration and usage
 * \snippet grid_unit_tests.hpp Sub-grid iterator test usage
 *
 */

template<unsigned int dim, typename warn>
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

		for (size_t i = 0 ; i < dim ; i++)
		{
			// if start smaller than 0
			if (gk_start.get(i) < 0)
			{
#ifdef DEBUG
				warn::pw1(i,gk_start);
#endif
				if (gk_start.get(i) < gk_stop.get(i))
				{
#ifdef DEBUG
					warn::pw2(i,gk_start);
#endif
					gk_start.set_d(i,0);
				}
				else
				{
#ifdef DEBUG
					warn::pw3();
#endif
					// No points are available
					gk_start.set_d(dim-1,gk_stop.get(dim-1)+1);
					break;
				}
			}

			// if stop bigger than the domain
			if (gk_stop.get(i) >= (long int)grid_base.size(i))
			{
#ifdef DEBUG
				warn::pw4(i,gk_stop,grid_base);
#endif

				if (gk_start.get(i) < (long int)grid_base.size(i))
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
	 * \warning extremely unsafe
	 * If you use this constructor before use the iterator you should call reinitialize first
	 *
	 */
	grid_key_dx_iterator_sub()
{}

	/*! \brief Constructor from another grid_key_dx_iterator_sub
	 *
	 * It construct an iterator over an hyper-cube defined by start and stop,
	 * \warning if start and stop are outside the domain defined by g the intersection
	 * will be considered
	 *
	 * \param g_s_it grid_key_dx_iterator_sub
	 *
	 */
	grid_key_dx_iterator_sub(const grid_key_dx_iterator_sub<dim> & g_s_it)
	: grid_key_dx_iterator_sub(g_s_it.grid_base,g_s_it.gk_start,g_s_it.gk_stop)
	{}


	/*! \brief Constructor require a grid grid<dim,T>
	 *
	 * It construct an iterator over an hyper-cube defined by start and stop,
	 * \warning if start and stop are outside the domain defined by g the intersection
	 * will be considered
	 *
	 * \tparam T type of object that the grid store
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
				warn::pw5();
			}
		}
#endif

		Initialize();
	}


	/*! \brief Constructor require a grid grid<dim,T>
	 *
	 * It construct an iterator over an hyper-cube defined by start and stop,
	 * \warning if start and stop are outside the domain defined by g the intersection
	 * will be considered
	 *
	 * \tparam T type of object that the grid store
	 *
	 * \param g info of the grid where we are iterating
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

	/*! \brief Constructor require a grid grid<dim,T>
	 *
	 * It construct an iterator over an hyper-cube defined by start and stop,
	 *
	 * \tparam T type of object that the grid store
	 *
	 * \param g Grid on which iterate
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	template<typename T> grid_key_dx_iterator_sub(const grid_sm<dim,T> & g, const size_t (& start)[dim], const size_t (& stop)[dim])
			: grid_key_dx_iterator<dim>(g),grid_base(g),gk_start(start), gk_stop(stop)
			  {
#ifdef DEBUG
		//! If we are on debug check that the stop grid_key id bigger than the start
		//! grid_key

		for (unsigned int i = 0 ; i < dim ; i++)
		{
			if (start[i] > stop[i])
			{
				warn::pw5();
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

		size_t i = 0;
		for ( ; i < dim-1 ; i++)
		{
			size_t id = this->gk.get(i);
			if ((long int)id > gk_stop.get(i))
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

	inline bool isNext()
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
	inline grid_key_dx<dim> get()
	{
#ifdef DEBUG
		if (initialized == false)
		{std::cerr << "Error: " << __FILE__ << __LINE__ << " using unitialized iterator" << "\n";}
#endif

		return grid_key_dx_iterator<dim>::get();
	}

	/*! \brief Reinitialize the iterator
	 *
	 * it re-initialize the iterator with the passed grid_key_dx_iterator_sub
	 * the actual position of the grid_key_dx_iterator_sub is ignored
	 *
	 * \param g_s_it grid_key_dx_iterator_sub
	 *
	 */

	inline void reinitialize(const grid_key_dx_iterator_sub<dim> & g_s_it)
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
				warn::pw5();
			}
		}

		initialized = true;
#endif

		Initialize();
	}

	/*! \brief Get the volume spanned by this sub-grid iterator
	 *
	 * \return the volume
	 *
	 */
	inline size_t getVolume()
	{
		return Box<dim,long int>::getVolumeKey(gk_start.k, gk_stop.k);
	}

	/*! \brief Reset the iterator (it restart from the beginning)
	 *
	 */
	inline void reset()
	{
		//! Initialize to 0 the index

		for (size_t i = 0 ; i < dim ; i++)
		{this->gk.set_d(i,gk_start.get(i));}
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

	/*! \brief constructor from another key
	 *
	 * \param key
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

	/*! \brief Return the actual key
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
