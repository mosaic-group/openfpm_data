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
#include "util/mathutil.hpp"
#include "iterators/stencil_type.hpp"

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
	/*! \brief Check if vertex exist
	 *
	 * Check if exist
	 *
	 * \param v_id Vertex id
	 * \param sz Size limit for the vertex id
	 *
	 * \return true if exist
	 *
	 */
	static bool valid(size_t v_id, size_t sz)
	{
		return v_id < sz;
	}
};

template<unsigned int dim>
struct ite_gpu
{
#ifdef CUDA_GPU

	dim3 thr;
	dim3 wthr;

	grid_key_dx<dim,int> start;
	grid_key_dx<dim,int> stop;

	size_t nblocks()
	{
		return wthr.x * wthr.y * wthr.z;
	}

#endif
};

//! Declaration grid_sm
template<unsigned int N, typename T> class grid_sm;

template<unsigned int dim, typename T2, typename T>
ite_gpu<dim> getGPUIterator_impl(const grid_sm<dim,T2> & g1, grid_key_dx<dim,T> & key1, grid_key_dx<dim,T> & key2, size_t n_thr = 1024);

//! Declaration print_warning_on_adjustment
template <unsigned int dim> class print_warning_on_adjustment;

//! Declaration grid_key_dx_iterator_sub
template<unsigned int dim,typename stencil=no_stencil,typename warn=print_warning_on_adjustment<dim>> class grid_key_dx_iterator_sub;

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

	/*! \brief Initialize the basic structure
	 *
	 * Initialize the basic structure
	 *
	 * \param sz vector that store the size of the grid on each
	 *           dimensions
	 *
	 */

	inline void Initialize(const size_t sz)
	{
		//! Initialize the basic structure for each dimension
		sz_s[0] = sz;
		this->sz[0] = sz;

		// set the box
		box.setHigh(0,sz);
		box.setLow(0,0);

		for (size_t i = 1 ;  i < N ; i++)
		{
			/* coverity[dead_error_begin] */
			sz_s[i] = sz*sz_s[i-1];
			this->sz[i] = sz;

			// set the box
			box.setHigh(i,sz);
			box.setLow(i,0);
		}
	}

	/*! \brief Initialize the basic structure
	 *
	 * Initialize the basic structure
	 *
	 * \param sz vector that store the size of the grid on each
	 *           dimensions
	 *
	 */

	inline void Initialize(const size_t (& sz)[N])
	{
		//! Initialize the basic structure for each dimension
		sz_s[0] = sz[0];
		this->sz[0] = sz[0];

		// set the box
		box.setHigh(0,sz[0]);
		box.setLow(0,0);

		for (size_t i = 1 ;  i < N ; i++)
		{
			/* coverity[dead_error_begin] */
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

	inline void Initialize()
	{
		//! Initialize the basic structure for each dimension
		sz_s[0] = 0;
		this->sz[0] = 0;

		// set the box
		box.setHigh(0,0);
		box.setLow(0,0);

		for (size_t i = 1 ;  i < N ; i++)
		{
			/* coverity[dead_error_begin] */
			sz_s[i] = sz[i]*sz_s[i-1];

			// set the box
			box.setHigh(i,sz[i]);
			box.setLow(i,0);
		}
	}

public:

	/*! \brief Return the box enclosing the grid
	 *
	 * \return the box
	 *
	 */
	inline Box<N,size_t> getBox() const
	{
		return box;
	}

	/*! \brief Return the box enclosing the grid
	 *
	 * While getBox return as P2 the size of the grid
	 * getBoxKey return the size - 1 equivalent to the maximum valid point
	 * that does not overflow the grid
	 *
	 * \return the box
	 *
	 */
	inline const Box<N,size_t> getBoxKey() const
	{
		Box<N,size_t> bx;

		for (size_t i = 0 ; i < N ; i++)
		{
			bx.setLow(i,box.getLow(i));
			bx.setHigh(i,box.getHigh(i) - 1);
		}

		return bx;
	}

	/*! \brief Reset the dimension of the grid
	 *
	 * \param dims store on each dimension the size of the grid
	 *
	 */
	inline void setDimensions(const size_t  (& dims)[N])
	{
		Initialize(dims);
		size_tot = totalSize(dims);
	}

	/*! \brief Default constructor
	 *
	 * It produce a grid of size 0 on each dimension
	 *
	 */

	inline grid_sm()
	:size_tot(0)
	{
		// Initialize sz
		for (size_t i = 0 ; i < N ; i++)
		{sz[i] = 0;}

		Initialize();
	}

	/*! \brief construct a grid from another grid
	 *
	 * \param g grid info
	 *
	 * construct a grid from another grid, type can be different
	 *
	 */

	template<typename S> inline grid_sm(const grid_sm<N,S> & g)
	{
		size_tot = g.size_tot;

		for (size_t i = 0 ; i < N ; i++)
		{
			sz[i] = g.sz[i];
			sz_s[i] = g.sz_s[i];
		}
	}

	// Static element to calculate total size

	inline size_t totalSize(const size_t sz)
	{
		size_t tSz = 1;

		for (size_t i = 0 ;  i < N ; i++)
		{
			tSz *= sz;
		}

		return tSz;
	}

	// Static element to calculate total size

	inline size_t totalSize(const size_t (& sz)[N])
	{
		size_t tSz = 1;

		for (size_t i = 0 ;  i < N ; i++)
		{
			tSz *= sz[i];
		}

		return tSz;
	}


	/*! \brief Construct a grid of a specified size
	 *
	 * Construct a grid of a specified size
	 *
	 * \param sz is an array that contain the size of the grid on each dimension
	 *
	 */

	inline grid_sm(const size_t & sz)
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

	inline grid_sm(const size_t (& sz)[N])
	: size_tot(totalSize(sz))
	{
		Initialize(sz);
	}

	/*! \brief Linearization of the grid_key_dx with a specified shift
	 *
	 * \tparam check class that check the linearization, if this check fail the function return -1
	 * \param gk grid_key_dx to linearize
	 * \param sum_id shift on each dimension
	 *
	 * \return The linearization of the gk key shifted by c, or -1 if the check fail
	 */

	template<typename check=NoCheck, typename ids_type>
	inline mem_id LinId(const grid_key_dx<N,ids_type> & gk, const char sum_id[N]) const
	{
		mem_id lid;

		// Check the sum produce a valid key

		if (check::valid(gk.get(0) + sum_id[0],sz[0]) == false)
			return -1;

		lid = gk.get(0) + sum_id[0];


		for (mem_id i = 1 ; i < N ; i++)
		{
			// Check the sum produce a valid key

			if (check::valid(gk.get(i) + sum_id[i],sz[i]) == false)
				return -1;

			lid += (gk.get(i) + sum_id[i]) * sz_s[i-1];
		}

		return lid;
	}

	/*! \brief Linearization of the grid_key_dx with a specified shift
	 *
	 * \tparam check class that check the linearization, if this check fail the function return -1
	 * \param gk grid_key_dx to linearize
	 * \param sum_id shift on each dimension
	 * \param bc boundary conditions
	 *
	 * \return The linearization of the gk key shifted by c, or -1 if the check fail
	 */

	template<typename check=NoCheck,typename ids_type>
	inline mem_id LinId(const grid_key_dx<N,ids_type> & gk, const char sum_id[N], const size_t (&bc)[N]) const
	{
		mem_id lid;

		// Check the sum produce a valid key

		if (bc[0] == NON_PERIODIC)
		{
			if (check::valid(gk.get(0) + sum_id[0],sz[0]) == false)
				return -1;

			lid = gk.get(0) + sum_id[0];
		}
		else
		{
			lid = openfpm::math::positive_modulo(gk.get(0) + sum_id[0],sz[0]);
		}

		for (mem_id i = 1 ; i < N ; i++)
		{
			// Check the sum produce a valid key

			/* coverity[dead_error_line] */
			if (bc[i] == NON_PERIODIC)
			{
				if (check::valid(gk.get(i) + sum_id[i],sz[i]) == false)
					return -1;

				lid += (gk.get(i) + sum_id[i]) * sz_s[i-1];
			}
			else
			{
				lid += (openfpm::math::positive_modulo(gk.get(i) + sum_id[i],sz[i])) * sz_s[i-1];
			}
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
	inline mem_id LinIdPtr(size_t * k) const
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

	__device__ __host__ inline mem_id LinId(const size_t (& k)[N]) const
	{
		mem_id lid = k[0];
		for (mem_id i = 1 ; i < N ; i++)
		{
			/* coverity[dead_error_line] */
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

	template<typename ids_type> __device__ __host__ inline mem_id LinId(const grid_key_dx<N,ids_type> & gk) const
	{
		mem_id lid = gk.get(0);
		for (mem_id i = 1 ; i < N ; i++)
		{
			/* coverity[dead_error_line] */
			lid += gk.get(i) * sz_s[i-1];
		}

		return lid;
	}

	/*! \brief linearize an arbitrary set of index
	 *
	 * linearize an arbitrary set of index
	 *
	 */
	template<typename a, typename ...lT> inline mem_id Lin(a v,lT...t) const
	{
#ifdef SE_CLASS1
		if (sizeof...(t)+1 > N)
		{
			std::cerr << "Error incorrect grid cannot linearize more index than its dimensionality" << "\n";
		}
#endif

		return v*sz_s[sizeof...(t)-1] + Lin(t...);
	}

	//! Linearize a set of index
	template<typename a> __device__ __host__ inline mem_id Lin(a v) const
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
	__device__ __host__ inline grid_key_dx<N> InvLinId(mem_id id) const
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
	template<unsigned int dim, unsigned int p> inline mem_id LinId(const grid_key_d<dim,p> & gk) const
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
	inline mem_id LinId(mem_id * id) const
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
	__device__ __host__ inline size_t size() const
	{
		return size_tot;
	};

	/*! \brief Copy the grid from another grid
	 *
	 * \param g grid from witch to copy
	 *
	 */

	__device__ __host__ inline grid_sm<N,T> & operator=(const grid_sm<N,T> & g)
	{
		size_tot = g.size_tot;

		for (size_t i = 0 ; i < N ; i++)
		{
			sz[i] = g.sz[i];
			sz_s[i] = g.sz_s[i];
		}

		box = g.box;

		return *this;
	}

	/*! \brief Check if the two grid_sm are the same
	 *
	 * \param g element to check
	 *
	 * \return true if they are the same
	 *
	 */

	inline bool operator==(const grid_sm<N,T> & g)
	{
		for (size_t i = 0 ; i < N ; i++)
		{
			if (sz[i] != g.sz[i])
				return false;
		}

#ifdef SE_CLASS1

		if (size_tot != g.size_tot)
			return false;

		for (size_t i = 0 ; i < N ; i++)
		{
			if (sz_s[i] != g.sz_s[i])
				return false;
		}

#endif
		return true;
	}

	/*! \brief Check if the two grid_sm are the same
	 *
	 * \param g element to check
	 *
	 */

	inline bool operator!=(const grid_sm<N,T> & g)
	{
		return ! this->operator==(g);
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

	inline size_t size_s(unsigned int i) const
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

	__device__ __host__ inline size_t size(unsigned int i) const
	{
		return sz[i];
	}

	/*! \brief Return the size of the grid as an array
	 *
	 * \return get the size of the grid as an array
	 *
	 */
	inline const size_t (& getSize() const)[N]
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

#ifdef CUDA_GPU

	/*! \brief Get an iterator for the GPU
	 *
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	template<typename T2>
	struct ite_gpu<N> getGPUIterator(grid_key_dx<N,T2> & key1, grid_key_dx<N,T2> & key2, size_t n_thr = 1024) const
	{
		return getGPUIterator_impl<N>(*this,key1,key2,n_thr);
	}

	/*! \brief Get an iterator for the GPU
	 *
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	struct ite_gpu<N> getGPUIterator(size_t n_thr = 1024) const
	{
		grid_key_dx<N> k1;
		grid_key_dx<N> k2;

		for (size_t i = 0 ; i < N ; i++)
		{
			k1.set_d(i,0);
			k2.set_d(i,size(i));
		}

		return getGPUIterator_impl<N>(*this,k1,k2,n_thr);
	}

#endif

	/*! \brief swap the grid_sm informations
	 *
	 * \param g grid to swap
	 *
	 */
	inline void swap(grid_sm<N,T> & g)
	{
		size_t tmp = size_tot;
		size_tot = g.size_tot;
		g.size_tot = tmp;

		for (size_t i = 0 ; i < N ; i++)
		{
			tmp = sz[i];
			sz[i] = g.sz[i];
			g.sz[i] = tmp;

			tmp = sz_s[i];
			sz_s[i] = g.sz_s[i];
			g.sz_s[i] = tmp;
		}
	}

	/*! \brief Produce a string from the object
	 *
	 * \return string
	 *
	 */
	std::string toString() const
	{
		std::stringstream str;

		for (size_t i = 0 ; i < N ; i++)
			str << "sz[" << i << "]=" << size(i) << " ";

		return str.str();
	}

	//!  It simply mean that all the classes grid are friend of all its specialization
	template <unsigned int,typename> friend class grid_sm;
};


template<unsigned int dim, typename T2, typename T>
ite_gpu<dim> getGPUIterator_impl(const grid_sm<dim,T2> & g1, grid_key_dx<dim,T> & key1, grid_key_dx<dim,T> & key2, size_t n_thr)
{
	size_t tot_work = 1;
	for (size_t i = 0 ; i < dim ; i++)
	{tot_work *= key2.get(i) - key1.get(i) + 1;}

	size_t n = (tot_work <= n_thr)?openfpm::math::round_big_2(tot_work):n_thr;

	// Work to do
	ite_gpu<dim> ig;

	if (tot_work == 0)
	{
		ig.thr.x = 0;
		ig.thr.y = 0;
		ig.thr.z = 0;

		ig.wthr.x = 0;
		ig.wthr.y = 0;
		ig.wthr.z = 0;

		return ig;
	}

	ig.thr.x = 1;
	ig.thr.y = 1;
	ig.thr.z = 1;

	int dir = 0;

	while (n != 1)
	{
		if (dir % 3 == 0)
		{ig.thr.x = ig.thr.x << 1;}
		else if (dir % 3 == 1)
		{ig.thr.y = ig.thr.y << 1;}
		else if (dir % 3 == 2)
		{ig.thr.z = ig.thr.z << 1;}

		n = n >> 1;

		dir++;
		dir %= dim;
	}

	if (dim >= 1)
	{ig.wthr.x = (key2.get(0) - key1.get(0) + 1) / ig.thr.x + (((key2.get(0) - key1.get(0) + 1)%ig.thr.x != 0)?1:0);}

	if (dim >= 2)
	{ig.wthr.y = (key2.get(1) - key1.get(1) + 1) / ig.thr.y + (((key2.get(1) - key1.get(1) + 1)%ig.thr.y != 0)?1:0);}
	else
	{ig.wthr.y = 1;}

	if (dim >= 3)
	{
		// Roll the other dimensions on z
		ig.wthr.z = 1;
		for (size_t i = 2 ; i < dim ; i++)
		{ig.wthr.z *= (key2.get(i) - key1.get(i) + 1) / ig.thr.z + (((key2.get(i) - key1.get(i) + 1)%ig.thr.z != 0)?1:0);}
	}
	else
	{ig.wthr.z = 1;}

	// crop if wthr == 1

	if (dim >= 1 && ig.wthr.x == 1)
	{ig.thr.x = (key2.get(0) - key1.get(0) + 1);}

	if (dim >= 2 && ig.wthr.y == 1)
	{ig.wthr.y = key2.get(1) - key1.get(1) + 1;}

	if (dim == 3 && ig.wthr.z == 1)
	{ig.wthr.z = key2.get(2) - key1.get(2) + 1;}

	for (size_t i = 0 ; i < dim ; i++)
	{
		ig.start.set_d(i,key1.get(i));
		ig.stop.set_d(i,key2.get(i));
	}

	return ig;
}


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
