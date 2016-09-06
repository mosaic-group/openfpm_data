/*
 * grid_key_dx_iterator_sub.hpp
 *
 *  Created on: Dec 15, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SUB_HPP_
#define OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SUB_HPP_

#include "grid_key_dx_iterator.hpp"
#include "Grid/grid_key.hpp"

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
	inline static void pw4(size_t i, const grid_key_dx<dim> & gk_stop, const grid_sm<dim,void> & grid_base) {std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " stop index bigger than cell domain x[" << i << "]=" << gk_stop.get(i) << " > " << grid_base.size(i) << "\n";}
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
			gk_start.set_d(i,m);
			gk_stop.set_d(i,g.size(i)-m-1);
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

		long int i = 0;
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
	inline grid_key_dx<dim> get() const
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

	inline void reinitialize(const grid_key_dx_iterator_sub<dim,warn> & g_s_it)
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

	/*! \brief Return the grid information related to this grid
	 *
	 *
	 *
	 */
	inline const grid_sm<dim,void> & getGridInfo()
	{
		return grid_base;
	}
};


//// Specialization for dimension 0

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

template<typename warn>
class grid_key_dx_iterator_sub<0,warn>
{

public:

	grid_key_dx_iterator_sub()
	{}

	grid_key_dx_iterator<0> & operator++()
	{return *this;}

	inline bool isNext()
	{return false;}

	inline size_t getVolume()
	{return 0;}

	inline void reset()
	{}
};

#endif /* OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SUB_HPP_ */
