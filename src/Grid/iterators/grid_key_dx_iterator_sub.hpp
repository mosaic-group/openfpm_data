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
template <unsigned int dim, typename linearizer>
class print_warning_on_adjustment
{
public:

	/*! \brief print warning type1
	 *
	 * \param gk_start start point
	 *
	 */
	inline static void pw1(size_t i, const grid_key_dx<dim> & gk_start) {std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " start with index smaller than zero x[" << i << "]=" << gk_start.get(i) << "\n";}

	/*! \brief print warning type2
	 *
	 * \param i component
	 * \param gk_start start point
	 *
	 */
	inline static void pw2(size_t i, const grid_key_dx<dim> & gk_start) {std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " Cropping start point x[" << i << "]=" << gk_start.get(i) << " to x[" << i << "]=0 \n"  ;}


	/*! \brief print warning type3
	 *
	 *
	 */
	inline static void pw3() {std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " stop with smaller index than start \n";}


	/*! \brief print warning type4
	 *
	 * \param i component
	 * \param gk_stop stop point
	 * \param grid_base grid information
	 *
	 */
	inline static void pw4(size_t i, const grid_key_dx<dim> & gk_stop, const linearizer & grid_base) {std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " stop index bigger than cell domain x[" << i << "]=" << gk_stop.get(i) << " > " << grid_base.size(i) << "\n";}

	/*! \brief print warning type5
	 *
	 *
	 */
	inline static void pw5() {std::cerr << "Warning grid_key_dx_iterator_sub: " << __FILE__ << ":" << __LINE__ << " the starting point is not smaller or equal than than the stop point in all the coordinates" << "\n";}
};

/* \brief grid_key_dx_iterator_sub can adjust the domain if the border go out-of-side
 *        in this case a warning is produced
 *
 * \tparam dim dimensionality
 * \tparam gb type of grid
 *
 */
template <unsigned int dim, typename linearizer>
class do_not_print_warning_on_adjustment
{
public:

	/*! \brief print warning type1
	 *
	 * \param i component
	 * \param gk_start start point
	 *
	 */
	inline static void pw1(size_t i, const grid_key_dx<dim> & gk_start) {}

	/*! \brief print warning type2
	 *
	 * \param i component
	 * \param gk_start start point
	 *
	 */
	inline static void pw2(size_t i, const grid_key_dx<dim> & gk_start) {}

	/*! \brief print warning type3
	 *
	 *
	 */
	inline static void pw3() {}

	/*! \brief print warning type4
	 *
	 * \param i component
	 * \param gk_stop stop point
	 * \param grid_base grid information
	 *
	 */
	inline static void pw4(size_t i, const grid_key_dx<dim> & gk_stop, const linearizer & grid_base) {}

	/*! \brief print warning type5
	 *
	 *
	 */
	inline static void pw5() {}
};

template<unsigned int dim, typename stl_type, typename linearizer>
struct post_increment_sub_impl
{
	static inline void inc(grid_key_dx<dim> & gk,
										grid_key_dx<dim> & gk_start,
										grid_key_dx<dim> & gk_stop,
										stl_type & stl_code,
										linearizer & grid_base)
	{
		long int i = 0;
		for ( ; i < dim-1 ; i++)
		{
			/* coverity[dead_error_begin] */
			size_t id = gk.get(i);
			if ((long int)id > gk_stop.get(i))
			{
				// ! overflow, increment the next index

				size_t idr = gk.get(i) - gk_start.get(i);
				gk.set_d(i,gk_start.get(i));
				id = gk.get(i+1);
				gk.set_d(i+1,id+1);

				stl_code.adjust_offset(i,idr,grid_base);
			}
			else
			{
				break;
			}
		}
	}
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

template<unsigned int dim, typename stencil, typename linearizer, typename warn>
class grid_key_dx_iterator_sub : public grid_key_dx_iterator<dim,stencil,linearizer>
{
#ifdef SE_CLASS1
	bool initialized = false;
#endif

	//! grid base where we are iterating
	linearizer grid_base;

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
#ifdef SE_CLASS1
				warn::pw1(i,gk_start);
#endif
				if (gk_start.get(i) < gk_stop.get(i))
				{
#ifdef SE_CLASS1
					warn::pw2(i,gk_start);
#endif
					gk_start.set_d(i,0);
				}
				else
				{
#ifdef SE_CLASS1
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
#ifdef SE_CLASS1
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

#ifdef SE_CLASS1
		initialized = true;
#endif
	}

	/*! \brief After incremented we have to check we did not overflo
	 *         any dimension and in case adjust the dimensions
	 *
	 */
	void post_increment()
	{
		//! check the overflow of all the index with exception of the last dimensionality

//		post_increment_sub_impl<dim,stencil>::inc(this->gk,gk_start,gk_stop,this->stl_code,grid_base);

		long int i = 0;
		for ( ; i < dim-1 ; i++)
		{
			/* coverity[dead_error_begin] */
			size_t id = this->gk.get(i);
			if ((long int)id > gk_stop.get(i))
			{
				// ! overflow, increment the next index

				size_t idr = this->gk.get(i) - gk_start.get(i);
				this->gk.set_d(i,gk_start.get(i));
				id = this->gk.get(i+1);
				this->gk.set_d(i+1,id+1);

				this->stl_code.adjust_offset(i,idr,grid_base);
			}
			else
			{
				break;
			}
		}
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
	 *
	 * \param g_s_it grid_key_dx_iterator_sub
	 *
	 */
	grid_key_dx_iterator_sub(const grid_key_dx_iterator_sub<dim,stencil,linearizer> & g_s_it)
	:grid_key_dx_iterator<dim,stencil>(g_s_it),grid_base(g_s_it.grid_base),gk_start(g_s_it.gk_start), gk_stop(g_s_it.gk_stop)
	{
#ifdef SE_CLASS1
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
	 * \param g Grid on which iterate
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	grid_key_dx_iterator_sub(const linearizer & g,
												  const grid_key_dx<dim> & start,
												  const grid_key_dx<dim> & stop)
	: grid_key_dx_iterator<dim,stencil,linearizer>(g),grid_base(g),gk_start(start), gk_stop(stop)
	{
#ifdef SE_CLASS1
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
	 * \param g Grid on which iterate
	 * \param start starting point
	 * \param stop end point
	 * \param stencil_pnt stencil points
	 *
	 */
	grid_key_dx_iterator_sub(const linearizer & g,
			                 const grid_key_dx<dim> & start,
							 const grid_key_dx<dim> & stop,
							 const grid_key_dx<dim> (& stencil_pnt)[stencil::nsp])
	:grid_key_dx_iterator<dim,stencil,linearizer>(g,stencil_pnt),grid_base(g),gk_start(start), gk_stop(stop)
	{
#ifdef SE_CLASS1
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
		grid_key_dx_iterator<dim,stencil>::calc_stencil_offset(this->gk);
	}

	/*! \brief Constructor
	 *
	 * It construct an iterator setting only the stencil point.
	 *  Require reinitialize to make it work
	 *
	 *
	 * \param stencil_pnt stencil points
	 *
	 */
	grid_key_dx_iterator_sub(const grid_key_dx<dim> (& stencil_pnt)[stencil::nsp])
	:grid_key_dx_iterator<dim,stencil>(stencil_pnt)
	{
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
	grid_key_dx_iterator_sub(const linearizer & g, const size_t m)
	:grid_key_dx_iterator<dim>(g),grid_base(g)
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
	 *
	 * \param g Grid on which iterate
	 * \param start starting point
	 * \param stop end point
	 *
	 */
	grid_key_dx_iterator_sub(const linearizer & g,
												  const size_t (& start)[dim],
												  const size_t (& stop)[dim])
	:grid_key_dx_iterator<dim>(g),grid_base(g),gk_start(start), gk_stop(stop)
	{
#ifdef SE_CLASS1
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
	grid_key_dx_iterator_sub<dim,stencil,linearizer,warn> & operator++()
	{
#ifdef SE_CLASS1
		if (initialized == false)
		{std::cerr << "Error: " << __FILE__ << __LINE__ << " using unitialized iterator" << "\n";}
#endif

		//! increment the first index

		size_t id = this->gk.get(0);
		this->gk.set_d(0,id+1);

		this->stl_code.increment();

		//! check the overflow of all the index with exception of the last dimensionality

		post_increment();

		return *this;
	}

	/*! \brief increment the operator by more than one
	 *
	 * \return the next grid_key
	 *
	 */
	grid_key_dx_iterator_sub<dim,stencil,warn> & operator+=(int nsteps)
	{
#ifdef SE_CLASS1
		if (initialized == false)
		{std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using unitialized iterator" << "\n";}
#endif

		for (size_t i = 0 ; i < nsteps ; i++)
		{
			this->operator ++();
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
#ifdef SE_CLASS1
		if (initialized == false)
		{std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using unitialized iterator" << "\n";}
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
	 * \return the actual key
	 *
	 */
	inline grid_key_dx<dim> get() const
	{
#ifdef SE_CLASS1
		if (initialized == false)
		{std::cerr << "Error: " << __FILE__ << __LINE__ << " using unitialized iterator" << "\n";}
#endif

		return grid_key_dx_iterator<dim,stencil,linearizer>::get();
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

		grid_key_dx_iterator<dim,stencil>::reinitialize(g_s_it);
		grid_base = g_s_it.getGridInfo();
		gk_start = g_s_it.getStart();
		gk_stop = g_s_it.getStop();


#ifdef SE_CLASS1
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
		grid_key_dx_iterator<dim,stencil>::calc_stencil_offset(this->gk);
	}

	/*! \brief Get the volume spanned by this sub-grid iterator
	 *
	 * \return the volume
	 *
	 */
	inline size_t getVolume()
	{
		return Box<dim,long int>::getVolumeKey(gk_start.get_k(), gk_stop.get_k());
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
	inline const linearizer & getGridInfo() const
	{
		return grid_base;
	}

	/*! \brief Starting point
	 *
	 * \return starting point
	 *
	 */
	const grid_key_dx<dim> & getStart() const
	{
		return gk_start;
	}

	/*! \brief Stop point
	 *
	 * \return stop point
	 *
	 */
	const grid_key_dx<dim> & getStop() const
	{
		return gk_stop;
	}

	/*! \brief Sum a template constant
	 *
	 * \tparam compile-time offset
	 *
	 */
	template<unsigned int tot_add>
	inline void private_sum()
	{
		this->stl_code.template private_sum<tot_add>();
	}

	/*! \brief Sum a template constant
	 *
	 * \param tot_add Add an offset to all the pointer
	 *
	 */
	inline void private_adjust(size_t tot_add)
	{
		this->stl_code.private_adjust(tot_add);
	}

	/* \brief Set the iterator in a way that isNext return false
	 *
	 */
	void invalidate()
	{
		this->gk.set_d(dim-1,1);
		this->gk_stop.set_d(dim-1,0);

#ifdef SE_CLASS1
		this->initialized = true;
#endif
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

	/*! \brief Next point
	 *
	 * stub version do nothing
	 *
	 * \return itself
	 *
	 */
	grid_key_dx_iterator<0> & operator++()
	{return *this;}

	/*! \brief Is there a next point
	 *
	 * stub version
	 *
	 * \return always false
	 *
	 */
	inline bool isNext()
	{return false;}

	/*! \brief return the volume
	 *
	 * stub version
	 *
	 * \return always 0
	 *
	 */
	inline size_t getVolume()
	{return 0;}

	/*! \brief Reset the iterator
	 *
	 * stub version, do nothing
	 *
	 *
	 */
	inline void reset()
	{}
};

#endif /* OPENFPM_DATA_SRC_GRID_ITERATORS_GRID_KEY_DX_ITERATOR_SUB_HPP_ */
