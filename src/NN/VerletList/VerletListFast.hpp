/*
 * VerletListFast.hpp
 *
 *  Created on: Aug 16, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTFAST_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTFAST_HPP_

#include "VerletNNIterator.hpp"

#define VERLET_STARTING_NSLOT 128

#define VL_NON_SYMMETRIC 0
#define VL_SYMMETRIC 1
#define VL_SYMMETRIC_RED 2

#define WITH_RADIUS 3

/*! \brief Get the heighborhood iterator based on type
 *
 * \tparam dim dimensionality
 * \tparam T type of space
 * \tparam CellListImpl Cell-list implementation
 * \tparam type of neighborhood
 *
 */
template<unsigned int dim, typename T, typename CellListImpl, int type>
struct NNType
{
	/*! \brief Get the neighborhood
	 *
	 * \param xp Position of the particle p
	 * \param p id of the particle p
	 * \param cl Cell-list type implementation
	 * \param r_cut Cutoff radius
	 *
	 */
	static inline auto get(Point<dim,T> & xp, size_t p, CellListImpl & cl, T r_cut, openfpm::vector<Point<dim,T>> & v_pos, size_t g_m, const Box<dim,T> & dom) -> decltype(cl.template getNNIterator<NO_CHECK>(0))
	{
		return cl.template getNNIterator<NO_CHECK>(cl.getCell(xp));
	}
};


/*! \brief Get the heighborhood iterator based on type
 *
 * spacialization for the case with NN with radius
 *
 * \tparam dim dimensionality
 * \tparam T type of space
 * \tparam CellListImpl Cell-list implementation
 *
 */
template<unsigned int dim, typename T, typename CellListImpl>
struct NNType<dim,T,CellListImpl,WITH_RADIUS>
{
	/*! \brief Get the neighborhood
	 *
	 * \param xp Position of the particle p
	 * \param p id of the particle p
	 * \param cl Cell-list type implementation
	 * \param r_cut Cutoff radius
	 *
	 */
	static inline auto get(Point<dim,T> & xp, size_t p, CellListImpl & cl, T r_cut, openfpm::vector<Point<dim,T>> & v_pos, size_t g_m, const Box<dim,T> & dom) -> decltype(cl.template getNNIteratorRadius<NO_CHECK>(0,0.0))
	{
		return cl.template getNNIteratorRadius<NO_CHECK>(cl.getCell(xp),r_cut);
	}
};


/*! \brief Get the heighborhood iterator based on type
 *
 * spacialization for the case with NN symmetric
 *
 * \tparam dim dimensionality
 * \tparam T type of space
 * \tparam CellListImpl Cell-list implementation
 *
 */
template<unsigned int dim, typename T, typename CellListImpl>
struct NNType<dim,T,CellListImpl,VL_SYMMETRIC>
{
	/*! \brief Get the neighborhood
	 *
	 * \param xp Position of the particle p
	 * \param p id of the particle p
	 * \param cl Cell-list type implementation
	 * \param r_cut Cutoff radius
	 *
	 */
	static inline auto get(Point<dim,T> & xp, size_t p, CellListImpl & cl, T r_cut, openfpm::vector<Point<dim,T>> & v_pos, size_t g_m, const Box<dim,T> & dom) -> decltype(cl.template getNNIteratorSym<NO_CHECK>(0,0))
	{
		return cl.template getNNIteratorSym<NO_CHECK>(cl.getCell(xp),p);
	}
};


/*! \brief Get the heighborhood iterator based on type
 *
 * spacialization for the case with NN symmetric
 *
 * \tparam dim dimensionality
 * \tparam T type of space
 * \tparam CellListImpl Cell-list implementation
 *
 */
template<unsigned int dim, typename T, typename CellListImpl>
struct NNType<dim,T,CellListImpl,VL_SYMMETRIC_RED>
{
	/*! \brief Get the neighborhood
	 *
	 * \param xp Position of the particle p
	 * \param p id of the particle p
	 * \param cl Cell-list type implementation
	 * \param r_cut Cutoff radius
	 * \param v_pos vector of particle positions
	 * \param dom Processor domain
	 *
	 */
	static inline auto get(Point<dim,T> & xp, size_t p, CellListImpl & cl, T r_cut, openfpm::vector<Point<dim,T>> & v_pos, size_t g_m, const Box<dim,T> & dom) -> decltype(cl.template getNNIteratorSymRed<NO_CHECK>(0,0,0,openfpm::vector<Point<dim,T>>(),Box<dim,T>()))
	{
		return cl.template getNNIteratorSymRed<NO_CHECK>(cl.getCell(xp),p,g_m,v_pos,dom);
	}
};

/*! \brief Class for FAST Verlet list implementation
 *
 * This class implement the FAST Verlet list, fast but memory
 * expensive. The memory allocation is (M * N_nn_max)*sizeof(ele)
 *
 * * M = number of particles
 * * N_nn_max = maximum number of neighborhood
 * * ele = element the structure is storing
 *
 * Example of a 2D Cell list 6x6 structure with padding 1 without shift, cell indicated with p are padding cell
 * the origin of the cell or point (0,0) is marked with cell number 9
 *
 * \tparam dim Dimensionality of the space
 * \tparam T type of the space float, double ...
 * \tparam base Base structure that store the information
 *
 * ### Declaration of a Verlet list [VerS == VerletList<3,double,FAST>]
 * \snippet VerletList_test.hpp create verlet
 * ### Declaration of a Verlet list from external Cell-list [VerS == CellList<3,double,FAST>]
 * \snippet VerletList_test.hpp Fill external cell list
 * \snippet VerletList_test.hpp create verlet cell
 * ### Usage of Verlet-list
 * \snippet VerletList_test.hpp usage of verlet
 *
 */
template<unsigned int dim, typename T, typename transform, typename CellListImpl>
class VerletList<dim,T,FAST,transform,CellListImpl>
{
protected:


private:

	//! Number of slot for each particle. Or maximum number of particles for each particle
	size_t slot;

	//! number of neighborhood particles for each particle
	openfpm::vector<size_t> cl_n;

	//! Neighborhood indexes for each particle store (each particle can store a number
	//! of elements == slot)
	openfpm::vector<size_t> cl_base;

	//! Interlal cell-list
	CellListImpl cli;

	//! Realloc the vectors
	void realloc()
	{
		// we do not have enough slots reallocate the basic structure with more
		// slots
		openfpm::vector<size_t> cl_base_(2*slot * cl_n.size());

		// copy cl_base
		for (size_t i = 0 ; i < cl_n.size() ; i++)
		{
			for (size_t j = 0 ; j < cl_n.get(i) ; j++)
				cl_base_.get(2*i*slot + j) = cl_base.get(slot * i + j);
		}

		// Double the number of slots
		slot *= 2;

		// swap the memory
		cl_base.swap(cl_base_);
	}

	/*! \brief Fill the cell-list with data
	 *
	 * \param cli Cell-list
	 * \param pos vector of positions
	 *
	 */
	void initCl(CellListImpl & cli, openfpm::vector<Point<dim,T>> & pos)
	{
		cli.clear();

		for (size_t i = 0; i < pos.size() ; i++)
		{
			cli.add(pos.get(i), i);
		}
	}

	/*! \brief Create the Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cl Cell-list elements to use to construct the verlet list
	 * \param dom Processor domain
	 * \param opt options to create the verlet list like VL_SYMMETRIC or VL_NON_SYMMETRIC
	 *
	 */
	inline void create(openfpm::vector<Point<dim,T>> & pos, T r_cut, size_t g_m, CellListImpl & cl, const Box<dim,T> & dom, size_t opt)
	{
		if (opt == VL_SYMMETRIC)
			create_<decltype(cli.template getNNIteratorSym<NO_CHECK>(0,0)),VL_SYMMETRIC>(pos,r_cut,g_m,cli,dom,opt);
		else if (opt == VL_SYMMETRIC_RED)
			create_<decltype(cli.template getNNIteratorSymRed<NO_CHECK>(0,0,0,openfpm::vector<Point<dim,T>>(),Box<dim,T>())),VL_SYMMETRIC_RED>(pos,r_cut,g_m,cli,dom,opt);
		else
			create_<decltype(cli.template getNNIterator<NO_CHECK>(0)),VL_NON_SYMMETRIC>(pos,r_cut,g_m,cli,dom,opt);
	}

	/*! \brief Create the Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cl Cell-list elements to use to construct the verlet list
	 * \param dom Processor domain
	 *
	 */
	template<typename NN_type, int type> inline void create_(openfpm::vector<Point<dim,T>> & pos, T r_cut, size_t g_m, CellListImpl & cli, const Box<dim,T> & dom, size_t opt)
	{
		size_t end;

		// resize verlet to store the number of particles
		if (opt == VL_SYMMETRIC || opt == VL_SYMMETRIC_RED)
			end = pos.size();
		else
			end = g_m;

		cl_n.resize(end);
		cl_base.resize(end*slot);
		cl_n.fill(0);

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		for (size_t i = 0 ; i < end ; i++)
		{
			Point<dim,T> p = pos.template get<0>(i);

			// Get the neighborhood of the particle
			NN_type NN = NNType<dim,T,CellListImpl,type>::get(p,i,cli,r_cut,pos,g_m,dom);
			while (NN.isNext())
			{
				auto nnp = NN.get();

				Point<dim,T> q = pos.template get<0>(nnp);

				if (p.distance2(q) < r_cut2)
					addPart(i,nnp);

				// Next particle
				++NN;
			}
		}
	}

	/*! \brief Create the Verlet list from a given cell-list with a particular cut-off radius
	 *
	 * \param pos vector of positions of particles
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cl Cell-list elements to use to construct the verlet list
	 *
	 */
	inline void createR(openfpm::vector<Point<dim,T>> & pos, T r_cut, size_t g_m, CellListImpl & cl)
	{
		// resize verlet to store the number of particles
		cl_n.resize(g_m);
		cl_n.fill(0);
		cl_base.resize(g_m*slot);

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		for (size_t i = 0 ; i < g_m ; i++)
		{
			Point<dim,T> p = pos.template get<0>(i);

			// Get the neighborhood of the particle
			auto NN = cl.template getNNIteratorRadius<NO_CHECK>(cl.getCell(p),r_cut);
			while (NN.isNext())
			{
				auto nnp = NN.get();

				Point<dim,T> q = pos.template get<0>(nnp);

				if (p.distance2(q) < r_cut2)
					addPart(i,nnp);

				// Next particle
				++NN;
			}
		}
	}

	/*! \brief Add to the cell
	 *
	 * \param part_id Cell id where to add
	 * \param ele element to add
	 *
	 */
	inline void addPart(size_t part_id, size_t ele)
	{
		// Get the number of element the cell is storing

		size_t nl = getNNPart(part_id);

		if (nl + 1 >= slot)
			realloc();

		// we have enough slot to store another neighbor element

		cl_base.get(slot * part_id + cl_n.get(part_id)) = ele;
		cl_n.get(part_id)++;
	}

public:

	//! Object type that the structure store
	typedef size_t value_type;

	/*! \brief Return for how many particles has been constructed this verlet list
	 *
	 * \return number of particles
	 *
	 */
	size_t size()
	{
		return cl_n.size();
	}


	/*! Initialize the verlet list
	 *
	 * \param box Domain where this cell list is living
	 * \param dom Processor domain
	 * \param r_cut cut-off radius
	 * \param pos vector of particle positions
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param opt option to generate Verlet list
	 *
	 */
	void Initialize(const Box<dim,T> & box, const Box<dim,T> & dom, T r_cut, openfpm::vector<Point<dim,T>> & pos, size_t g_m, size_t opt = VL_NON_SYMMETRIC)
	{
		// Number of divisions
		size_t div[dim];

		Box<dim,T> bt = box;

		// Calculate the divisions for the Cell-lists
		cl_param_calculate(bt,div,r_cut,Ghost<dim,T>(0.0));

		// Initialize a cell-list
		cli.Initialize(bt,div);
		initCl(cli,pos);

		// create verlet
		create(pos,r_cut,g_m,cli,dom,opt);
	}

	/*! \brief update the Verlet list
	 *
	 * \param r_cut cutoff radius
	 * \param dom Processor domain
	 * \param pos vector of particle positions
	 * \param g_m ghost marker
	 * \param opt option to create the Verlet list
	 *
	 */
	void update(const Box<dim,T> & dom, T r_cut, openfpm::vector<Point<dim,T>> & pos, size_t & g_m, size_t opt)
	{
		initCl(cli,pos);
		create(pos,r_cut,g_m,cli,dom,opt);
	}

	/*! Initialize the verlet list from an already filled cell-list
	 *
	 * \param cli external Cell-list
	 * \param r_cut cutoff-radius
	 * \param pos vector of particle positions
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	void Initialize(const Box<dim,T> & dom, CellListImpl & cli, T r_cut, openfpm::vector<Point<dim,T>> & pos, size_t g_m)
	{
		cl_n.resize(g_m);
		cl_base.resize(g_m*slot);

		Point<dim,T> spacing = cli.getCellBox().getP2();

		// Create with radius or not
		bool wr = true;

		for (size_t i = 0 ; i < dim ; i++)
			wr &= r_cut <= spacing.get(i);

		if (wr == true)
			create_<decltype(cli.template getNNIterator<NO_CHECK>(0)),VL_NON_SYMMETRIC>(pos,r_cut,g_m,cli,dom);
		else
			create_<decltype(cli.template getNNIteratorRadius<NO_CHECK>(0,0.0)),WITH_RADIUS>(pos,r_cut,g_m,cli,dom);
	}

	//! Default Constructor
	VerletList()
	:slot(VERLET_STARTING_NSLOT)
	{};

	//! Copy constructor
	VerletList(const VerletList<dim,T,FAST,transform,CellListImpl> & cell)
	:slot(VERLET_STARTING_NSLOT)
	{
		this->operator=(cell);
	}

	//! Copy constructor
	VerletList(VerletList<dim,T,FAST,transform,CellListImpl> && cell)
	:slot(VERLET_STARTING_NSLOT)
	{
		this->operator=(cell);
	}


	/*! \brief Verlet-list constructor
	 *
	 * \param box Domain where this verlet-list is living
	 * \param r_cut cutoff radius
	 * \param mat Matrix transformation
	 * \param pad padding for the internal Cell-list padding
	 * \param slot maximum number of slots or maximum number of neighborhood per particle
	 *
	 */
	VerletList(Box<dim,T> & box, T r_cut, Matrix<dim,T> mat, const size_t pad = 1, size_t slot=STARTING_NSLOT)
	:slot(VERLET_STARTING_NSLOT),CellDecomposer_sm<dim,T,transform>(box,div,mat,box.getP1(),pad)
	{
		SpaceBox<dim,T> sbox(box);
		Initialize(sbox,r_cut,pad,slot);
	}

	/*! \brief Verlet-list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param r_cut cut-off radius
	 * \param pos vector position of particles
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param slot maximum number of slots (or maximum number each particle can have)
	 *
	 * \note the maximum number of particle per slot if just an indication for performance
	 *
	 */
	VerletList(Box<dim,T> & box, T r_cut, openfpm::vector<Point<dim,T>> & pos, size_t g_m, size_t slot=VERLET_STARTING_NSLOT)
	:slot(slot)
	{
		SpaceBox<dim,T> sbox(box);
		Initialize(sbox,r_cut,pos,g_m);
	}

	/*! \brief Cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param dom Simulation domain
	 * \param r_cut cut-off radius
	 * \param pos vector position of particles
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param slot maximum number of slots (or maximum number each particle can have)
	 *
	 * \note the maximum number of particle per slot if just an indication for performance
	 *
	 */
	VerletList(SpaceBox<dim,T> & box, Box<dim,T> & dom, T r_cut, openfpm::vector<Point<dim,T>> & pos, size_t g_m, size_t slot=VERLET_STARTING_NSLOT)
	:slot(slot)
	{
		Initialize(box,r_cut,pos);
	}

	/*! \brief Destructor
	 *
	 *
	 */
	~VerletList()
	{}

	/*! \brief Copy the verlet list
	 *
	 * \param vl verlet list to copy
	 *
	 * \return itself
	 *
	 */
	VerletList<dim,T,FAST,transform,CellListImpl> & operator=(VerletList<dim,T,FAST,transform,CellListImpl> && vl)
	{
		slot = vl.slot;

		cl_n.swap(vl.cl_n);
		cl_base.swap(vl.cl_base);

		return *this;
	}

	/*! \brief Copy a verlet list
	 *
	 * \param vl verlet-list to copy
	 *
	 */
	VerletList<dim,T,FAST,transform,CellListImpl> & operator=(const VerletList<dim,T,FAST,transform,CellListImpl> & vl)
	{
		slot = vl.slot;

		cl_n = vl.cl_n;
		cl_base = vl.cl_base;

		return *this;
	}

	/*! \brief Return the number of neighborhood particles for the particle id
	 *
	 * \param part_id id of the particle
	 *
	 * \return number of neighborhood particles for a particular particle id
	 *
	 */
	inline size_t getNNPart(size_t  part_id) const
	{
		return cl_n.get(part_id);
	}

	/*! \brief Get the neighborhood element j for the particle i
	 *
	 * \param i particle id
	 * \param j neighborhood j
	 *
	 * \return The element value
	 *
	 */
	inline size_t get(size_t i, size_t j) const
	{
		return cl_base.get(i * slot + j);
	}

	/*! \brief Swap the memory
	 *
	 * \param vl Verlet list with witch you swap the memory
	 *
	 */
	inline void swap(VerletList<dim,T,FAST,transform,CellListImpl> & vl)
	{
		cl_n.swap(vl.cl_n);
		cl_base.swap(vl.cl_base);

		size_t vl_slot_tmp = vl.slot;
		vl.slot = slot;
		slot = vl_slot_tmp;
	}

	/*! \brief Get the Neighborhood iterator
	 *
	 * It iterate across all the neighborhood particles of a selected particle
	 *
	 * \param part_id particle id
	 *
	 * \return an interator across the neighborhood particles
	 *
	 */
	template<unsigned int impl=NO_CHECK> inline VerletNNIterator<dim,VerletList<dim,T,FAST,transform,CellListImpl>> getNNIterator(size_t part_id)
	{
		VerletNNIterator<dim,VerletList<dim,T,FAST,transform,CellListImpl>> vln(part_id,*this);

		return vln;
	}

	/*! \brief Clear the cell list
	 *
	 */
	void clear()
	{
		slot = VERLET_STARTING_NSLOT;
		for (size_t i = 0 ; i < cl_n.size() ; i++)
			cl_n.get(i) = 0;
	}
};



#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTFAST_HPP_ */
