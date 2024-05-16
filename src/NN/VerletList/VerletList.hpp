/*
 * VerletList.hpp
 *
 *  Created on: Aug 16, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_

#include "Vector/map_vector.hpp"

#include "VerletNNIterator.hpp"
#include "NN/CellList/CellList.hpp"
#include "NN/CellList/CellList_util.hpp"
#include "NN/Mem_type/MemFast.hpp"
#include "NN/Mem_type/MemBalanced.hpp"
#include "NN/Mem_type/MemMemoryWise.hpp"

#define VERLETLIST_FAST(dim,St) VerletList<dim,St,Mem_fast<>,shift<dim,St> >
#define VERLETLIST_BAL(dim,St) VerletList<dim,St,Mem_bal<>,shift<dim,St> >
#define VERLETLIST_MEM(dim,St) VerletList<dim,St,Mem_mem<>,shift<dim,St> >

#define VERLET_STARTING_NSLOT 128

#ifdef LOCAL_INDEX64
typedef size_t local_index_;
#else
typedef unsigned int local_index_;
#endif


/*! \brief Class for Verlet list implementation
 *
 * * M = number of particles
 * * N_nn_max = maximum number of neighborhood
 * * ele = element the structure is storing
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
template<unsigned int dim,
		 typename T,
		 typename Mem_type = Mem_fast<HeapMemory,local_index_>,
		 typename transform = no_transform<dim,T>,
		 typename vector_pos_type = openfpm::vector<Point<dim,T>>,
		 typename CellListImpl = CellList<dim,T,Mem_fast<HeapMemory,typename Mem_type::local_index_type>,transform,vector_pos_type> >
class VerletList: public Mem_type
{
protected:

	//! Number of slot for each particle. Or maximum number of particles for each particle
	typename Mem_type::local_index_type slot;

	//! Domain particles
	openfpm::vector<typename Mem_type::local_index_type> domainParticlesCRS;

	//! Option flags
	size_t opt;

private:

	//! decomposition counter
	size_t n_dec;

	//! Interlal cell-list
	CellListImpl cli;


	/*! \brief Fill the cell-list with data
	 *
	 * \param cellList Cell-list
	 * \param vPos vector of positions
	 * \param ghostMarker marker
	 * \param opt VL_SYMMETRIC or VL_NON_SYMMETRIC
	 *
	 */
	void initCl(CellListImpl & cellList, vector_pos_type & vPos, size_t ghostMarker)
	{
		// CellList_gpu receives a property vector to potentially reorder it during cell list construction
		// stub because of legacy naming
		openfpm::vector<aggregate<int>> vPropStub;

		if (opt & VL_SYMMETRIC || opt & VL_CRS_SYMMETRIC)
			cellList.setOpt(CL_SYMMETRIC);
		else
			cellList.setOpt(CL_NON_SYMMETRIC);

		cellList.fill(vPos, vPropStub, ghostMarker);
	}

	/*! \brief Create the CRS Symmetric Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param pos2 vector of position for the neighborhood
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cli Cell-list elements to use to construct the verlet list
	 * \param dom list of domain cells with normal neighborhood
	 * \param anom list of domain cells with non-normal neighborhood
	 *
	 */
	inline void createCRSSymmetric(
		const vector_pos_type & pos,
		const vector_pos_type & pos2,
		const openfpm::vector<size_t> & dom,
		const openfpm::vector<subsub_lin<dim>>& anom,
		T r_cut,
		size_t ghostMarker,
		CellListImpl& cli)
	{
		size_t end = pos.size();

		Mem_type::init_to_zero(slot,end);
		domainParticlesCRS.clear();

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		auto it = ParticleItCRS_Cells<dim,CellListImpl,vector_pos_type>(cli,dom,anom,cli.getNNc_sym());
		while (it.isNext())
		{
			typename Mem_type::local_index_type i = it.get();
			Point<dim,T> xp = pos.template get<0>(i);

			domainParticlesCRS.add(i);

			// Get the neighborhood of the particle
			auto NN = it.getNNIteratorCSR(pos);
			while (NN.isNext())
			{
				auto nnp = NN.get();

				Point<dim,T> xq = pos2.template get<0>(nnp);

				if (xp.distance2(xq) < r_cut2)
					addPart(i,nnp);

				// Next particle
				++NN;
			}

			++it;
		}
	}

	/*! \brief Create the Symmetric Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param pos2 vector of position for the neighborhood
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cli Cell-list elements to use to construct the verlet list
	 *
	 */
	inline void createSymmetric(
		const vector_pos_type & pos,
		const vector_pos_type & pos2,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cli)
	{
		size_t end = ghostMarker;

		Mem_type::init_to_zero(slot,end);
		domainParticlesCRS.clear();

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		auto it = pos.getIteratorTo(end);
		while (it.isNext())
		{
			typename Mem_type::local_index_type i = it.get();
			Point<dim,T> xp = pos.template get<0>(i);

			// Get the neighborhood of the particle
			auto NN = cli.getNNIteratorSym(cli.getCell(xp),i,pos);
			while (NN.isNext())
			{
				auto nnp = NN.get();

				Point<dim,T> xq = pos2.template get<0>(nnp);

				if (xp.distance2(xq) < r_cut2)
					addPart(i,nnp);

				// Next particle
				++NN;
			}

			++it;
		}
	}

	/*! \brief Create the Non-symmetric Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param pos2 vector of position for the neighborhood
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cli Cell-list elements to use to construct the verlet list
	 *
	 */
	inline void createNonSymmetric(
		const vector_pos_type & pos,
		const vector_pos_type & pos2,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cli)
	{
		size_t end = ghostMarker;

		Mem_type::init_to_zero(slot,end);

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		auto it = pos.getIteratorTo(end);
		while (it.isNext())
		{
			typename Mem_type::local_index_type i = it.get();
			Point<dim,T> xp = pos.template get<0>(i);

			// Get the neighborhood of the particle
			auto NN = cli.getNNIterator(cli.getCell(xp));
			while (NN.isNext())
			{
				auto nnp = NN.get();

				Point<dim,T> xq = pos2.template get<0>(nnp);

				if (xp.distance2(xq) < r_cut2)
					addPart(i,nnp);

				// Next particle
				++NN;
			}

			++it;
		}
	}

	/*! \brief Create the Non-symmetric Verlet list from a given cell-list (Radius)
	 *
	 * \param pos vector of positions
	 * \param pos2 vector of position for the neighborhood
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cli Cell-list elements to use to construct the verlet list
	 *
	 */
	inline void createNonSymmetricRadius(
		const vector_pos_type & pos,
		const vector_pos_type & pos2,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cli)
	{
		size_t end = ghostMarker;

		Mem_type::init_to_zero(slot,end);

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		auto it = pos.getIteratorTo(end);
		while (it.isNext())
		{
			typename Mem_type::local_index_type i = it.get();
			Point<dim,T> xp = pos.template get<0>(i);

			// Get the neighborhood of the particle
			auto NN = cli.getNNIteratorRadius(cli.getCell(xp),r_cut);
			while (NN.isNext())
			{
				auto nnp = NN.get();

				Point<dim,T> xq = pos2.template get<0>(nnp);

				if (xp.distance2(xq) < r_cut2)
					addPart(i,nnp);

				// Next particle
				++NN;
			}

			++it;
		}
	}

public:

	//! type for the local index
	typedef Mem_type Mem_type_type;

	//! Object type that the structure store
	typedef size_t value_type;

	/*! \brief Return for how many particles has been constructed this verlet list
	 *
	 * \return number of particles
	 *
	 */
	size_t size()
	{
		return Mem_type::size();
	}

	/*! \brief Add a neighborhood particle to a particle
	 *
	 * \param part_id part id where to add
	 * \param ele element to add
	 *
	 */
	inline void addPart(size_t part_id, size_t ele)
	{
		Mem_type::addCell(part_id,ele);
	}

	/*! Initialize the verlet list for Non-Symmetric case
	 *
	 * \param box Domain where this cell list is living
	 * \param dom Processor domain
	 * \param r_cut cut-off radius
	 * \param pos vector of particle positions
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param opt option to generate Verlet list
	 *
	 */
	void Initialize(
		const Box<dim,T> & box,
		const Box<dim,T> & dom,
		T r_cut,
		vector_pos_type & pos,
		size_t ghostMarker)
	{
		// Number of divisions
		size_t div[dim];

		Box<dim,T> bt = box;

		// Calculate the divisions for the Cell-lists
		cl_param_calculate(bt,div,r_cut,Ghost<dim,T>(0.0));

		// Initialize a cell-list
		cli.Initialize(bt,div);

		initCl(cli,pos,ghostMarker);

		createNonSymmetric(pos,pos,r_cut,ghostMarker,cli);
	}

	/*! \brief Initialize the symmetric Verlet-list
	 *
	 * \param box Simulation domain
	 * \param dom Processor domain
	 * \param g ghost size
	 * \param r_cut cut-off radius
	 * \param pos vector of particle positions
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	void InitializeSym(const Box<dim,T> & box, const Box<dim,T> & dom, const Ghost<dim,T> & g, T r_cut, vector_pos_type & pos, size_t ghostMarker)
	{
		// Padding
		size_t pad = 0;

		// Cell decomposer
		CellDecomposer_sm<dim,T,shift<dim,T>> cd_sm;

		// Calculate the divisions for the Cell-lists
		cl_param_calculateSym<dim,T>(box,cd_sm,g,r_cut,pad);

		// Initialize a cell-list
		cli.Initialize(cd_sm,dom,pad);
		initCl(cli,pos,ghostMarker);

		// create verlet
		createSymmetric(pos, pos, r_cut, ghostMarker, cli);
	}


	/*! \brief Initialize the symmetric Verlet-list CRS scheme
	 *
	 * \param box Simulation domain
	 * \param dom Processor domain
	 * \param g ghost size
	 * \param r_cut cut-off radius
	 * \param pos vector of particle positions
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	void InitializeCrs(const Box<dim,T> & box, const Box<dim,T> & dom, const Ghost<dim,T> & g, T r_cut, vector_pos_type & pos, size_t ghostMarker)
	{
		// Padding
		size_t pad = 0;

		// Cell decomposer
		CellDecomposer_sm<dim,T,shift<dim,T>> cd_sm;

		// Calculate the divisions for the Cell-lists
		cl_param_calculateSym<dim,T>(box,cd_sm,g,r_cut,pad);

		// Initialize a cell-list
		cli.Initialize(cd_sm,dom,pad);
		initCl(cli,pos,ghostMarker);
	}

	/*! \brief Create the Verlet-list with the crossing scheme
	 *
	 * \param pos vector with the particle positions
	 * \param ghostMarker ghost marker
	 * \param pos vector with the particle positions
	 * \param r_cut cut-off radius
	 * \param dom_c domain cells
	 * \param anom_c cells with anomalos neighborhood
	 *
	 */
	void createVerletCrs(
		T r_cut,
		size_t ghostMarker,
		vector_pos_type & pos,
		openfpm::vector<size_t> & dom_c,
		openfpm::vector<subsub_lin<dim>> & anom_c)
	{
		createCRSSymmetric(pos, pos,dom_c,anom_c,r_cut,ghostMarker,cli);
	}

	/*! \brief update the Verlet list
	 *
	 * \param r_cut cutoff radius
	 * \param dom Processor domain
	 * \param pos vector of particle positions
	 * \param ghostMarker ghost marker
	 * \param opt option to create the Verlet list
	 *
	 */
	void update(const Box<dim,T> & dom, T r_cut, vector_pos_type & pos, size_t & ghostMarker)
	{
		initCl(cli,pos,ghostMarker);

		if (opt & VL_SYMMETRIC)
			createSymmetric(pos,pos,r_cut,ghostMarker,cli);
		else
			createNonSymmetric(pos,pos,r_cut,ghostMarker,cli);
	}

	/*! \brief update the Verlet list
	 *
	 * \param r_cut cutoff radius
	 * \param dom Processor domain
	 * \param pos vector of particle positions
	 * \param ghostMarker ghost marker
	 * \param dom_c list of cells with normal neighborhood
	 * \param anom_c list of cells with anormal neighborhood
	 *
	 */
	void updateCrs(
		const Box<dim,T> & dom,
		T r_cut,
		vector_pos_type & pos,
		size_t & ghostMarker,
		const openfpm::vector<size_t> & dom_c,
		const openfpm::vector<subsub_lin<dim>> & anom_c)
	{
		initCl(cli,pos,ghostMarker);
		createCRSSymmetric(pos,pos,dom_c,anom_c,r_cut,ghostMarker,cli);
	}

	/*! Initialize the verlet list from an already filled cell-list
	 *
	 * \param cli external Cell-list
	 * \param r_cut cutoff-radius
	 * \param pos vector of particle positions
	 * \param pos2 vector of particle position for the neighborhood
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * 	\param opt options for the Verlet-list creation
	 *
	 */
	void Initialize(CellListImpl & cli,
		T r_cut,
		const vector_pos_type & pos,
		const vector_pos_type & pos2,
		size_t ghostMarker)
	{
		Point<dim,T> spacing = cli.getCellBox().getP2();

		// Create with radius or not
		bool wr = true;

		for (size_t i = 0 ; i < dim ; i++)
			wr &= r_cut <= spacing.get(i);

		if (wr == true)
			createNonSymmetricRadius(pos,pos2,r_cut,ghostMarker,cli);
		else
			createNonSymmetric(pos,pos2,r_cut,ghostMarker,cli);
	}

	//! Default Constructor
	VerletList(size_t opt=VL_NON_SYMMETRIC)
	:Mem_type(VERLET_STARTING_NSLOT),slot(VERLET_STARTING_NSLOT),n_dec(0),opt(opt)
	{};

	//! Copy constructor
	VerletList(const VerletList<dim,T,Mem_type,transform,vector_pos_type,CellListImpl> & cell)
	:Mem_type(VERLET_STARTING_NSLOT),slot(VERLET_STARTING_NSLOT)
	{
		this->operator=(cell);
	}

	//! Copy constructor
	VerletList(VerletList<dim,T,Mem_type,transform,vector_pos_type,CellListImpl> && cell)
	:Mem_type(VERLET_STARTING_NSLOT),slot(VERLET_STARTING_NSLOT),n_dec(0)
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
	VerletList(Box<dim,T> & box, T r_cut, Matrix<dim,T> mat, const size_t pad = 1, size_t opt=VL_NON_SYMMETRIC, size_t slot=STARTING_NSLOT)
	:slot(VERLET_STARTING_NSLOT),CellDecomposer_sm<dim,T,transform>(box,div,mat,box.getP1(),pad),opt(opt)
	{
		Box<dim,T> sbox(box);
		Initialize(sbox,r_cut,pad,slot);
	}

	/*! \brief Verlet-list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param r_cut cut-off radius
	 * \param pos vector position of particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param slot maximum number of slots (or maximum number each particle can have)
	 *
	 * \note the maximum number of particle per slot if just an indication for performance
	 *
	 */
	VerletList(Box<dim,T> & box, T r_cut, vector_pos_type & pos, size_t ghostMarker, size_t opt=VL_NON_SYMMETRIC, size_t slot=VERLET_STARTING_NSLOT)
	:slot(slot),opt(opt)
	{
		Box<dim,T> sbox(box);
		Initialize(sbox,r_cut,pos,ghostMarker);
	}

	/*! \brief Cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param dom Simulation domain
	 * \param r_cut cut-off radius
	 * \param pos vector position of particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param slot maximum number of slots (or maximum number each particle can have)
	 *
	 * \note the maximum number of particle per slot if just an indication for performance
	 *
	 */
	VerletList(Box<dim,T> & box, Box<dim,T> & dom, T r_cut, vector_pos_type & pos, size_t ghostMarker, size_t opt=VL_NON_SYMMETRIC, size_t slot=VERLET_STARTING_NSLOT)
	:slot(slot),opt(VL_NON_SYMMETRIC),opt(opt)
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
	VerletList<dim,T,Mem_type,transform,vector_pos_type,CellListImpl> &
	operator=(VerletList<dim,T,Mem_type,transform,vector_pos_type,CellListImpl> && vl)
	{
		slot = vl.slot;

		Mem_type::operator=(vl);
		domainParticlesCRS.swap(vl.domainParticlesCRS);

		n_dec = vl.n_dec;
		opt = vl.opt;

		return *this;
	}

	/*! \brief Copy a verlet list
	 *
	 * \param vl verlet-list to copy
	 *
	 * \return itself
	 *
	 */
	VerletList<dim,T,Mem_type,transform,vector_pos_type,CellListImpl> &
	operator=(const VerletList<dim,T,Mem_type,transform,vector_pos_type,CellListImpl> & vl)
	{
		slot = vl.slot;
		opt = vl.opt;

		Mem_type::operator=(vl);

		cli = vl.cli;

		domainParticlesCRS = vl.domainParticlesCRS;
		n_dec = vl.n_dec;

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
		return Mem_type::getNelements(part_id);
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
		return Mem_type::get(i,j);
	}

	/*! \brief Swap the memory
	 *
	 * \param vl Verlet list with witch you swap the memory
	 *
	 */
	inline void swap(VerletList<dim,T,Mem_type,transform,vector_pos_type,CellListImpl> & vl)
	{
		Mem_type::swap(vl);
		domainParticlesCRS.swap(vl.domainParticlesCRS);

		size_t vl_slot_tmp = vl.slot;
		vl.slot = slot;
		slot = vl_slot_tmp;

		cli.swap(vl.cli);

		size_t n_dec_tmp = vl.n_dec;
		vl.n_dec = n_dec;
		n_dec = n_dec_tmp;

		size_t optTmp = vl.opt;
		vl.opt = opt;
		opt = optTmp;
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
	inline VerletNNIterator<dim,VerletList<dim,T,Mem_type,transform,vector_pos_type,CellListImpl>>
	getNNIterator(size_t part_id)
	{
		VerletNNIterator<dim,VerletList<dim,T,Mem_type,transform,vector_pos_type,CellListImpl>> vln(part_id,*this);

		return vln;
	}

	/*! \brief Clear the cell list
	 *
	 */
	void clear()
	{
		Mem_type::clear();
	}

	/*! \brief Return the starting point of the neighborhood for the particle p
	 *
	 * \param part_id particle id
	 *
	 * \return the index
	 *
	 */
	inline const typename Mem_type::local_index_type &
	getStart(typename Mem_type::local_index_type part_id)
	{
		return Mem_type::getStartId(part_id);
	}

	/*! \brief Return the end point of the neighborhood for the particle p
	 *
	 * \param part_id particle id
	 *
	 * \return the stop index
	 *
	 */
	inline const typename Mem_type::local_index_type &
	getStop(typename Mem_type::local_index_type part_id)
	{
		return Mem_type::getStopId(part_id);
	}

	/*! \brief Return the neighborhood id
	 *
	 * \param part_id particle id
	 *
	 * \return the neighborhood id
	 *
	 */
	inline const typename Mem_type::local_index_type &
	get_lin(const typename Mem_type::local_index_type * part_id)
	{
		return Mem_type::get_lin(part_id);
	}

	/*! \brief Get the internal cell-list used to construct the Verlet-list
	 *
	 * \return the internal Cell-list
	 *
	 */
	CellListImpl & getInternalCellList()
	{
		return cli;
	}

	/*! \brief Set the n_dec number
	 *
	 * \param n_dec
	 *
	 */
	void set_ndec(size_t n_dec)
	{
		this->n_dec = n_dec;

		cli.set_ndec(n_dec);
	}

	/*! \brief Set the n_dec number
	 *
	 * \return n_dec
	 *
	 */
	size_t get_ndec()
	{
		return n_dec;
	}

	/*! \brief Returns the option flags that control the Verlet list
	 *
	 *
	 * \return option flags
	 *
	 */
	size_t getOpt() const
	{
		return opt;
	}

	/*! \brief Sets the option flags that control the Verlet list
	 *
	 * \param opt option flags
	 *
	 */
	size_t setOpt(size_t opt)
	{
		this->opt = opt;
	}

	/*! \brief Return the domain particle sequence
	 *
	 * \return the particle sequence
	 *
	 */
	openfpm::vector<typename Mem_type::local_index_type> & getParticleSeq()
	{
		return domainParticlesCRS;
	}
};



#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_ */
