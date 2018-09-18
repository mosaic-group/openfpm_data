/*
 * VerletListFast.hpp
 *
 *  Created on: Aug 16, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTFAST_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTFAST_HPP_

#include "VerletNNIterator.hpp"
#include "NN/CellList/CellList_util.hpp"
#include "NN/Mem_type/MemFast.hpp"
#include "NN/Mem_type/MemBalanced.hpp"
#include "NN/Mem_type/MemMemoryWise.hpp"

#define VERLET_STARTING_NSLOT 128


#define WITH_RADIUS 3

/*! \brief Get the neighborhood iterator based on type
 *
 * \tparam dim dimensionality
 * \tparam T type of space
 * \tparam CellListImpl Cell-list implementation
 * \tparam type of neighborhood
 * \tparam PartIt Particle iterator
 *
 */
template<unsigned int dim, typename T, typename CellListImpl, typename PartIt, int type, typename local_index>
struct NNType
{
	/*! \brief Get the neighborhood
	 *
	 * \param v particle positions
	 * \param xp Position of the particle p
	 * \param p id of the particle p
	 * \param cl Cell-list type implementation
	 * \param r_cut Cutoff radius
	 *
	 * \return the NN iterator
	 *
	 */
	static inline auto get(const PartIt & it,
			               const openfpm::vector<Point<dim,T>> & v,
						   Point<dim,T> & xp,
						   size_t p,
						   CellListImpl & cl,
						   T r_cut) -> decltype(cl.template getNNIterator<NO_CHECK>(0))
	{
		return cl.template getNNIterator<NO_CHECK>(cl.getCell(xp));
	}

	/*! \brief Add particle in the list of the domain particles
	 *
	 * \param p particle id
	 * \param pc domain particle sequence
	 *
	 */
	static inline void add(local_index p, openfpm::vector<local_index> & pc)
	{
	}
};


/*! \brief Get the neighborhood iterator based on type
 *
 * specialization for the case with NN with radius
 *
 * \tparam dim dimensionality
 * \tparam T type of space
 * \tparam CellListImpl Cell-list implementation
 * \tparam PartIt Particle iterator
 *
 */
template<unsigned int dim, typename T, typename CellListImpl, typename PartIt, typename local_index>
struct NNType<dim,T,CellListImpl,PartIt,WITH_RADIUS,local_index>
{
	/*! \brief Get the neighborhood
	 *
	 * \param v particle position vector
	 * \param xp Position of the particle p
	 * \param p id of the particle p
	 * \param cl Cell-list type implementation
	 * \param r_cut Cutoff radius
	 *
	 * \return the NN iterator
	 *
	 */
	static inline auto get(const PartIt & it,
			               const openfpm::vector<Point<dim,T>> & v,
						   Point<dim,T> & xp,
						   size_t p,
						   CellListImpl & cl,
						   T r_cut) -> decltype(cl.template getNNIteratorRadius<NO_CHECK>(0,0.0))
	{
		return cl.template getNNIteratorRadius<NO_CHECK>(cl.getCell(xp),r_cut);
	}

	/*! \brief Add particle in the list of the domain particles
	 *
	 * \param p particle id
	 * \param pc domain particle sequence
	 *
	 */
	static inline void add(local_index p, openfpm::vector<local_index> & pc)
	{
	}
};


/*! \brief Get the neighborhood iterator based on type
 *
 * specialization for the case with NN symmetric
 *
 * \tparam dim dimensionality
 * \tparam T type of space
 * \tparam CellListImpl Cell-list implementation
 * \tparam PartIt particle iterator
 *
 */
template<unsigned int dim, typename T, typename CellListImpl, typename PartIt, typename local_index>
struct NNType<dim,T,CellListImpl,PartIt,VL_SYMMETRIC,local_index>
{
	/*! \brief Get the neighborhood
	 *
	 * \param v vector of the particles positions
	 * \param xp Position of the particle p
	 * \param p id of the particle p
	 * \param cl Cell-list type implementation
	 * \param r_cut Cutoff radius
	 *
	 * \return the NN iterator
	 *
	 */
	static inline auto get(const PartIt & it, const openfpm::vector<Point<dim,T>> & v, Point<dim,T> & xp, size_t p, CellListImpl & cl, T r_cut) -> decltype(cl.template getNNIteratorSym<NO_CHECK>(0,0,openfpm::vector<Point<dim,T>>()))
	{
		return cl.template getNNIteratorSym<NO_CHECK>(cl.getCell(xp),p,v);
	}

	/*! \brief Add particle in the list of the domain particles
	 *
	 * \param p particle id
	 * \param pc domain particle sequence
	 *
	 */
	static inline void add(local_index p, openfpm::vector<local_index> & pc)
	{
	}
};


/*! \brief Get the neighborhood iterator based on type
 *
 * specialization for the case with NN CRS symmetric
 *
 * \tparam dim dimensionality
 * \tparam T type of space
 * \tparam CellListImpl Cell-list implementation
 * \tparam PartIt Particle iterator
 *
 */
template<unsigned int dim, typename T, typename CellListImpl, typename PartIt, typename local_index>
struct NNType<dim,T,CellListImpl,PartIt,VL_CRS_SYMMETRIC,local_index>
{
	/*! \brief Get the neighborhood
	 *
	 * \param it particle iterator
	 * \param v vector of the particles positions
	 * \param xp Position of the particle p
	 * \param p id of the particle p
	 * \param cl Cell-list type implementation
	 * \param r_cut Cutoff radius
	 *
	 * \return the NN iterator
	 *
	 */
	static inline auto get(const PartIt & it,const openfpm::vector<Point<dim,T>> & v, Point<dim,T> & xp, size_t p, CellListImpl & cl, T r_cut) -> decltype(it.getNNIteratorCSR(v))
	{
		return it.getNNIteratorCSR(v);
	}

	/*! \brief Add particle in the list of the domain particles
	 *
	 * \param p particle id
	 * \param pc domain particle sequence
	 *
	 */
	static inline void add(local_index p, openfpm::vector<local_index> & pc)
	{
		pc.add(p);
	}
};

/*! \brief In general different NN scheme like full symmetric or CRS require different
 *         iterators over particles this class select the proper one
 *
 * This is the full or symmetric case
 *
 */
template<unsigned int type, unsigned int dim, typename vector, typename CellList>
class PartItNN
{
public:

	/*! \brief It return the particle iterator
	 *
	 * \param pos vector with position of the particles
	 * \param dom list of cells with normal neighborhood
	 * \param anom list of cells with not-normal neighborhood
	 * \param cli Cell-list used for Verlet-list construction
	 * \param g_m ghost marker
	 * \param end final particle id
	 *
	 * \return the particle iterator
	 *
	 */
	static inline auto get(const vector & pos, const openfpm::vector<size_t> & dom, const openfpm::vector<subsub_lin<dim>> & anom, CellList & cli, size_t g_m, size_t & end) -> decltype(pos.getIteratorTo(0))
	{
		end = g_m;
		return pos.getIteratorTo(end);
	}
};

/*! \brief In general different NN scheme like full symmetric or CRS require different
 *         iterators over particles this class select the proper one
 *
 * This is the CRS scheme
 *
 */
template<unsigned int dim, typename vector, typename CellList>
class PartItNN<VL_CRS_SYMMETRIC,dim,vector,CellList>
{
public:

	/*! \brief It return the particle iterator
	 *
	 * \param pos vector with position of the particles
	 * \param dom list of cells with normal neighborhood
	 * \param anom list of cells with not-normal neighborhood
	 * \param cli Cell-list used for Verlet-list construction
	 * \param g_m ghost marker
	 * \param end final particle id
	 *
	 * \return the particle iterator
	 *
	 */
	static inline ParticleItCRS_Cells<dim,CellList> get(const vector & pos, const openfpm::vector<size_t> & dom, const openfpm::vector<subsub_lin<dim>> & anom, CellList & cli, size_t g_m, size_t & end)
	{
		end = pos.size();
		return ParticleItCRS_Cells<dim,CellList>(cli,dom,anom,cli.getNNc_sym());
	}
};

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
		 typename CellListImpl = CellList<dim,T,Mem_fast<HeapMemory,typename Mem_type::local_index_type>,transform> >
class VerletList: public Mem_type
{
protected:

	//! Number of slot for each particle. Or maximum number of particles for each particle
	typename Mem_type::local_index_type slot;

	//! Domain particles
	openfpm::vector<typename Mem_type::local_index_type> dp;

private:

	//! decomposition counter
	size_t n_dec;

	//! Interlal cell-list
	CellListImpl cli;


	/*! \brief Fill the cell-list with data
	 *
	 * \param cli Cell-list
	 * \param pos vector of positions
	 * \param g_m marker
	 * \param opt VL_SYMMETRIC or VL_NON_SYMMETRIC
	 *
	 */
	void initCl(CellListImpl & cli, openfpm::vector<Point<dim,T>> & pos, size_t g_m, size_t opt)
	{
		if (opt & VL_SYMMETRIC || opt & VL_CRS_SYMMETRIC)
		{populate_cell_list(pos,cli,g_m,CL_SYMMETRIC);}
		else
		{populate_cell_list(pos,cli,g_m,CL_NON_SYMMETRIC);}
	}

	/*! \brief Create the Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param pos2 vector of positions of neighborhood particles
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cl Cell-list elements to use to construct the verlet list
	 * \param opt options to create the verlet list like VL_SYMMETRIC or VL_NON_SYMMETRIC
	 *
	 */
	inline void create(const openfpm::vector<Point<dim,T>> & pos, const openfpm::vector<Point<dim,T>> & pos2, const openfpm::vector<size_t> & dom, const openfpm::vector<subsub_lin<dim>> & anom, T r_cut, size_t g_m, CellListImpl & cl, size_t opt)
	{
		if (opt == VL_CRS_SYMMETRIC)
		{
			create_<CellNNIteratorSym<dim,CellListImpl,RUNTIME,NO_CHECK>,VL_CRS_SYMMETRIC>(pos,pos2,dom,anom,r_cut,g_m,cl,opt);
		}
		else if (opt == VL_SYMMETRIC)
		{
			create_<decltype(cl.template getNNIteratorSym<NO_CHECK>(0,0,pos)),VL_SYMMETRIC>(pos,pos2,dom,anom,r_cut,g_m,cl,opt);
		}
		else
		{
			create_<decltype(cl.template getNNIterator<NO_CHECK>(0)),VL_NON_SYMMETRIC>(pos,pos2,dom,anom,r_cut,g_m,cl,opt);
		}
	}

	/*! \brief Create the Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param pos2 vector of position for the neighborhood
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cli Cell-list elements to use to construct the verlet list
	 * \param dom list of domain cells with normal neighborhood
	 * \param anom list of domain cells with non-normal neighborhood
	 * \param opt options
	 *
	 */
	template<typename NN_type, int type> inline void create_(const openfpm::vector<Point<dim,T>> & pos, const openfpm::vector<Point<dim,T>> & pos2 , const openfpm::vector<size_t> & dom, const openfpm::vector<subsub_lin<dim>> & anom, T r_cut, size_t g_m, CellListImpl & cli, size_t opt)
	{
		size_t end;

		auto it = PartItNN<type,dim,openfpm::vector<Point<dim,T>>,CellListImpl>::get(pos,dom,anom,cli,g_m,end);

		Mem_type::init_to_zero(slot,end);

		dp.clear();

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		while (it.isNext())
		{
			typename Mem_type::local_index_type i = it.get();
			Point<dim,T> xp = pos.template get<0>(i);

			// Get the neighborhood of the particle
			NN_type NN = NNType<dim,T,CellListImpl,decltype(it),type,typename Mem_type::local_index_type>::get(it,pos,xp,i,cli,r_cut);
			NNType<dim,T,CellListImpl,decltype(it),type,typename Mem_type::local_index_type>::add(i,dp);

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

		Mem_type::init_to_zero(slot,g_m);

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

public:

	//! type for the local index
	typedef Mem_type Mem_type_type;

	//! Object type that the structure store
	typedef size_t value_type;

	//! CellList implementation used for Verlet list construction
	typedef CellListImpl CellListImpl_;

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
		initCl(cli,pos,g_m,opt);

		// Unuseful empty vector
		openfpm::vector<subsub_lin<dim>> anom_c;
		openfpm::vector<size_t> dom_c;

		// create verlet
		create(pos, pos,dom_c,anom_c,r_cut,g_m,cli,opt);
	}

	/*! \brief Initialize the symmetric Verlet-list
	 *
	 * \param box Simulation domain
	 * \param dom Processor domain
	 * \param g ghost size
	 * \param r_cut cut-off radius
	 * \param pos vector of particle positions
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	void InitializeSym(const Box<dim,T> & box, const Box<dim,T> & dom, const Ghost<dim,T> & g, T r_cut, openfpm::vector<Point<dim,T>> & pos, size_t g_m)
	{
		// Padding
		size_t pad = 0;

		// Cell decomposer
		CellDecomposer_sm<dim,T,shift<dim,T>> cd_sm;

		// Calculate the divisions for the Cell-lists
		cl_param_calculateSym<dim,T>(box,cd_sm,g,r_cut,pad);

		// Initialize a cell-list
		cli.Initialize(cd_sm,dom,pad);
		initCl(cli,pos,g_m,VL_SYMMETRIC);

		// Unused
		openfpm::vector<subsub_lin<dim>> anom_c;
		openfpm::vector<size_t> dom_c;

		// create verlet
		create(pos, pos,dom_c,anom_c,r_cut,g_m,cli,VL_SYMMETRIC);
	}


	/*! \brief Initialize the symmetric Verlet-list CRS scheme
	 *
	 * \param box Simulation domain
	 * \param dom Processor domain
	 * \param g ghost size
	 * \param r_cut cut-off radius
	 * \param pos vector of particle positions
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	void InitializeCrs(const Box<dim,T> & box, const Box<dim,T> & dom, const Ghost<dim,T> & g, T r_cut, openfpm::vector<Point<dim,T>> & pos, size_t g_m)
	{
		// Padding
		size_t pad = 0;

		// Cell decomposer
		CellDecomposer_sm<dim,T,shift<dim,T>> cd_sm;

		// Calculate the divisions for the Cell-lists
		cl_param_calculateSym<dim,T>(box,cd_sm,g,r_cut,pad);

		// Initialize a cell-list
		cli.Initialize(cd_sm,dom,pad);
		initCl(cli,pos,g_m,VL_SYMMETRIC);
	}

	/*! \brief Create the Verlet-list with the crossing scheme
	 *
	 * \param pos vector with the particle positions
	 * \param g_m ghost marker
	 * \param pos vector with the particle positions
	 * \param r_cut cut-off radius
	 * \param dom_c domain cells
	 * \param anom_c cells with anomalos neighborhood
	 *
	 */
	void createVerletCrs(T r_cut, size_t g_m, openfpm::vector<Point<dim,T>> & pos, openfpm::vector<size_t> & dom_c, openfpm::vector<subsub_lin<dim>> & anom_c)
	{
		// create verlet
		create(pos, pos,dom_c,anom_c,r_cut,g_m,cli,VL_CRS_SYMMETRIC);
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
		initCl(cli,pos,g_m,opt);

		// Unused
		openfpm::vector<subsub_lin<dim>> anom_c;
		openfpm::vector<size_t> dom_c;

		create(pos, pos,dom_c,anom_c,r_cut,g_m,cli,opt);
	}

	/*! \brief update the Verlet list
	 *
	 * \param r_cut cutoff radius
	 * \param dom Processor domain
	 * \param pos vector of particle positions
	 * \param g_m ghost marker
	 * \param dom_c list of cells with normal neighborhood
	 * \param anom_c list of cells with anormal neighborhood
	 *
	 */
	void updateCrs(const Box<dim,T> & dom, T r_cut, openfpm::vector<Point<dim,T>> & pos, size_t & g_m, const openfpm::vector<size_t> & dom_c, const openfpm::vector<subsub_lin<dim>> & anom_c)
	{
		initCl(cli,pos,g_m,VL_CRS_SYMMETRIC);

		create(pos,pos,dom_c,anom_c,r_cut,g_m,cli,VL_CRS_SYMMETRIC);
	}

	/*! Initialize the verlet list from an already filled cell-list
	 *
	 * \param cli external Cell-list
	 * \param r_cut cutoff-radius
	 * \param pos vector of particle positions
	 * \param pos2 vector of particle position for the neighborhood
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * 	\param opt options for the Verlet-list creation
	 *
	 */
	void Initialize(CellListImpl & cli,
					T r_cut,
					const openfpm::vector<Point<dim,T>> & pos,
					const openfpm::vector<Point<dim,T>> & pos2,
					size_t g_m, size_t opt = VL_NON_SYMMETRIC)
	{
		Point<dim,T> spacing = cli.getCellBox().getP2();

		// Create with radius or not
		bool wr = true;

		for (size_t i = 0 ; i < dim ; i++)
			wr &= r_cut <= spacing.get(i);

		if (wr == true || opt == VL_SYMMETRIC)
		{
			openfpm::vector<subsub_lin<dim>> anom_c;
			openfpm::vector<size_t> dom_c;

			create(pos,pos2,dom_c,anom_c,r_cut,g_m,cli,opt);
		}
		else
		{
			openfpm::vector<subsub_lin<dim>> anom_c;
			openfpm::vector<size_t> dom_c;

			create_<decltype(cli.template getNNIteratorRadius<NO_CHECK>(0,0.0)),WITH_RADIUS>(pos,pos2,dom_c,anom_c,r_cut,g_m,cli,VL_NON_SYMMETRIC);
		}
	}

	//! Default Constructor
	VerletList()
	:Mem_type(VERLET_STARTING_NSLOT),slot(VERLET_STARTING_NSLOT),n_dec(0)
	{};

	//! Copy constructor
	VerletList(const VerletList<dim,T,Mem_type,transform,CellListImpl> & cell)
	:Mem_type(VERLET_STARTING_NSLOT),slot(VERLET_STARTING_NSLOT)
	{
		this->operator=(cell);
	}

	//! Copy constructor
	VerletList(VerletList<dim,T,Mem_type,transform,CellListImpl> && cell)
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
	VerletList<dim,T,Mem_type,transform,CellListImpl> &
	operator=(VerletList<dim,T,Mem_type,transform,CellListImpl> && vl)
	{
		slot = vl.slot;

		Mem_type::operator=(vl);
		dp.swap(vl.dp);

		n_dec = vl.n_dec;

		return *this;
	}

	/*! \brief Copy a verlet list
	 *
	 * \param vl verlet-list to copy
	 *
	 * \return itself
	 *
	 */
	VerletList<dim,T,Mem_type,transform,CellListImpl> & operator=(const VerletList<dim,T,Mem_type,transform,CellListImpl> & vl)
	{
		slot = vl.slot;

		Mem_type::operator=(vl);

		cli = vl.cli;

		dp = vl.dp;
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
	inline void swap(VerletList<dim,T,Mem_type,transform,CellListImpl> & vl)
	{
		Mem_type::swap(vl);
		dp.swap(vl.dp);

		size_t vl_slot_tmp = vl.slot;
		vl.slot = slot;
		slot = vl_slot_tmp;

		cli.swap(vl.cli);

		size_t n_dec_tmp = vl.n_dec;
		vl.n_dec = n_dec;
		n_dec = n_dec_tmp;
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
	template<unsigned int impl=NO_CHECK>
	inline VerletNNIterator<dim,VerletList<dim,T,Mem_type,transform,CellListImpl>>
	getNNIterator(size_t part_id)
	{
		VerletNNIterator<dim,VerletList<dim,T,Mem_type,transform,CellListImpl>> vln(part_id,*this);

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

	/*! \brief Return the domain particle sequence
	 *
	 * \return the particle sequence
	 *
	 */
	openfpm::vector<typename Mem_type::local_index_type> & getParticleSeq()
	{
		return dp;
	}
};



#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTFAST_HPP_ */
