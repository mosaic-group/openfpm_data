#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_

#include "Vector/map_vector.hpp"

#include "VerletNNIterator.hpp"
#include "NN/CellList/CellList.hpp"
#include "NN/CellList/CellList_util.hpp"
#include "NN/Mem_type/MemFast.hpp"
#include "NN/Mem_type/MemBalanced.hpp"
#include "NN/Mem_type/MemMemoryWise.hpp"


// Verlet list config options
constexpr int VL_NON_SYMMETRIC = 1;
constexpr int VL_SYMMETRIC = 2;
constexpr int VL_CRS_SYMMETRIC = 4;
constexpr int VL_ADAPTIVE_RCUT = 8;
constexpr int VL_NMAX_NEIGHBOR = 16;
constexpr int VL_SKIP_REF_PART = 32;

constexpr int VERLET_STARTING_NSLOT = 128;

#ifdef LOCAL_INDEX64
typedef size_t local_index_;
#else
typedef unsigned int local_index_;
#endif


/*! Functor class for Verlet list particle neighborhood iteration on initialization
 * Partial template specializations implement different methods of initialization
 *
 * isNMax = maximum number of neighbors when VL_NMAX_NEIGHBOR is enabled
 * skipRefPart = skip reference particle when VL_SKIP_REF_PART is enabled
 *
*/

template<bool isNMax, bool skipRefPart>
struct iteratePartNeighbor;

/*! Functor class for Verlet list particle neighborhood iteration on initialization
 * Partial template specializations implement different methods of initialization
 *
 * isNMax = maximum number of neighbors when VL_NMAX_NEIGHBOR is enabled
 * skipRefPart = skip reference particle when VL_SKIP_REF_PART is enabled
 *
 * \tparam verletList Verlet List to be filled
 * \tparam it particle neighborhood iterator
 * \tparam pos position of reference particle p
 * \tparam p index (id) of reference particle p
 * \tparam xp position of reference particle p
 * \tparam r_cut cut-off radius for particle interaction
 *
 *
 */
template<>
struct iteratePartNeighbor<false, false> {

	template<typename NeighborIter_type, typename VerletList_type, typename vPos_type, unsigned int dim, typename T>
	inline void operator() (
		VerletList_type& verletList,
		NeighborIter_type & it,
		const vPos_type & pos,
		size_t p, Point<dim,T> xp,
		T r_cut, size_t neighborMaxNum)
	{
		T r_cut2 = r_cut * r_cut;

		while (it.isNext())
		{
			auto q = it.get();

			Point<dim,T> xq = pos.template get<0>(q);

			if (xp.distance2(xq) < r_cut2)
				verletList.addPart(p,q);

			++it;
		}
	}
};

/*! Functor class for Verlet list particle neighborhood iteration on initialization
 * Partial template specializations implement different methods of initialization
 *
 * isNMax = maximum number of neighbors when VL_NMAX_NEIGHBOR is enabled
 * skipRefPart = skip reference particle when VL_SKIP_REF_PART is enabled
 *
 * \tparam verletList Verlet List to be filled
 * \tparam it particle neighborhood iterator
 * \tparam pos position of reference particle p
 * \tparam p index (id) of reference particle p
 * \tparam xp position of reference particle p
 * \tparam r_cut cut-off radius for particle interaction
*
 *
 */
template<>
struct iteratePartNeighbor<false, true> {

	template<typename NeighborIter_type, typename VerletList_type, typename vPos_type, unsigned int dim, typename T>
	inline void operator() (
		VerletList_type& verletList,
		NeighborIter_type & it,
		const vPos_type & pos,
		size_t p, Point<dim,T> xp,
		T r_cut, size_t neighborMaxNum)
	{
		T r_cut2 = r_cut * r_cut;

		while (it.isNext())
		{
			auto q = it.get();

			Point<dim,T> xq = pos.template get<0>(q);

			if (xp.distance2(xq) < r_cut2 && p != q)
				verletList.addPart(p,q);
			++it;
		}
	}
};

/*! Functor class for Verlet list particle neighborhood iteration on initialization
 * Partial template specializations implement different methods of initialization
 *
 * isNMax = maximum number of neighbors when VL_NMAX_NEIGHBOR is enabled
 * skipRefPart = skip reference particle when VL_SKIP_REF_PART is enabled
 *
 * \tparam verletList Verlet List to be filled
 * \tparam it particle neighborhood iterator
 * \tparam pos position of reference particle p
 * \tparam p index (id) of reference particle p
 * \tparam xp position of reference particle p
 * \tparam r_cut cut-off radius for particle interaction
 * \tparam neighborMaxNum maximum numberof neighbors
 *
 *
 */
template<>
struct iteratePartNeighbor<true, false> {
	template<typename NeighborIter_type, typename VerletList_type, typename vPos_type, unsigned int dim, typename T>
	inline void operator() (
		VerletList_type& verletList,
		NeighborIter_type & it,
		const vPos_type & pos,
		size_t p, Point<dim,T> xp,
		T r_cut, size_t neighborMaxNum)
	{
		struct ReorderType {
			T dist;
			size_t id;
			bool operator<(const ReorderType &other) const { return this->dist < other.dist; }
		};

		openfpm::vector<ReorderType> neighborList;

		while (it.isNext())
		{
			auto q = it.get();

			Point<dim,T> xq = pos.template get<0>(q);

			if (xp.distance(xq) < r_cut) {
				ReorderType qReorder;
				qReorder.dist = xp.distance(xq);
				qReorder.id = q;
				neighborList.add(qReorder);
			}

			++it;
		}

		neighborList.sort();
		for(size_t i = 0; i < neighborList.size(); ++i) {
			if (i < neighborMaxNum)
			verletList.addPart(p, neighborList.get(i).id);
		}
	}
};

/*! Functor class for Verlet list particle neighborhood iteration on initialization
 * Partial template specializations implement different methods of initialization
 *
 * isNMax = maximum number of neighbors when VL_NMAX_NEIGHBOR is enabled
 * skipRefPart = skip reference particle when VL_SKIP_REF_PART is enabled
 *
 * \tparam verletList Verlet List to be filled
 * \tparam it particle neighborhood iterator
 * \tparam pos position of reference particle p
 * \tparam p index (id) of reference particle p
 * \tparam xp position of reference particle p
 * \tparam r_cut cut-off radius for particle interaction
 * \tparam neighborMaxNum maximum number of neighbors
 *
 *
 */
template<>
struct iteratePartNeighbor<true, true> {

	template<typename NeighborIter_type, typename VerletList_type, typename vPos_type, unsigned int dim, typename T>
	inline void operator() (
		VerletList_type& verletList,
		NeighborIter_type & it,
		const vPos_type & pos,
		size_t p, Point<dim,T> xp,
		T r_cut, size_t neighborMaxNum)
	{
		struct ReorderType {
			T dist;
			size_t id;
			bool operator<(const ReorderType &other) const { return this->dist < other.dist; }
		};

		openfpm::vector<ReorderType> neighborList;

		while (it.isNext())
		{
			auto q = it.get();

			Point<dim,T> xq = pos.template get<0>(q);

			if (xp.distance(xq) < r_cut) {
				ReorderType qReorder;
				qReorder.dist = xp.distance(xq);
				qReorder.id = q;
				neighborList.add(qReorder);
			}

			++it;
		}

		neighborList.sort();
		for(size_t i = 0; i < neighborList.size(); ++i) {
			if (i < neighborMaxNum && p != neighborList.get(i).id)
				verletList.addPart(p, neighborList.get(i).id);
		}
	}
};


/*! \brief Class for Verlet list implementation
 *
 * * M = number of particles
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
	unsigned int opt = VL_NON_SYMMETRIC,
	typename Mem_type = Mem_fast<HeapMemory,local_index_>,
	typename transform = no_transform<dim,T>,
	typename vPos_type = openfpm::vector<Point<dim,T>>,
	typename CellListImpl = CellList<dim,T,Mem_fast<HeapMemory,typename Mem_type::local_index_type>,transform,vPos_type> >
class VerletList: public Mem_type
{
protected:

	//! Number of slot for each particle. Or maximum number of particles for each particle
	typename Mem_type::local_index_type slot;

	//! Domain particles
	openfpm::vector<typename Mem_type::local_index_type> domainParticlesCRS;

	//! Max number of Nieghbors. Only with opt |= VL_NMAX_NEIGHBOR
	size_t neighborMaxNum=0;

private:

	//! decomposition counter
	size_t n_dec;

	//! Internal cell-list
	CellListImpl cellList;


	/*! \brief Fill the cell-list with data
	 *
	 * \param cellList Cell-list
	 * \param vPos vector of positions
	 * \param ghostMarker marker
	 * \param opt VL_SYMMETRIC or VL_NON_SYMMETRIC
	 *
	 */
	void initCl(CellListImpl & cellList, vPos_type & vPos, size_t ghostMarker)
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


public:

	/*! \brief Fill CRS Symmetric Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cellList Cell-list elements to use to construct the verlet list
	 * \param dom list of domain cells with normal neighborhood
	 * \param anom list of domain cells with non-normal neighborhood
	 *
	 */
	inline void fillCRSSymmetric(
		const vPos_type & pos,
		T r_cut,
		size_t ghostMarker,
		const openfpm::vector<size_t> & dom,
		const openfpm::vector<subsub_lin<dim>>& anom)
	{
		size_t end = pos.size();

		Mem_type::init_to_zero(slot,end);
		domainParticlesCRS.clear();

		// iterate the particles
		auto it = ParticleItCRS_Cells<dim,CellListImpl,vPos_type>(cellList,dom,anom,cellList.getNNc_sym());
		while (it.isNext())
		{
			typename Mem_type::local_index_type p = it.get();
			Point<dim,T> xp = pos.template get<0>(p);

			domainParticlesCRS.add(p);

			// Get the neighborhood of the particle
			auto NN = it.getNNIteratorCSR(pos);

			iteratePartNeighbor<opt&VL_NMAX_NEIGHBOR,opt&VL_SKIP_REF_PART>{}(*this, NN, pos, p, xp, r_cut, neighborMaxNum);
			++it;
		}
	}

	/*! \brief Fill Symmetric Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cellList Cell-list elements to use to construct the verlet list
	 *
	 */
	inline void fillSymmetric(
		const vPos_type & pos,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cellList)
	{
		size_t end = ghostMarker;

		Mem_type::init_to_zero(slot,end);
		domainParticlesCRS.clear();

		// iterate the particles
		auto it = pos.getIteratorTo(end);
		while (it.isNext())
		{
			typename Mem_type::local_index_type p = it.get();
			Point<dim,T> xp = pos.template get<0>(p);

			// Get the neighborhood of the particle
			auto NN = cellList.getNNIteratorBoxSym(cellList.getCell(xp),p,pos);

			iteratePartNeighbor<opt&VL_NMAX_NEIGHBOR,opt&VL_SKIP_REF_PART>{}(*this, NN, pos, p, xp, r_cut, neighborMaxNum);
			++it;
		}
	}

	/*! \brief Fill Non-symmetric Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cellList Cell-list elements to use to construct the verlet list
	 *
	 */
	inline void fillNonSymmetric(
		const vPos_type & pos,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cellList)
	{
		size_t end = ghostMarker;

		Mem_type::init_to_zero(slot,end);

		// iterate the particles
		auto it = pos.getIteratorTo(end);
		while (it.isNext())
		{
			typename Mem_type::local_index_type p = it.get();
			Point<dim,T> xp = pos.template get<0>(p);

			// Get the neighborhood of the particle
			auto NN = cellList.getNNIteratorBox(cellList.getCell(xp));

			iteratePartNeighbor<opt&VL_NMAX_NEIGHBOR,opt&VL_SKIP_REF_PART>{}(*this, NN, pos, p, xp, r_cut, neighborMaxNum);
			++it;
		}
	}

	/*! \brief Fill Non-symmetric Verlet list from a given cell-list
	 *
	 * \param it domain iterator of particle ids
	 * \param pos vector of particle positions
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cellList Cell-list elements to use to construct the verlet list
	 *
	 */
	template <typename domainIterator_type>
	inline void fillNonSymmetricIterator(
		domainIterator_type& it,
		const vPos_type & pos,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cellList)
	{
		size_t end = ghostMarker;

		Mem_type::init_to_zero(slot,end);

		while (it.isNext())
		{
			typename Mem_type::local_index_type p = it.get();
			Point<dim,T> xp = pos.template get<0>(p);

			// Get the neighborhood of the particle
			auto NN = cellList.getNNIteratorBox(cellList.getCell(xp));

			iteratePartNeighbor<opt&VL_NMAX_NEIGHBOR,opt&VL_SKIP_REF_PART>{}(*this, NN, pos, p, xp, r_cut, neighborMaxNum);
			++it;
		}
	}

	/*! \brief Fill non-symmetric adaptive r-cut Verlet list from a list of cut-off radii
	 *
	 * \param pos vector of positions
	 * \param rCuts list of cut-off radii for every particle in pos
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	inline void fillNonSymmAdaptive(
		const vPos_type & pos,
		openfpm::vector<T> &rCuts,
		size_t ghostMarker)
	{
		if (rCuts.size() != ghostMarker)
		{
			std::cerr << __FILE__ << ":" << __LINE__
				<< " ERROR: when constructing adaptive cut-off Verlet list, pos.size_local() != rCuts.size(), ["
				<< rCuts.size() << "!=" << pos.size_local() << "]" << std::endl;
			std::runtime_error("Runtime adaptive cut-off Verlet list error");
		}

		size_t end = ghostMarker;

		Mem_type::init_to_zero(slot,end);

		// iterate the particles
		auto it = pos.getIteratorTo(end);
		while (it.isNext())
		{
			typename Mem_type::local_index_type p = it.get();
			Point<dim,T> xp = pos.template get<0>(p);

			// iterate the whole domain instead of neighborhood iteration
			// no auxillary cell list is needed
			auto NN = pos.getIteratorTo(pos.size_local());
			T r_cut = rCuts.get(p);

			iteratePartNeighbor<opt&VL_NMAX_NEIGHBOR,opt&VL_SKIP_REF_PART>{}(*this, NN, pos, p, xp, r_cut, neighborMaxNum);
			++it;
		}
	}

	/*! \brief Fill Non-symmetric Verlet list from a given cell-list (Radius)
	 *
	 * \param pos vector of positions
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cellList Cell-list elements to use to construct the verlet list
	 *
	 */
	inline void fillNonSymmetricRadius(
		const vPos_type & pos,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cellList)
	{
		size_t end = ghostMarker;

		Mem_type::init_to_zero(slot,end);

		// iterate the particles
		auto it = pos.getIteratorTo(end);
		while (it.isNext())
		{
			typename Mem_type::local_index_type p = it.get();
			Point<dim,T> xp = pos.template get<0>(p);

			// Get the neighborhood of the particle
			auto NN = cellList.getNNIteratorRadius(cellList.getCell(xp),r_cut);

			iteratePartNeighbor<opt&VL_NMAX_NEIGHBOR,opt&VL_SKIP_REF_PART>{}(*this, NN, pos, p, xp, r_cut, neighborMaxNum);
			++it;
		}
	}

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

	/*! \brief Replace the neighborhood particles for the particle id part_id with the given buffer
	 *
	 * \param part_id id of the particle
	 * \param buffer used to update the neighborhood
	 *
	 */
	inline void replaceNeighbParts(size_t part_id, const openfpm::vector<size_t> &buffer)
	{
		Mem_type::clear(part_id);
		for (size_t i = 0; i < buffer.size(); i++)
		{
			Mem_type::addCell(part_id,buffer.get(i));
		}
	}

	/*! Initialize the verlet list for Non-Symmetric case
	 *
	 * \param box Domain where this cell list is living
	 * \param r_cut cut-off radius
	 * \param pos vector of particle positions
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	void Initialize(
		const Box<dim,T> & box,
		T r_cut,
		vPos_type & pos,
		size_t ghostMarker)
	{
		// Number of divisions
		size_t div[dim];

		Box<dim,T> bt = box;

		// Calculate the divisions for the Cell-lists
		cl_param_calculate(bt,div,r_cut,Ghost<dim,T>(0.0));

		// Initialize a cell-list
		cellList.Initialize(bt,div);

		initCl(cellList,pos,ghostMarker);

		fillNonSymmetric(pos,r_cut,ghostMarker,cellList);
	}

	/*! Initialize the verlet list for Non-Symmetric case from an already filled cell-list
	 *
	 * \param cellList external Cell-list
	 * \param r_cut cutoff-radius
	 * \param it domain iterator of particle ids
	 * \param pos vector of particle positions
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	template <typename domainIterator_type>
	void Initialize(
		CellListImpl& cellList,
		T r_cut,
		domainIterator_type& it,
		const vPos_type& pos,
		size_t ghostMarker)
	{
		Point<dim,T> spacing = cellList.getCellBox().getP2();

		// Create with radius or not
		bool wr = true;

		for (size_t i = 0 ; i < dim ; i++)
			wr &= r_cut <= spacing.get(i);

		fillNonSymmetricIterator(it,pos,r_cut,ghostMarker,cellList);
	}

	/*! \brief Initialize the symmetric Verlet-list
	 *
	 * \param box Simulation domain
	 * \param dom Processor domain
	 * \param ghostSize ghost size
	 * \param r_cut cut-off radius
	 * \param pos vector of particle positions
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	void InitializeSym(
		const Box<dim,T> & box,
		const Box<dim,T> & dom,
		const Ghost<dim,T> & ghostSize,
		T r_cut,
		vPos_type & pos,
		size_t ghostMarker)
	{
		// Padding
		size_t pad = 0;

		// Cell decomposer
		CellDecomposer_sm<dim,T,shift<dim,T>> cd_sm;

		// Calculate the divisions for the Cell-lists
		cl_param_calculateSym<dim,T>(box,cd_sm,ghostSize,r_cut,pad);

		// Initialize a cell-list
		cellList.Initialize(cd_sm,dom,pad);
		initCl(cellList,pos,ghostMarker);

		// create verlet
		fillSymmetric(pos, r_cut, ghostMarker, cellList);
	}


	/*! \brief Initialize the symmetric Verlet-list CRS scheme
	 *
	 * \param box Simulation domain
	 * \param dom Processor domain
	 * \param ghostSize ghost size
	 * \param r_cut cut-off radius
	 * \param pos vector of particle positions
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 *
	 */
	void initializeCrs(
		const Box<dim,T> & box,
		const Box<dim,T> & dom,
		const Ghost<dim,T> & ghostSize,
		T r_cut,
		vPos_type & pos,
		size_t ghostMarker)
	{
		// Padding
		size_t pad = 0;

		// Cell decomposer
		CellDecomposer_sm<dim,T,shift<dim,T>> cd_sm;

		// Calculate the divisions for the Cell-lists
		cl_param_calculateSym<dim,T>(box,cd_sm,ghostSize,r_cut,pad);

		// Initialize a cell-list
		cellList.Initialize(cd_sm,dom,pad);
		initCl(cellList,pos,ghostMarker);
	}

	/*! \brief update the Verlet list
	 *
	 * \param r_cut cutoff radius
	 * \param dom Processor domain
	 * \param pos vector of particle positions
	 * \param ghostMarker ghost marker
	 *
	 */
	void update(
		const Box<dim,T> & dom,
		T r_cut,
		vPos_type & pos,
		size_t & ghostMarker)
	{
		initCl(cellList,pos,ghostMarker);

		if (opt & VL_SYMMETRIC)
			fillSymmetric(pos,r_cut,ghostMarker,cellList);
		else
			fillNonSymmetric(pos,r_cut,ghostMarker,cellList);
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
		vPos_type & pos,
		T r_cut,
		size_t & ghostMarker,
		const openfpm::vector<size_t> & dom_c,
		const openfpm::vector<subsub_lin<dim>> & anom_c)
	{
		initCl(cellList,pos,ghostMarker);
		fillCRSSymmetric(pos,r_cut,ghostMarker,dom_c,anom_c);
	}

	//! Default Constructor
	VerletList()
	:Mem_type(VERLET_STARTING_NSLOT),slot(VERLET_STARTING_NSLOT),n_dec(0)
	{};

	//! Copy constructor
	VerletList(const VerletList<dim,T,opt,Mem_type,transform,vPos_type,CellListImpl> & cell)
	:Mem_type(VERLET_STARTING_NSLOT),slot(VERLET_STARTING_NSLOT)
	{
		this->operator=(cell);
	}

	//! Copy constructor
	VerletList(VerletList<dim,T,opt,Mem_type,transform,vPos_type,CellListImpl> && cell)
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
	VerletList(Box<dim,T> & box, T r_cut, vPos_type & pos, size_t ghostMarker, size_t slot=VERLET_STARTING_NSLOT)
	:slot(slot)
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
	VerletList(Box<dim,T> & box, Box<dim,T> & dom, T r_cut, vPos_type & pos, size_t ghostMarker, size_t slot=VERLET_STARTING_NSLOT)
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
	VerletList<dim,T,opt,Mem_type,transform,vPos_type,CellListImpl> &
	operator=(VerletList<dim,T,opt,Mem_type,transform,vPos_type,CellListImpl> && vl)
	{
		slot = vl.slot;

		Mem_type::operator=(vl);
		domainParticlesCRS.swap(vl.domainParticlesCRS);

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
	VerletList<dim,T,opt,Mem_type,transform,vPos_type,CellListImpl> &
	operator=(const VerletList<dim,T,opt,Mem_type,transform,vPos_type,CellListImpl> & vl)
	{
		slot = vl.slot;

		Mem_type::operator=(vl);

		cellList = vl.cellList;

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
	inline void swap(VerletList<dim,T,opt,Mem_type,transform,vPos_type,CellListImpl> & vl)
	{
		Mem_type::swap(vl);
		domainParticlesCRS.swap(vl.domainParticlesCRS);

		size_t vl_slot_tmp = vl.slot;
		vl.slot = slot;
		slot = vl_slot_tmp;

		cellList.swap(vl.cellList);

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
	inline VerletNNIterator<dim,VerletList<dim,T,opt,Mem_type,transform,vPos_type,CellListImpl>>
	getNNIterator(size_t part_id)
	{
		VerletNNIterator<dim,VerletList<dim,T,opt,Mem_type,transform,vPos_type,CellListImpl>> vln(part_id,*this);

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
		return cellList;
	}

	/*! \brief Set the n_dec number
	 *
	 * \param n_dec
	 *
	 */
	void set_ndec(size_t n_dec)
	{
		this->n_dec = n_dec;

		cellList.set_ndec(n_dec);
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

	/*! \brief Returns the max number of neighbors per particle. Only with opt | VL_NMAX_NEIGHBOR
	 *
	 *
	 * \return max number of neighbors per particle
	 *
	 */
	size_t getNeighborMaxNum() const
	{
		return neighborMaxNum;
	}

	/*! \brief Sets the max number of neighbors per particle. Only with opt | VL_NMAX_NEIGHBOR
	 *
	 *
	 * \return max number of neighbors per particle
	 *
	 */
	void setNeighborMaxNum(size_t neighborMaxNum)
	{
		this->neighborMaxNum = neighborMaxNum;
	}

	/*! \brief Returns the option flag template parameter opt that controls the Verlet list
	 *
	 *
	 * \return option flags
	 *
	 */
	size_t getOpt() const
	{
		return opt;
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

template<unsigned int dim, typename St, unsigned int opt> using VERLET_MEMFAST = VerletList<dim,St,opt,Mem_fast<>,shift<dim,St>>;
template<unsigned int dim, typename St, unsigned int opt> using VERLET_MEMBAL = VerletList<dim,St,opt,Mem_bal<>,shift<dim,St>>;
template<unsigned int dim, typename St, unsigned int opt> using VERLET_MEMMW = VerletList<dim,St,opt,Mem_mw<>,shift<dim,St>>;

template<unsigned int dim, typename St, unsigned int opt> using VERLET_MEMFAST_INT = VerletList<dim,St,opt,Mem_fast<HeapMemory,unsigned int>,shift<dim,St>>;
template<unsigned int dim, typename St, unsigned int opt> using VERLET_MEMBAL_INT = VerletList<dim,St,opt,Mem_bal<unsigned int>,shift<dim,St>>;
template<unsigned int dim, typename St, unsigned int opt> using VERLET_MEMMW_INT = VerletList<dim,St,opt,Mem_mw<unsigned int>,shift<dim,St>>;


#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLIST_HPP_ */
