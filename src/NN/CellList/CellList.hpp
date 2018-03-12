/*
 * CellList.hpp
 *
 *  Created on: Mar 21, 2015
 *      Author: Pietro Incardona
 */

#ifndef CELLLIST_HPP_
#define CELLLIST_HPP_

#include "Vector/map_vector.hpp"
#include "CellDecomposer.hpp"
#include "Space/SpaceBox.hpp"
#include "util/mathutil.hpp"
#include "CellNNIterator.hpp"
#include "Space/Shape/HyperCube.hpp"
#include "CellListNNIteratorRadius.hpp"
#include <unordered_map>

#include "CellListIterator.hpp"
#include "ParticleIt_Cells.hpp"
#include "ParticleItCRS_Cells.hpp"
#include "util/common.hpp"

#include "NN/Mem_type/MemFast.hpp"
#include "NN/Mem_type/MemBalanced.hpp"
#include "NN/Mem_type/MemMemoryWise.hpp"
#include "NN/CellList/NNc_array.hpp"

//! Wrapper of the unordered map
template<typename key,typename val>
class wrap_unordered_map: public std::unordered_map<key,val>
{
};

#ifdef HAVE_LIBQUADMATH

#include <boost/multiprecision/float128.hpp>


//! Wrapper of the unordered map
template<typename val>
class wrap_unordered_map<boost::multiprecision::float128,val>
{
};

#endif

//! Point at witch the cell do a reallocation (it should the the maximum for all the implementations)
#define CELL_REALLOC 16ul

#define STARTING_NSLOT 16

/*! \brief Calculate the the Neighborhood for symmetric interactions CSR scheme
 *
 * \param cNN calculated cross neighborhood
 * \param div Number of divisions in each direction
 *
 */
template<unsigned int dim> void NNcalc_csr(openfpm::vector<std::pair<grid_key_dx<dim>,grid_key_dx<dim>>> & cNN)
{
	// Calculate the NNc_full array, it is a structure to get the neighborhood array

	// compile-time array {0,0,0,....}  {2,2,2,...} {1,1,1,...}

	typedef typename generate_array<size_t,dim, Fill_zero>::result NNzero;
	typedef typename generate_array<size_t,dim, Fill_two>::result NNtwo;

	// Generate the sub-grid iterator

	size_t div[dim];

	// Calculate the divisions

	for (size_t i = 0 ; i < dim ; i++)
		div[i] = 4;

	grid_sm<dim,void> gs(div);
	grid_key_dx_iterator_sub<dim> gr_sub3(gs,NNzero::data,NNtwo::data);

	grid_key_dx<dim> src_; // Source cell
	for (size_t i = 0; i < dim; i++)
		src_.set_d(i,1);

	size_t middle = gs.LinId(src_);

	// Calculate the symmetric crs array
	while (gr_sub3.isNext())
	{
		auto dst = gr_sub3.get();
		grid_key_dx<dim> src = src_;

		if ((long int)middle > gs.LinId(dst))
		{
			++gr_sub3;
			continue;
		}

		// Here we adjust src and dst to be in the positive quadrant

		for (size_t i = 0 ; i < dim; i++)
		{
			if (dst.get(i) == 0)
			{
				src.set_d(i,src.get(i) + 1);
				dst.set_d(i,dst.get(i) + 1);
			}
		}

		src -= src_;
		dst -= src_;

		cNN.add(std::pair<grid_key_dx<dim>,grid_key_dx<dim>>(src,dst));

		++gr_sub3;

	}
};

/*! \brief Calculate the the Neighborhood for symmetric interactions
 *
 * \param cNN calculated cross neighborhood
 * \param div Number of divisions in each direction
 *
 */
template<unsigned int dim> void NNcalc_sym(openfpm::vector<grid_key_dx<dim>> & cNN)
{
	// Calculate the NNc_full array, it is a structure to get the neighborhood array

	// compile-time array {0,0,0,....}  {2,2,2,...} {1,1,1,...}

	typedef typename generate_array<size_t,dim, Fill_zero>::result NNzero;
	typedef typename generate_array<size_t,dim, Fill_two>::result NNtwo;

	// Generate the sub-grid iterator

	size_t div[dim];

	// Calculate the divisions

	for (size_t i = 0 ; i < dim ; i++)
		div[i] = 4;

	grid_sm<dim,void> gs(div);
	grid_key_dx_iterator_sub<dim> gr_sub3(gs,NNzero::data,NNtwo::data);

	grid_key_dx<dim> src_; // Source cell
	for (size_t i = 0; i < dim; i++)
		src_.set_d(i,1);

	size_t middle = gs.LinId(src_);

	// Calculate the symmetric array
	while (gr_sub3.isNext())
	{
		auto dst = gr_sub3.get();

		if ((long int)middle > gs.LinId(dst))
		{
			++gr_sub3;
			continue;
		}

		cNN.add(dst - src_);

		++gr_sub3;

	}
};

/*! \brief Calculate the Neighborhood cells
 *
 * \param cNN calculated cross neighborhood
 * \param div Number of divisions in each direction
 *
 */
template<unsigned int dim> void NNcalc_full(openfpm::vector<grid_key_dx<dim>> & cNN)
{
	// Calculate the NNc_full array, it is a structure to get the neighborhood array

	// compile-time array {0,0,0,....}  {2,2,2,...} {1,1,1,...}

	typedef typename generate_array<size_t,dim, Fill_zero>::result NNzero;
	typedef typename generate_array<size_t,dim, Fill_two>::result NNtwo;

	// Generate the sub-grid iterator

	size_t div[dim];

	// Calculate the divisions

	for (size_t i = 0 ; i < dim ; i++)
	{div[i] = 4;}

	grid_sm<dim,void> gs(div);
	grid_key_dx_iterator_sub<dim> gr_sub3(gs,NNzero::data,NNtwo::data);

	grid_key_dx<dim> src_; // Source cell
	for (size_t i = 0; i < dim; i++)
		src_.set_d(i,1);

	// Calculate the symmetric crs array
	while (gr_sub3.isNext())
	{
		auto dst = gr_sub3.get();

		cNN.add(dst - src_);

		++gr_sub3;
	}
};


/* NOTE all the implementations
 *
 * has complexity O(1) in getting the cell id and the elements in a cell
 * but with different constants
 *
 */

/*! \brief Class for FAST cell list implementation
 *
 * This class implement the FAST cell list, fast but memory
 * expensive. The memory allocation is (M * N_cell_max)*sizeof(ele) + M*8
 *
 * * M = number of cells
 * * N_cell_max = maximum number of elements in a cell
 * * ele = element the structure is storing
 *
 * \note Because N_cell_max >= N/M then M * N_cell_max >= O(N)
 *
 * \warning Not not use for high asymmetric distribution
 *
 * Example of a 2D Cell list 6x6 structure with padding 1 without shift, cell indicated with p are padding cell
 * the origin of the cell or point (0,0) is marked with cell number 9
 *
 * \verbatim
 * +-----------------------+
 * |p |p |p |p |p |p |p |p |
 * +-----------------------+
 * |p |  |  |  |  |  |  |p |
 * +-----------------------+
 * |p |  |  |  |  |  |  |p |
 * +-----------------------+
 * |p |  |  |  |  |  |  |p |
 * +-----------------------+
 * |p |9 |  |  |  |  |  |p |
 * +-----------------------+
 * |p |p |p |p |p |p |p |p |
 * +-----------------------+
 * \endverbatim
 *
 *
 * \tparam dim Dimensionality of the space
 * \tparam T type of the space float, double ...
 * \tparam base Base structure that store the information
 *
 * ### Declaration of a cell list
 * \snippet CellList_test.hpp Declare a cell list
 * ### Usage of cell list [CellS == CellList<3,double,FAST>]
 * \snippet CellList_test.hpp Usage of cell list
 * ### Remove one particle from each cell
 * \snippet CellList_test.hpp remove one particle from each cell
 * ### Usage of the neighborhood iterator
 * \snippet CellList_test.hpp Usage of the neighborhood iterator
 *
 */
template<unsigned int dim, typename T,  typename Mem_type, typename transform = no_transform<dim,T>, typename base=openfpm::vector<size_t>>
class CellList : public CellDecomposer_sm<dim,T,transform>, public Mem_type
{
protected:
	//! The array contain the neighborhood of the cell-id in case of asymmetric interaction
	//
	//    * * *
	//    * x *
	//    * * *


	NNc_array<dim,(unsigned int)openfpm::math::pow(3,dim)> NNc_full;
//	long int NNc_full[openfpm::math::pow(3,dim)];

	//! The array contain the neighborhood of the cell-id in case of symmetric interaction
	//
	//   * * *
	//     x *
	//
//	long int NNc_sym[openfpm::math::pow(3,dim)/2+1];
	NNc_array<dim,(unsigned int)openfpm::math::pow(3,dim)/2+1> NNc_sym;

private:

	//! Caching of r_cutoff radius
	wrap_unordered_map<T,openfpm::vector<long int>> rcache;

	//! True if has been initialized from CellDecomposer
	bool from_cd;

	//! Additional information in general (used to understand if the cell-list)
	//! has been constructed from an old decomposition
	size_t n_dec;

	/*! Calculate the neighborhood cells based on the radius
	 *
	 * \note To the calculated neighborhood cell you have to add the id of the central cell
	 *
		\verbatim
       +-----------------------+
       |p |p |p |p |p |p |p |p |
       +-----------------------+
       |p |  |  |  |  |  |  |p |
       +-----------------------+
       |p |  |  |7 |8 |9 |  |p |
       +-----------------------+
       |p |  |  |-1|0 |1 |  |p |
       +-----------------------+
       |p |9 |  |-9|-8|-7|  |p |
       +-----------------------+
       |p |p |p |p |p |p |p |p |
       +-----------------------+
		\endverbatim
	 *
	 * The number indicate the cell id calculated
	 *
	 * -9,-8,-7,-1,0,1,7,8,9
	 *
	 * The cell 0 has id = 22 in the big cell matrix, so to calculate the
	 * neighborhood cells you have to sum the id of the center cell
	 *
	 * 13,14,15,21,22,23,29,30,31
	 *
	 * \param r_cut Cutoff-radius
	 * \param NNcell vector containing the neighborhood cells ids
	 *
	 */
	void NNcalc(T r_cut, openfpm::vector<long int> & NNcell)
	{
		size_t n_cell[dim];
		size_t n_cell_mid[dim];

		Point<dim,T> spacing = this->getCellBox().getP2();
		const grid_sm<dim,void> & gs = this->getGrid();

		for (size_t i = 0 ; i < dim ; i++)
		{
			n_cell[i] = 2*(std::ceil(r_cut / spacing.get(i)))+1;
			n_cell_mid[i] = n_cell[i] / 2;
		}

		grid_sm<dim,void> gsc(n_cell);
		grid_key_dx_iterator<dim> gkdi(gsc);

		while (gkdi.isNext())
		{
			auto key = gkdi.get();

			for (size_t i = 0 ; i < dim ; i++)
				key.set_d(i,key.get(i) - n_cell_mid[i]);

			NNcell.add(gs.LinId(key));

			++gkdi;
		}
	}

	//! Initialize the structures of the data structure
	void InitializeStructures(const size_t (& div)[dim], size_t tot_n_cell, size_t slot=STARTING_NSLOT)
	{
		Mem_type::init_to_zero(slot,tot_n_cell);

		NNc_full.set_size(div);
		NNc_full.init_full();

		NNc_sym.set_size(div);
		NNc_sym.init_sym();
	}

	void setCellDecomposer(CellDecomposer_sm<dim,T,transform> & cd, const CellDecomposer_sm<dim,T,transform> & cd_sm, const Box<dim,T> & dom_box, size_t pad) const
	{
		size_t bc[dim];

		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = NON_PERIODIC;

		Box<dim,long int> bx = cd_sm.convertDomainSpaceIntoCellUnits(dom_box,bc);

		size_t div[dim];
		size_t div_big[dim];

		for (size_t i = 0 ; i < dim ; i++)
		{
			div[i] = bx.getHigh(i) - bx.getLow(i);
			div_big[i] = cd_sm.getDiv()[i] - 2*cd_sm.getPadding(i);
		}

		cd.setDimensions(cd_sm.getDomain(),div_big,div, pad, bx.getP1());
	}

public:

	//! Type of internal memory structure
	typedef Mem_type Mem_type_type;

	typedef CellNNIteratorSym<dim,CellList<dim,T,Mem_type,transform,base>,RUNTIME,NO_CHECK> SymNNIterator;

	//! Object type that the structure store
	typedef typename base::value_type value_type;

	//! Type of the coordinate space (double float)
	typedef T stype;

	/*! \brief Return the underlying grid information of the cell list
	 *
	 * \return the grid infos
	 *
	 */
	const grid_sm<dim,void> & getGrid()
	{
		return CellDecomposer_sm<dim,T,transform>::getGrid();
	}

	/*! Initialize the cell list from a well-define Cell-decomposer
	 *
	 * In some cases is needed to have a Cell-list with Cells consistent
	 * with a well predefined CellDecomposer. In this case we use this function.
	 * Using this initialization the Cell-list maintain the Cells defined by this
	 * Cell-decomposer consistently
	 *
	 * \param cd_sm Cell-Decomposer
	 * \param dom_box domain box (carefully this is going to be adjusted)
	 * \param pad cell-list padding
	 * \param slot slots for each cell
	 *
	 */
	void Initialize(CellDecomposer_sm<dim,T,transform> & cd_sm, const Box<dim,T> & dom_box,const size_t pad = 1, size_t slot=STARTING_NSLOT)
	{
		size_t bc[dim];
		for (size_t i = 0 ; i < dim ; i++)	{bc[i] = NON_PERIODIC;}

		Box<dim,long int> bx = cd_sm.convertDomainSpaceIntoCellUnits(dom_box,bc);

		setCellDecomposer(*this,cd_sm,dom_box,pad);

		size_t div_w_pad[dim];
		size_t tot_cell = 1;

		for (size_t i = 0 ; i < dim ; i++)
		{
			div_w_pad[i] = bx.getHigh(i) - bx.getLow(i) + 2*pad;
			tot_cell *= div_w_pad[i];
		}

		// here we set the cell-shift for the CellDecomposer

		InitializeStructures(div_w_pad,tot_cell);

		// Initialized from CellDecomposer
		from_cd = true;
	}

	/*! Initialize the cell list
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param pad padding cell
	 * \param slot maximum number of slot
	 *
	 */
	void Initialize(const Box<dim,T> & box, const size_t (&div)[dim], const size_t pad = 1, size_t slot=STARTING_NSLOT)
	{
		SpaceBox<dim,T> sbox(box);

		// Initialize point transformation

		Initialize(sbox,div,pad,slot);
	}

	/*! Initialize the cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param pad padding cell
	 * \param slot maximum number of slot
	 *
	 */
	void Initialize(const SpaceBox<dim,T> & box, const size_t (&div)[dim], const size_t pad = 1, size_t slot=STARTING_NSLOT)
	{
		Matrix<dim,T> mat;

		CellDecomposer_sm<dim,T,transform>::setDimensions(box,div, mat, pad);
		Mem_type::set_slot(slot);

		// create the array that store the number of particle on each cell and se it to 0
		InitializeStructures(this->gr_cell.getSize(),this->gr_cell.size());

		from_cd = false;
	}

	//! Default Constructor
	CellList()

	:Mem_type(STARTING_NSLOT)
	{};

	//! Copy constructor
	CellList(const CellList<dim,T,Mem_type,transform,base> & cell)
	:Mem_type(STARTING_NSLOT)
	{
		this->operator=(cell);
	}

	//! Copy constructor
	CellList(CellList<dim,T,Mem_type,transform,base> && cell)
	:Mem_type(STARTING_NSLOT)
	{
		this->operator=(cell);
	}


	/*! \brief Cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param mat Matrix transformation
	 * \param pad Cell padding
	 * \param slot maximum number of slot
	 *
	 */
	CellList(Box<dim,T> & box, const size_t (&div)[dim], Matrix<dim,T> mat, const size_t pad = 1, size_t slot=STARTING_NSLOT)
	:Mem_type(slot),CellDecomposer_sm<dim,T,transform>(box,div,mat,box.getP1(),pad)
	{
		SpaceBox<dim,T> sbox(box);
		Initialize(sbox,div,pad,slot);
	}

	/*! \brief Cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param pad Cell padding
	 * \param slot maximum number of slot
	 *
	 */
	CellList(Box<dim,T> & box, const size_t (&div)[dim], const size_t pad = 1, size_t slot=STARTING_NSLOT)
	:Mem_type(slot),n_dec(0)
	{
		SpaceBox<dim,T> sbox(box);
		Initialize(sbox,div,pad,slot);
	}

	/*! \brief Cell list constructor
	 *
	 * \param box Domain where this cell list is living
	 * \param div grid size on each dimension
	 * \param pad Cell padding
	 * \param slot maximum number of slot
	 *
	 */
	CellList(SpaceBox<dim,T> & box, const size_t (&div)[dim], const size_t pad = 1, size_t slot=STARTING_NSLOT)
	:Mem_type(slot)
	{
		Initialize(box,div,pad,slot);
	}

	/*! \brief Cell list constructor from a cell decomposer
	 *
	 * \see Initialize
	 *
	 * \param cd_sm Cell-Decomposer
	 * \param box domain box (carefully this is going to be adjusted)
	 * \param pad Cell list padding
	 * \param slot number of slot for each cell
	 *
	 */
	CellList(CellDecomposer_sm<dim,T,transform> & cd_sm, const Box<dim,T> & box, const size_t pad = 1, size_t slot=STARTING_NSLOT)
	:Mem_type(slot),n_dec(0)
	{
		Initialize(cd_sm,box,pad,slot);
	}

	/*! \brief Destructor
	 *
	 *
	 */
	~CellList()
	{}

	/*! \brief Constructor from a temporal object
	 *
	 * \param cell Cell list structure
	 *
	 * \return itself
	 *
	 */
	CellList<dim,T,Mem_type,transform,base> & operator=(CellList<dim,T,Mem_type,transform,base> && cell)
	{
		std::copy(&cell.NNc_full[0],&cell.NNc_full[openfpm::math::pow(3,dim)],&NNc_full[0]);
		std::copy(&cell.NNc_sym[0],&cell.NNc_sym[openfpm::math::pow(3,dim)/2+1],&NNc_sym[0]);

		Mem_type::swap(static_cast<Mem_type &&>(cell));

		static_cast<CellDecomposer_sm<dim,T,transform> &>(*this).swap(cell);

		n_dec = cell.n_dec;
		from_cd = cell.from_cd;

		return *this;
	}

	/*! \brief Constructor from a temporal object
	 *
	 * \param cell Cell list structure
	 *
	 * \return itself
	 *
	 */
	CellList<dim,T,Mem_type,transform,base> & operator=(const CellList<dim,T,Mem_type,transform,base> & cell)
	{
		NNc_full = cell.NNc_full;
		NNc_sym = cell.NNc_sym;

		Mem_type::operator=(static_cast<const Mem_type &>(cell));

		static_cast<CellDecomposer_sm<dim,T,transform> &>(*this) = static_cast<const CellDecomposer_sm<dim,T,transform> &>(cell);

		n_dec = cell.n_dec;
		from_cd = cell.from_cd;

		return *this;
	}

	/*! \brief Get an iterator over particles following the cell structure
	 *
	 * \param dom_cells cells in the domain
	 *
	 * \return a particle iterator
	 *
	 */
	ParticleIt_Cells<dim,CellList<dim,T,Mem_fast<>,transform,base>> getDomainIterator(openfpm::vector<size_t> & dom_cells)
	{
		ParticleIt_Cells<dim,CellList<dim,T,Mem_fast<>,transform,base>> it(*this,dom_cells);

		return it;
	}

	/*! \brief Add to the cell
	 *
	 * \param cell_id Cell id where to add
	 * \param ele element to add
	 *
	 */
	inline void addCell(size_t cell_id, typename base::value_type ele)
	{
		Mem_type::addCell(cell_id,ele);
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void add(const T (& pos)[dim], typename base::value_type ele)
	{
		// calculate the Cell id
		size_t cell_id = this->getCell(pos);

		Mem_type::add(cell_id,ele);
	}

	/*! \brief Add an element in the cell list
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void add(const Point<dim,T> & pos, typename base::value_type ele)
	{
		// calculate the Cell id
		size_t cell_id = this->getCell(pos);

		Mem_type::add(cell_id,ele);
	}


	/*! \brief Add an element in the cell list forcing to be in the domain cells
	 *
	 * \warning careful is intended to be used ONLY to avoid round-off problems
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void addDom(const T (& pos)[dim], typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = this->getCellDom(pos);

		// add the element to the cell

		addCell(cell_id,ele);
	}

	/*! \brief Add an element in the cell list forcing to be in the domain cells
	 *
	 * \warning careful is intended to be used ONLY to avoid round-off problems
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void addDom(const Point<dim,T> & pos, typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = this->getCellDom(pos);

		// add the element to the cell

		addCell(cell_id,ele);
	}

	/*! \brief Add an element in the cell list forcing to be in the padding cells
	 *
	 * \warning careful is intended to be used ONLY to avoid round-off problems
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void addPad(const T (& pos)[dim], typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = this->getCellPad(pos);

		// add the element to the cell

		addCell(cell_id,ele);
	}

	/*! \brief Add an element in the cell list forcing to be in the padding cells
	 *
	 * \warning careful is intended to be used ONLY to avoid round-off problems
	 *
	 * \param pos array that contain the coordinate
	 * \param ele element to store
	 *
	 */
	inline void addPad(const Point<dim,T> & pos, typename base::value_type ele)
	{
		// calculate the Cell id

		size_t cell_id = this->getCell(pos);

		// add the element to the cell

		addCell(cell_id,ele);
	}

	/*! \brief remove an element from the cell
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 */
	inline void remove(size_t cell, size_t ele)
	{
		Mem_type::remove(cell,ele);
	}

	/*! \brief Return the number of elements in the cell
	 *
	 * \param cell_id id of the cell
	 *
	 * \return number of elements in the cell
	 *
	 */
	inline size_t getNelements(const size_t cell_id) const
	{
		return Mem_type::getNelements(cell_id);
	}

	/*! \brief Get an element in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline auto get(size_t cell, size_t ele) -> decltype(this->Mem_type::get(cell,ele))
	{
		return Mem_type::get(cell,ele);
	}

	/*! \brief Get an element in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param cell cell id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline auto get(size_t cell, size_t ele) const -> decltype(this->Mem_type::get(cell,ele))
	{
		return Mem_type::get(cell,ele);
	}

	/*! \brief Swap the memory
	 *
	 * \param cl Cell list with witch you swap the memory
	 *
	 */
	inline void swap(CellList<dim,T,Mem_type,transform,base> & cl)
	{
		NNc_full.swap(cl.NNc_full);
		NNc_sym.swap(cl.NNc_sym);

		Mem_type::swap(static_cast<Mem_type &>(cl));

		static_cast<CellDecomposer_sm<dim,T,transform> &>(*this) = static_cast<const CellDecomposer_sm<dim,T,transform> &>(cl);
	}

	/*! \brief Get the Cell iterator
	 *
	 * \param cell
	 *
	 * \return the iterator to the elements inside cell
	 *
	 */
	CellIterator<CellList<dim,T,Mem_type,transform,base>> getCellIterator(size_t cell)
	{
		return CellIterator<CellList<dim,T,Mem_type,transform,base>>(cell,*this);
	}

	/*! \brief Get the Neighborhood iterator
	 *
	 * It iterate across all the element of the selected cell and the near cells
	 *
	 *  \verbatim

	     * * *
	     * x *
	     * * *

	   \endverbatim
	 *
	 * * x is the selected cell
	 * * * are the near cell
	 *
	 * \param cell cell id
	 *
	 * \return An iterator across the neighhood particles
	 *
	 */
	template<unsigned int impl=NO_CHECK> inline CellNNIterator<dim,CellList<dim,T,Mem_type,transform,base>,(int)FULL,impl> getNNIterator(size_t cell)
	{
		CellNNIterator<dim,CellList<dim,T,Mem_type,transform,base>,(int)FULL,impl> cln(cell,NNc_full,*this);
		return cln;

	}

	/*! \brief Get the symmetric Neighborhood iterator
	 *
	 * It iterate across all the element of the selected cell and the near cells up to some selected radius
	 *
	 * \param cell cell id
	 * \param r_cut radius
	 *
	 * \return An iterator across the neighborhood particles
	 *
	 */
	template<unsigned int impl=NO_CHECK> inline CellNNIteratorRadius<dim,CellList<dim,T,Mem_type,transform,base>,impl> getNNIteratorRadius(size_t cell, T r_cut)
	{
		openfpm::vector<long int> & NNc = rcache[r_cut];

		if (NNc.size() == 0)
			NNcalc(r_cut,NNc);

		CellNNIteratorRadius<dim,CellList<dim,T,Mem_type,transform,base>,impl> cln(cell,NNc,*this);

		return cln;
	}



	/*! \brief Get the symmetric Neighborhood iterator
	 *
	 * It iterate across all the element of the selected cell and the near cells
	 *
	 *  \verbatim

	   * * *
	     x *

	   \endverbatim
	 *
	 * * x is the selected cell
	 * * * are the near cell
	 *
	 * \param cell cell id
	 * \param p particle id
	 *
	 * \return An aiterator across the neighborhood particles
	 *
	 */
	template<unsigned int impl>
	inline CellNNIteratorSym<dim,CellList<dim,T,Mem_type,transform,base>,(unsigned int)SYM,impl>
	getNNIteratorSym(size_t cell, size_t p, const openfpm::vector<Point<dim,T>> & v)
	{
#ifdef SE_CLASS1
		if (from_cd == false)
		{std::cerr << __FILE__ << ":" << __LINE__ << " Warning when you try to get a symmetric neighborhood iterator, you must construct the Cell-list in a symmetric way" << std::endl;}
#endif

		CellNNIteratorSym<dim,CellList<dim,T,Mem_type,transform,base>,SYM,impl> cln(cell,p,NNc_sym,*this,v);
		return cln;
	}

	/*! \brief Get the symmetric Neighborhood iterator
	 *
	 * It iterate across all the element of the selected cell and the near cells
	 *
	 *  \verbatim

	   * * *
	     x *

	   \endverbatim
	 *
	 * * x is the selected cell
	 * * * are the near cell
	 *
	 * \param cell cell id
	 * \param p particle id
	 * \param v_p1 first phase for particle p
	 * \param v_p2 second phase for particle q
	 *
	 * \return An aiterator across the neighborhood particles
	 *
	 */
	template<unsigned int impl>
	inline CellNNIteratorSymMP<dim,CellList<dim,T,Mem_type,transform,base>,(unsigned int)SYM,impl>
	getNNIteratorSymMP(size_t cell, size_t p, const openfpm::vector<Point<dim,T>> & v_p1, const openfpm::vector<Point<dim,T>> & v_p2)
	{
#ifdef SE_CLASS1
		if (from_cd == false)
			std::cerr << __FILE__ << ":" << __LINE__ << " Warning when you try to get a symmetric neighborhood iterator, you must construct the Cell-list in a symmetric way" << std::endl;
#endif

		CellNNIteratorSymMP<dim,CellList<dim,T,Mem_type,transform,base>,SYM,impl> cln(cell,p,NNc_sym,*this,v_p1,v_p2);
		return cln;
	}

	/*! \brief Get the symmetric neighborhood
	 *
	 * \return the symmetric neighborhood
	 *
	 */
	const NNc_array<dim,(unsigned int)openfpm::math::pow(3,dim)/2+1> & getNNc_sym() const
	{
		return NNc_sym;
	}

	/*! \brief Return the number of padding cells of the Cell decomposer
	 *
	 * \param i dimension
	 *
	 * \return the number of padding cells
	 *
	 */
	size_t getPadding(size_t i) const
	{
		return CellDecomposer_sm<dim,T,transform>::getPadding(i);
	}


	/*! \brief Clear the cell list
	 *
	 */
	void clear()
	{
		Mem_type::clear();
	}

	/*! \brief Return the starting point of the cell p
	 *
	 * \param cell_id cell id
	 *
	 * \return the index
	 *
	 */
	inline const typename Mem_type::loc_index & getStartId(typename Mem_type::loc_index cell_id) const
	{
		return Mem_type::getStartId(cell_id);
	}

	/*! \brief Return the end point of the cell p
	 *
	 * \param cell_id cell id
	 *
	 * \return the stop index
	 *
	 */
	inline const typename Mem_type::loc_index & getStopId(typename Mem_type::loc_index cell_id) const
	{
		return Mem_type::getStopId(cell_id);
	}

	/*! \brief Return the neighborhood id
	 *
	 * \param part_id particle id
	 *
	 * \return the neighborhood id
	 *
	 */
	inline const typename Mem_type::loc_index & get_lin(const typename Mem_type::loc_index * part_id) const
	{
		return Mem_type::get_lin(part_id);
	}

//////////////////////////////// POINTLESS BUT REQUIRED TO RESPECT THE INTERFACE //////////////////

	//! Ghost marker
	size_t g_m = 0;

	/*! \brief return the ghost marker
	 *
	 * \return ghost marker
	 *
	 */
	inline size_t get_gm()
	{
		return g_m;
	}

	/*! \brief Set the ghost marker
	 *
	 * \param g_m marker
	 *
	 */
	inline void set_gm(size_t g_m)
	{
		this->g_m = g_m;
	}

/////////////////////////////////////

	/////////////////////////////////////

		/*! \brief Set the n_dec number
		 *
		 * \param n_dec
		 *
		 */
		void set_ndec(size_t n_dec)
		{
			this->n_dec = n_dec;
		}

		/*! \brief Set the n_dec number
		 *
		 * \return n_dec
		 *
		 */
		size_t get_ndec() const
		{
			return n_dec;
		}

	/////////////////////////////////////
};

/*! \brief Calculate parameters for the cell list
 *
 * \param div Division array
 * \param r_cut interation radius or size of each cell
 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost. This parameter says how much must be enlarged
 *
 * \return the processor bounding box
 */
template<unsigned int dim, typename St> static inline void cl_param_calculate(Box<dim, St> & pbox, size_t (&div)[dim], St r_cut, const Ghost<dim, St> & enlarge)
{
	// calculate the parameters of the cell list

	// extend by the ghost
	pbox.enlarge(enlarge);

	// Calculate the division array and the cell box
	for (size_t i = 0; i < dim; i++)
	{
		div[i] = static_cast<size_t>((pbox.getP2().get(i) - pbox.getP1().get(i)) / r_cut);
		div[i]++;
		pbox.setHigh(i,pbox.getLow(i) + div[i]*r_cut);
	}
}

/*! \brief Calculate parameters for the symmetric cell list
 *
 * \param[in] dom Simulation domain
 * \param[output] cd_sm This cell-decomposer is set according to the needed division
 * \param[in] g Ghost dimensions
 * \param[output] pad required padding for the cell-list
 *
 * \return the processor bounding box
 */
template<unsigned int dim, typename St> static inline void cl_param_calculateSym(const Box<dim,St> & dom, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, Ghost<dim,St> g, St r_cut, size_t & pad)
{
	size_t div[dim];

	for (size_t i = 0 ; i < dim ; i++)
		div[i] = (dom.getHigh(i) - dom.getLow(i)) / r_cut;

	g.magnify(1.013);

	// Calculate the maximum padding
	for (size_t i = 0 ; i < dim ; i++)
	{
		size_t tmp = std::ceil(fabs(g.getLow(i)) / r_cut);
		pad = (pad > tmp)?pad:tmp;

		tmp = std::ceil(fabs(g.getHigh(i)) / r_cut);
		pad = (pad > tmp)?pad:tmp;
	}

	cd_sm.setDimensions(dom,div,pad);
}


/*! \brief Calculate parameters for the symmetric cell list
 *
 * \param[in] dom Simulation domain
 * \param[out] cd_sm This cell-decomposer is set according to the needed division
 * \param[in] g Ghost part extension
 * \param[out] pad required padding for the cell-list
 *
 * \return the processor bounding box
 */
template<unsigned int dim, typename St> static inline void cl_param_calculateSym(const Box<dim,St> & dom, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, Ghost<dim,St> g, size_t & pad)
{
	size_t div[dim];

	for (size_t i = 0 ; i < dim ; i++)
		div[i] = (dom.getHigh(i) - dom.getLow(i)) / g.getHigh(i);

	g.magnify(1.013);

	pad = 1;

	cd_sm.setDimensions(dom,div,pad);
}


#endif /* CELLLIST_HPP_ */
