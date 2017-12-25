/*
 * VerletListM.hpp
 *
 *  Created on: Oct 15, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTM_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTM_HPP_

#include "NN/VerletList/VerletNNIteratorM.hpp"



/*! \brief Get the neighborhood iterator based on type
 *
 * \tparam dim dimensionality
 * \tparam T type of space
 * \tparam CellListImpl Cell-list implementation
 * \tparam type of neighborhood
 * \tparam PartIt Particle iterator
 *
 */
template<unsigned int dim, typename T, typename CellListImpl, typename PartIt, int type>
struct NNTypeM
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
						   const openfpm::vector<Point<dim,T>> & pos,
			               const openfpm::vector<pos_v<dim,T>> & v,
						   Point<dim,T> & xp,
						   size_t pp,
						   size_t p,
						   CellListImpl & cl,
						   T r_cut) -> decltype(cl.template getNNIterator<NO_CHECK>(0))
	{
		return cl.template getNNIterator<NO_CHECK>(cl.getCell(xp));
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
template<unsigned int dim, typename T, typename CellListImpl, typename PartIt>
struct NNTypeM<dim,T,CellListImpl,PartIt,VL_CRS_SYMMETRIC>
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
	static inline auto get(const PartIt & it,
						   const openfpm::vector<Point<dim,T>> & pos,
			               const openfpm::vector<pos_v<dim,T>> & v,
						   Point<dim,T> & xp,
						   size_t pp,
						   size_t p,
						   CellListImpl & cl,
						   T r_cut) -> decltype(it.getNNIteratorCSRM(pos,v))
	{
		return it.getNNIteratorCSRM(pos,v);
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
template<unsigned int dim, typename T, typename CellListImpl, typename PartIt>
struct NNTypeM<dim,T,CellListImpl,PartIt,VL_SYMMETRIC>
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
	static inline auto get(const PartIt & it,
			           const openfpm::vector<Point<dim,T>> & pos,
					   const openfpm::vector<pos_v<dim,T>> & v,
					   Point<dim,T> & xp,
					   size_t pp,
					   size_t p,
					   CellListImpl & cl,
					   T r_cut) -> decltype(cl.template getNNIteratorSym<NO_CHECK>(0,0,0,openfpm::vector<Point<dim,T>>(),openfpm::vector<pos_v<dim,T>>()))
	{
		return cl.template getNNIteratorSym<NO_CHECK>(cl.getCell(xp),pp,p,pos,v);
	}
};



/*! \brief Class for Verlet list implementation with Multiphase
 *
 * \tparam dim Dimensionality of the space
 * \tparam T type of the space float, double ...
 * \tparam CellListImpl Base structure that store the information
 *
 */
template<unsigned int dim,
		 typename T,
		 unsigned int sh_byte ,
		 typename CellListImpl=CellListM<dim,T,sh_byte>,
		 typename transform = shift<dim,T>,
		 typename VerletBase=VerletList<dim,T,Mem_fast,transform, size_t> >
class VerletListM : public VerletBase
{

	//! Mask to get the high bits of a number
	typedef boost::high_bit_mask_t<sh_byte>  mask_high;

	//! Mask to get the low bits of a number
	typedef boost::low_bits_mask_t<sizeof(size_t)*8-sh_byte>  mask_low;

	/*! \brief Create the Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param pos2 vector of neighborhood particles position
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cl Cell-list elements to use to construct the verlet list
	 * \param opt options to create the verlet list like VL_SYMMETRIC or VL_NON_SYMMETRIC
	 *
	 */
	inline void create(const openfpm::vector<Point<dim,T>> & pos,
			           const openfpm::vector<pos_v<dim,T>> & pos2 ,
					   const openfpm::vector<size_t> & dom,
					   const openfpm::vector<subsub_lin<dim>> & anom,
					   size_t pp,
					   T r_cut,
					   size_t g_m,
					   CellListImpl & cl,
					   size_t opt)
	{
		if (opt == VL_CRS_SYMMETRIC)
		{
			create_<CellNNIteratorSymM<dim,CellListImpl,sh_byte,RUNTIME,NO_CHECK>,VL_CRS_SYMMETRIC>(pos,pos2,dom,anom,pp,r_cut,g_m,cl,opt);
		}
		else if (opt == VL_SYMMETRIC)
		{
			create_<decltype(cl.template getNNIteratorSym<NO_CHECK>(0,0,0,pos,pos2)),VL_SYMMETRIC>(pos,pos2,dom,anom,pp,r_cut,g_m,cl,opt);
		}
		else
		{
			create_<decltype(cl.template getNNIterator<NO_CHECK>(0)),VL_NON_SYMMETRIC>(pos,pos2,dom,anom,pp,r_cut,g_m,cl,opt);
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
	template<typename NN_type, int type> inline void create_(const openfpm::vector<Point<dim,T>> & pos,
			                                                 const openfpm::vector<pos_v<dim,T>> & pos2 ,
															 const openfpm::vector<size_t> & dom,
															 const openfpm::vector<subsub_lin<dim>> & anom,
															 size_t pp,
															 T r_cut,
															 size_t g_m,
															 CellListImpl & cli,
															 size_t opt)
	{
		size_t end;

		auto it = PartItNN<type,dim,openfpm::vector<Point<dim,T>>,CellListImpl>::get(pos,dom,anom,cli,g_m,end);

		this->cl_n.resize(end);
		this->cl_base.resize(end*this->slot);
		this->cl_n.fill(0);

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		while (it.isNext())
		{
			size_t i = it.get();
			Point<dim,T> xp = pos.template get<0>(i);

			// Get the neighborhood of the particle
			NN_type NN = NNTypeM<dim,T,CellListImpl,decltype(it),type>::get(it,pos,pos2,xp,pp,i,cli,r_cut);

			while (NN.isNext())
			{
				size_t nnp = NN.getP();
				size_t v = NN.getV();
				size_t nnp_a = NN.get();

				Point<dim,T> xq = pos2.get(v).pos.template get<0>(nnp);

				if (xp.distance2(xq) < r_cut2)
					this->addPart(i,nnp_a);

				// Next particle
				++NN;
			}

			++it;
		}
	}

public:

	/*! Initialize the verlet list from an already filled cell-list
	 *
	 * \param cli external Cell-list
	 * \param pp phase of pos
	 * \param r_cut cutoff-radius
	 * \param pos vector of particle positions
	 * \param pos2 vector of particle position for the neighborhood
	 * \param g_m Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and g_m = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * 	\param opt options for the Verlet-list creation
	 *
	 */
	void Initialize(CellListImpl & cli, size_t pp, T r_cut, const openfpm::vector<Point<dim,T>> & pos, const openfpm::vector<struct pos_v<dim,T>> & pos2, size_t g_m, size_t opt = VL_NON_SYMMETRIC)
	{
		this->cl_n.resize(g_m);
		this->cl_base.resize(g_m*this->slot);

		Point<dim,T> spacing = cli.getCellBox().getP2();

		// Create with radius or not
		bool wr = true;

		for (size_t i = 0 ; i < dim ; i++)
			wr &= r_cut <= spacing.get(i);

		if (wr == true || opt == VL_SYMMETRIC)
		{
			openfpm::vector<subsub_lin<dim>> anom_c;
			openfpm::vector<size_t> dom_c;

			create(pos,pos2,dom_c,anom_c,pp,r_cut,g_m,cli,opt);
		}
		else
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " error iterator with radius is not implemented yet " << std::endl;
		}
	}

	/*! \brief Get the element-id in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param part part id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline size_t getP(size_t part, size_t ele) const
	{
		return VerletBase::get(part,ele) & mask_low::sig_bits_fast;
	}

	/*! \brief Get the element vector in the cell
	 *
	 * \tparam i property to get
	 *
	 * \param part particle id
	 * \param ele element id
	 *
	 * \return The element value
	 *
	 */
	inline size_t getV(size_t part, size_t ele) const
	{
		return (VerletBase::get(part,ele)) >> (sizeof(size_t)*8-sh_byte);
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
	template<unsigned int impl=NO_CHECK> inline VerletNNIteratorM<dim,VerletListM<dim,T,sh_byte,CellListImpl,transform,VerletBase>,sh_byte> getNNIterator(size_t part_id)
	{
		VerletNNIteratorM<dim,VerletListM<dim,T,sh_byte,CellListImpl,transform,VerletBase>,sh_byte> vln(part_id,*this);

		return vln;
	}
};


#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTM_HPP_ */
