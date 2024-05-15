/*
 * VerletListM.hpp
 *
 *  Created on: Oct 15, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTM_HPP_
#define OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTM_HPP_

#include "NN/VerletList/VerletNNIteratorM.hpp"
#include "VerletList.hpp"


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
		 typename vector_pos_type = openfpm::vector<Point<dim,T>>,
		 typename VerletBase=VerletList<dim,T,Mem_fast<>,transform, vector_pos_type> >
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
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cl Cell-list elements to use to construct the verlet list
	 * \param opt options to create the verlet list like VL_SYMMETRIC or VL_NON_SYMMETRIC
	 *
	 */
	inline void create(
		const vector_pos_type & pos,
		const openfpm::vector<pos_v<vector_pos_type>> & pos2 ,
		size_t pp,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cl,
		size_t opt)
	{
		if (opt == VL_CRS_SYMMETRIC)
		{
			openfpm::vector<subsub_lin<dim>> anom_c;
			openfpm::vector<size_t> dom_c;
			createCRSSymmetric(pos,pos2,dom_c,anom_c,pp,r_cut,ghostMarker,cl,opt);
		}
		else if (opt == VL_SYMMETRIC)
		{
			createSymmetric(pos,pos2,pp,r_cut,ghostMarker,cl,opt);
		}
		else
		{
			createNonSymmetric(pos,pos2,pp,r_cut,ghostMarker,cl,opt);
		}
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
	 * \param opt options
	 *
	 */
	inline void createCRSSymmetric(
		const vector_pos_type & pos,
		const openfpm::vector<pos_v<vector_pos_type>> & pos2 ,
		const openfpm::vector<size_t> & dom,
		const openfpm::vector<subsub_lin<dim>> & anom,
		size_t pp,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cli,
		size_t opt)
	{
		size_t end = pos.size();

		this->init_to_zero(this->slot,end);

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		auto it = ParticleItCRS_Cells<dim,CellListImpl,vector_pos_type>(cli,dom,anom,cli.getNNc_sym());
		// iterate the particles
		while (it.isNext())
		{
			size_t i = it.get();
			Point<dim,T> xp = pos.template get<0>(i);

			// Get the neighborhood of the particle
			auto NN = it.getNNIteratorCSRM(pos,pos2);

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

	/*! \brief Create the Symmetric Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param pos2 vector of position for the neighborhood
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cli Cell-list elements to use to construct the verlet list
	 * \param opt options
	 *
	 */
	inline void createSymmetric(
		const vector_pos_type & pos,
		const openfpm::vector<pos_v<vector_pos_type>> & pos2 ,
		size_t pp,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cli,
		size_t opt)
	{
		size_t end = ghostMarker;
		this->init_to_zero(this->slot,end);

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		auto it = pos.getIteratorTo(end);
		while (it.isNext())
		{
			size_t i = it.get();
			Point<dim,T> xp = pos.template get<0>(i);

			// Get the neighborhood of the particle
			auto NN = cli.getNNIteratorSym(cli.getCell(xp),pp,i,pos,pos2);

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

	/*! \brief Create the Non-symmetric Verlet list from a given cell-list
	 *
	 * \param pos vector of positions
	 * \param pos2 vector of position for the neighborhood
	 * \param r_cut cut-off radius to get the neighborhood particles
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * \param cli Cell-list elements to use to construct the verlet list
	 * \param opt options
	 *
	 */
	inline void createNonSymmetric(
		const vector_pos_type & pos,
		const openfpm::vector<pos_v<vector_pos_type>> & pos2 ,
		size_t pp,
		T r_cut,
		size_t ghostMarker,
		CellListImpl & cli,
		size_t opt)
	{
		size_t end = ghostMarker;

		this->init_to_zero(this->slot,end);

		// square of the cutting radius
		T r_cut2 = r_cut * r_cut;

		// iterate the particles
		auto it = pos.getIteratorTo(end);
		while (it.isNext())
		{
			size_t i = it.get();
			Point<dim,T> xp = pos.template get<0>(i);

			// Get the neighborhood of the particle
			auto NN = cli.getNNIterator(cli.getCell(xp));

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
	 * \param ghostMarker Indicate form which particles to construct the verlet list. For example
	 * 			if we have 120 particles and ghostMarker = 100, the Verlet list will be constructed only for the first
	 * 			100 particles
	 * 	\param opt options for the Verlet-list creation
	 *
	 */
	void Initialize(CellListImpl & cli, size_t pp,
		T r_cut, const vector_pos_type & pos,
		const openfpm::vector<struct pos_v<vector_pos_type>> & pos2, size_t ghostMarker,
		size_t opt = VL_NON_SYMMETRIC)
	{
//		this->cl_n.resize(ghostMarker);
//		this->cl_base.resize(ghostMarker*this->slot);

		Point<dim,T> spacing = cli.getCellBox().getP2();

		// Create with radius or not
		bool wr = true;

		for (size_t i = 0 ; i < dim ; i++)
		{wr &= r_cut <= spacing.get(i);}

		if (wr == true || opt == VL_SYMMETRIC)
		{
			create(pos,pos2,pp,r_cut,ghostMarker,cli,opt);
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
	inline VerletNNIteratorM<dim,VerletListM<dim,T,sh_byte,CellListImpl,transform,vector_pos_type,VerletBase>,sh_byte> getNNIterator(size_t part_id)
	{
		VerletNNIteratorM<dim,VerletListM<dim,T,sh_byte,CellListImpl,transform,vector_pos_type,VerletBase>,sh_byte> vln(part_id,*this);

		return vln;
	}
};


#endif /* OPENFPM_DATA_SRC_NN_VERLETLIST_VERLETLISTM_HPP_ */
