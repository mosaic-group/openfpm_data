/*
 * CellList_util.hpp
 *
 *  Created on: Oct 20, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_
#define OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_

#define CL_SYMMETRIC 1
#define CL_NON_SYMMETRIC 2

enum cl_construct_opt
{
	Full,
	Only_reorder
};

#include "util/ofp_context.hpp"


/*! \brief populate the Cell-list with particles non symmetric case on GPU
 *
 * \tparam dim dimensionality of the space
 * \tparam T type of the space
 * \tparam CellList type of cell-list
 *
 * \param pos vector of positions
 * \param cli Cell-list
 * \param g_m marker (particle below this marker must be inside the domain, particles outside this marker must be outside the domain)
 *
 */
template<unsigned int dim, typename T, typename CellList> void populate_cell_list_no_sym_gpu(openfpm::vector<Point<dim,T>> & pos, CellList & cli, size_t g_m)
{
	// First we have to set the counter to zero
	cli.PopulateOnGPU(pos);
}

#include "Vector/map_vector.hpp"

template<bool is_gpu>
struct populate_cell_list_no_sym_impl
{
	template<unsigned int dim, typename T, typename prop, typename Memory, template <typename> class layout_base , typename CellList,unsigned int ... prp>
	static void populate(openfpm::vector<Point<dim,T>,Memory,layout_base > & pos,
						   openfpm::vector<Point<dim,T>,Memory,layout_base > & v_pos_out,
						   openfpm::vector<prop,Memory,layout_base > & v_prp,
						   openfpm::vector<prop,Memory,layout_base > & v_prp_out,
			   	   	   	   CellList & cli,
						   gpu::ofp_context_t & context,
			   	   	   	   size_t g_m,
			   	   	   	   cl_construct_opt optc)
	{
		cli.clear();

		for (size_t i = 0; i < pos.size() ; i++)
		{
			cli.add(pos.get(i), i);
		}
	}
};

template<>
struct populate_cell_list_no_sym_impl<true>
{
	template<unsigned int dim, typename T, typename prop, typename Memory, template <typename> class layout_base , typename CellList, unsigned int ... prp>
	static void populate(openfpm::vector<Point<dim,T>,Memory,layout_base > & pos,
						 openfpm::vector<Point<dim,T>,Memory,layout_base > & v_pos_out,
						 openfpm::vector<prop,Memory,layout_base > & v_prp,
						 openfpm::vector<prop,Memory,layout_base > & v_prp_out,
			   	   	   	   CellList & cli,
						   gpu::ofp_context_t & context,
			   	   	   	   size_t g_m,
			   	   	   	   cl_construct_opt optc)
	{
		v_prp_out.resize(pos.size());
		v_pos_out.resize(pos.size());

		cli.template construct<decltype(pos),decltype(v_prp),prp ...>(pos,v_pos_out,v_prp,v_prp_out,context,g_m,0,pos.size(),optc);
	}
};

template<bool is_gpu>
struct populate_cell_list_sym_impl
{
	template<unsigned int dim, typename T, typename Memory, template <typename> class layout_base , typename CellList>
	static void populate(openfpm::vector<Point<dim,T>,Memory,layout_base > & pos,
			   	   	   	   CellList & cli,
			   	   	   	   size_t g_m)
	{
		cli.clear();

		for (size_t i = 0; i < g_m ; i++)
		{
			cli.addDom(pos.get(i), i);
		}

		for (size_t i = g_m; i < pos.size() ; i++)
		{
			cli.addPad(pos.get(i), i);
		}
	}
};

template<>
struct populate_cell_list_sym_impl<true>
{
	template<unsigned int dim, typename T, typename Memory, template <typename> class layout_base , typename CellList>
	static void populate(openfpm::vector<Point<dim,T>,Memory,layout_base > & pos,
			   	   	   	   CellList & cli,
			   	   	   	   size_t g_m)
	{
		std::cout << __FILE__ << ":" << __LINE__ << " symmetric cell list on GPU is not implemented. (And will never be, race conditions make them non suitable for GPU)" << std::endl;
	}
};

/*! \brief populate the Cell-list with particles non symmetric case
 *
 * \tparam dim dimensionality of the space
 * \tparam T type of the space
 * \tparam CellList type of cell-list
 *
 * \param pos vector of positions
 * \param cli Cell-list
 * \param g_m marker (particle below this marker must be inside the domain, particles outside this marker must be outside the domain)
 *
 */
template<unsigned int dim,
		 typename T,
		 typename prop,
		 typename Memory,
		 template <typename> class layout_base ,
		 typename CellList,
		 unsigned int ... prp>
void populate_cell_list_no_sym(openfpm::vector<Point<dim,T>,Memory,layout_base > & pos,
		 	 	 	 	 	   openfpm::vector<Point<dim,T>,Memory,layout_base > & v_pos_out,
		 	 	 	 	 	   openfpm::vector<prop,Memory,layout_base > & v_prp,
		 	 	 	 	 	   openfpm::vector<prop,Memory,layout_base > & v_prp_out,
							   CellList & cli,
							   gpu::ofp_context_t & gpu,
							   size_t g_m,
							   cl_construct_opt optc)
{
	populate_cell_list_no_sym_impl<is_gpu_celllist<CellList>::value>
								  ::template populate<dim,T,prop,Memory,layout_base,CellList, prp ...>(pos,v_pos_out,v_prp,v_prp_out,cli,gpu,g_m,optc);
}

/*! \brief populate the Cell-list with particles symmetric case
 *
 * \tparam dim dimensionality of the space
 * \tparam T type of the space
 * \tparam CellList type of cell-list
 *
 * \param pos vector of positions
 * \param cli Cell-list
 * \param g_m marker (particle below this marker must be inside the domain, particles outside this marker must be outside the domain)
 *
 */
template<unsigned int dim, typename T, typename Memory, template <typename> class layout_base ,typename CellList>
void populate_cell_list_sym(openfpm::vector<Point<dim,T>,Memory,layout_base > & pos,
		      	  	  	    CellList & cli,
		      	  	  	    size_t g_m)
{
	populate_cell_list_sym_impl<is_gpu_celllist<CellList>::value>::populate(pos,cli,g_m);
}

/*! \brief populate the Cell-list with particles generic case
 *
 * \tparam dim dimensionality of the space
 * \tparam T type of the space
 * \tparam CellList type of cell-list
 *
 * \param pos vector of positions
 * \param cli Cell-list
 * \param opt option like CL_SYMMETRIC or CL_NON_SYMMETRIC
 * \param g_m marker (particle below this marker must be inside the domain, particles outside this marker must be outside the domain)
 *
 */
template<unsigned int dim,
		 typename T,
		 typename prop,
		 typename Memory,
		 template <typename> class layout_base,
		 typename CellList,
		 unsigned int ... prp>
void populate_cell_list(openfpm::vector<Point<dim,T>,Memory,layout_base> & pos,
 	 	   	   	   	    openfpm::vector<Point<dim,T>,Memory,layout_base > & v_pos_out,
 	 	   	   	   	    openfpm::vector<prop,Memory,layout_base > & v_prp,
 	 	   	   	   	    openfpm::vector<prop,Memory,layout_base > & v_prp_out,
						CellList & cli,
						gpu::ofp_context_t & context,
						size_t g_m,
						size_t opt,
						cl_construct_opt optc)
{
	if (opt == CL_NON_SYMMETRIC)
	{populate_cell_list_no_sym<dim,T,prop,Memory,layout_base,CellList,prp ...>(pos,v_pos_out,v_prp,v_prp_out,cli,context,g_m,optc);}
	else
	{populate_cell_list_sym(pos,cli,g_m);}
}

/*! \brief populate the Cell-list with particles generic case
 *
 * \note this function remain for backward compatibility it supposed to be remove once the verlet-list use the new populate-cell list form
 *
 * \tparam dim dimensionality of the space
 * \tparam T type of the space
 * \tparam CellList type of cell-list
 *
 * \param pos vector of positions
 * \param cli Cell-list
 * \param opt option like CL_SYMMETRIC or CL_NON_SYMMETRIC
 * \param g_m marker (particle below this marker must be inside the domain, particles outside this marker must be outside the domain)
 *
 */
template<unsigned int dim,
		 typename T,
		 typename Memory,
		 template <typename> class layout_base,
		 typename CellList,
		 unsigned int ... prp>
void populate_cell_list(openfpm::vector<Point<dim,T>,Memory,layout_base> & pos,
						CellList & cli,
						gpu::ofp_context_t & context,
						size_t g_m,
						size_t opt,
						cl_construct_opt optc)
{
	typedef openfpm::vector<aggregate<int>,Memory,layout_base> stub_prop_type;

	stub_prop_type stub1;
	stub_prop_type stub2;

	openfpm::vector<Point<dim,T>,Memory,layout_base> stub3;

	populate_cell_list<dim,T,aggregate<int>,Memory,layout_base,CellList,prp...>(pos,stub3,stub1,stub2,cli,context,g_m,opt,optc);
}

/*! \brief Structure that contain a reference to a vector of particles
 *
 *
 */
template<typename vector_pos_type>
struct pos_v
{
	vector_pos_type & pos;

	pos_v(vector_pos_type & pos)
	:pos(pos)
	{}
};

#endif /* OPENFPM_DATA_SRC_NN_CELLLIST_CELLLIST_UTIL_HPP_ */
