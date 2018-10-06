/*
 * CellList_gpu_ker.cuh
 *
 *  Created on: Jul 30, 2018
 *      Author: i-bird
 */

#ifndef CELLLIST_GPU_KER_CUH_
#define CELLLIST_GPU_KER_CUH_

template<unsigned int dim, typename cnt_type, typename ids_type>
class NN_gpu_it
{
	grid_key_dx<dim,ids_type> cell_act;

	grid_key_dx<dim,ids_type> cell_start;
	grid_key_dx<dim,ids_type> cell_stop;

	const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & starts;

	const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt;

	const openfpm::array<ids_type,dim,cnt_type> & div_c;

	const openfpm::array<ids_type,dim,cnt_type> & off;

	cnt_type p_id;
	cnt_type c_id;

	__device__ void SelectValid()
	{
		while (p_id >= starts.template get<0>(c_id+1) && isNext())
		{
			cnt_type id = cell_act.get(0);
			cell_act.set_d(0,id+1);

			//! check the overflow of all the index with exception of the last dimensionality

			int i = 0;
			for ( ; i < dim-1 ; i++)
			{
				size_t id = cell_act.get(i);
				if ((int)id > cell_stop.get(i))
				{
					// ! overflow, increment the next index

					cell_act.set_d(i,cell_start.get(i));
					id = cell_act.get(i+1);
					cell_act.set_d(i+1,id+1);
				}
				else
				{
					break;
				}
			}

			c_id = cid_<dim,cnt_type,ids_type,int>::get_cid(div_c,cell_act);
			p_id = starts.template get<0>(c_id);
		}
	}

	__device__ void SelectValid_debug()
	{
		while (p_id >= starts.template get<0>(c_id+1) && isNext())
		{
			cnt_type id = cell_act.get(0);
			cell_act.set_d(0,id+1);

			//! check the overflow of all the index with exception of the last dimensionality

			int i = 0;
			for ( ; i < dim-1 ; i++)
			{
				size_t id = cell_act.get(i);
				if ((int)id > cell_stop.get(i))
				{
					// ! overflow, increment the next index

					cell_act.set_d(i,cell_start.get(i));
					id = cell_act.get(i+1);
					cell_act.set_d(i+1,id+1);
				}
				else
				{
					break;
				}
			}

			c_id = cid_<dim,cnt_type,ids_type,int>::get_cid(div_c,cell_act);
			p_id = starts.template get<0>(c_id);
		}
	}

public:

	__device__ NN_gpu_it(const grid_key_dx<dim,ids_type> & cell_pos,
			             const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & starts,
			             const openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & srt,
			             const openfpm::array<ids_type,dim,cnt_type> & div_c,
			             const openfpm::array<ids_type,dim,cnt_type> & off)
	:starts(starts),srt(srt),div_c(div_c),off(off)
	{
		// calculate start and stop

		for (size_t i = 0 ; i < dim ; i++)
		{
			cell_start.set_d(i,cell_pos.get(i) - 1);
			cell_stop.set_d(i,cell_pos.get(i) + 1);
			cell_act.set_d(i,cell_pos.get(i) - 1);
		}

		c_id = cid_<dim,cnt_type,ids_type,int>::get_cid(div_c,cell_start);
		p_id = starts.template get<0>(c_id);

		SelectValid();
	}

	__device__ cnt_type get_sort()
	{
		return p_id;
	}

	__device__ cnt_type get()
	{
		return srt.template get<0>(p_id);
	}

	__device__ NN_gpu_it<dim,cnt_type,ids_type> & operator++()
	{
		++p_id;

		SelectValid();

		return *this;
	}

	__device__ NN_gpu_it<dim,cnt_type,ids_type> & plusplus()
	{
		++p_id;

		SelectValid_debug();

		return *this;
	}

	__device__ cnt_type get_start(unsigned int ce_id)
	{
		return starts.template get<0>(ce_id);
	}

	__device__ cnt_type get_cid()
	{
		return c_id;
	}

	__device__ bool isNext()
	{
		return cell_act.get(dim-1) <= cell_stop.get(dim-1);
	}
};

template<unsigned int dim, typename T, typename cnt_type, typename ids_type, typename transform>
class CellList_gpu_ker
{
	//! starting point for each cell
	openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> starts;

	//! Sorted to non sorted ids conversion
	openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> srt;

	//! Domain particles ids
	openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> dprt;

	//! Spacing
	openfpm::array<T,dim,cnt_type> spacing_c;

	//! \brief number of sub-divisions in each direction
	openfpm::array<ids_type,dim,cnt_type> div_c;

	//! \brief cell offset
	openfpm::array<ids_type,dim,cnt_type> off;

	//! transformation
	transform t;

public:

	CellList_gpu_ker(openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> starts,
					 openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> srt,
					 openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> dprt,
					 openfpm::array<T,dim,cnt_type> & spacing_c,
			         openfpm::array<ids_type,dim,cnt_type> & div_c,
			         openfpm::array<ids_type,dim,cnt_type> & off,
			         const transform & t)
	:starts(starts),srt(srt),dprt(dprt),spacing_c(spacing_c),div_c(div_c),off(off),t(t)
	{}

	inline __device__ grid_key_dx<dim,ids_type> getCell(const Point<dim,T> & xp) const
	{
		return cid_<dim,cnt_type,ids_type,transform>::get_cid_key(spacing_c,off,t,xp);
	}

	inline __device__ NN_gpu_it<dim,cnt_type,ids_type> getNNIterator(const grid_key_dx<dim,ids_type> & cid)
	{
		NN_gpu_it<dim,cnt_type,ids_type> ngi(cid,starts,srt,div_c,off);

		return ngi;
	}

	inline __device__ openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & getDomainSortIds()
	{
		return dprt;
	}

	inline __device__ openfpm::vector_gpu_ker<aggregate<cnt_type>,memory_traits_inte> & getSortToNonSort()
	{
		return srt;
	}
};


#endif /* CELLLIST_GPU_KER_CUH_ */
