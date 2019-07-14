/*
 * copy_grid_fast.hpp
 *
 *  Created on: Nov 29, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_GRID_COPY_GRID_FAST_HPP_
#define OPENFPM_DATA_SRC_GRID_COPY_GRID_FAST_HPP_

#include "Grid/iterators/grid_key_dx_iterator.hpp"

template<unsigned int dim>
struct striding
{
	size_t striding_src[dim - 1];
	size_t striding_dst[dim - 1];
	size_t n_cpy;
	size_t tot_y;
};

//////////////////////////////////// Functor to copy 1D grid in device memory ////////////

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one device memory into another device memory
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<bool lin_or_inte,typename data_type,typename S>
struct copy_fast_1d_device_memory
{
	//! set of pointers
	data_type & data_src;

	data_type & data_dst;

	/*! \brief constructor
	 *
	 * \param v set of pointer buffers to set
	 *
	 */
	inline copy_fast_1d_device_memory(data_type & data_src, data_type & data_dst)
	:data_src(data_src),data_dst(data_dst)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		static_cast<S *>(boost::fusion::at_c<T::value>(data_dst).mem)->copyDeviceToDevice(*static_cast<S *>(boost::fusion::at_c<T::value>(data_src).mem));
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one device memory into another device memory
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<typename data_type,typename S>
struct copy_fast_1d_device_memory<true,data_type,S>
{
	//! set of pointers
	data_type & data_src;

	data_type & data_dst;

	/*! \brief constructor
	 *
	 * \param v set of pointer buffers to set
	 *
	 */
	inline copy_fast_1d_device_memory(data_type & data_src, data_type & data_dst)
	:data_src(data_src),data_dst(data_dst)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		if (T::value == 0)
		{
			static_cast<S *>(data_dst.mem)->copyDeviceToDevice(*static_cast<S *>(data_src.mem));
		}
	}
};

///////////////////////////////////////////////////////////////////////////////////////////

/*! \brief This is a way to quickly copy a grid into another grid
 *
 *
 */
template<bool is_complex, unsigned int N, typename grid, typename ginfo>
struct copy_grid_fast
{
	static void copy(ginfo & gs_src,
				   ginfo & gs_dst,
				   const Box<N,size_t> & bx_src,
				   const Box<N,size_t> & bx_dst,
				   const grid & gd_src,
				   grid & gd_dst,
				   grid_key_dx<N> (& cnt)[1] )
	{
		grid_key_dx_iterator_sub<N,stencil_offset_compute<N,1>> sub_src(gs_src,bx_src.getKP1(),bx_src.getKP2(),cnt);
		grid_key_dx_iterator_sub<N,stencil_offset_compute<N,1>> sub_dst(gs_dst,bx_dst.getKP1(),bx_dst.getKP2(),cnt);

		while (sub_src.isNext())
		{
			// Option 1
			gd_dst.set(sub_dst.template getStencil<0>(),gd_src,sub_src.template getStencil<0>());

			++sub_src;
			++sub_dst;
		}
	}
};

/*! \brief This is a way to quickly copy a grid into another grid
 *
 *
 */
template<typename grid, typename ginfo>
struct copy_grid_fast<true,3,grid,ginfo>
{
	static void copy(ginfo & gs_src,
				   ginfo & gs_dst,
				   const Box<3,size_t> & bx_src,
				   const Box<3,size_t> & bx_dst,
				   const grid & gd_src,
				   grid & gd_dst,
				   grid_key_dx<3> (& cnt)[1] )
	{
		size_t lin_src = 0;
		size_t lin_dst = 0;

		lin_src += bx_src.getLow(2) * gs_src.size_s(1);
		lin_dst += bx_dst.getLow(2) * gs_dst.size_s(1);
		for (size_t i = bx_src.getLow(2) ; i <= bx_src.getHigh(2) ; i++)
		{
			lin_src += bx_src.getLow(1) * gs_src.size_s(0);
			lin_dst += bx_dst.getLow(1) * gs_dst.size_s(0);
			for (size_t j = bx_src.getLow(1) ; j <= bx_src.getHigh(1) ; j++)
			{
				lin_src += bx_src.getLow(0);
				lin_dst += bx_dst.getLow(0);
				for (size_t k = bx_src.getLow(0) ; k <= bx_src.getHigh(0) ; k++)
				{
					gd_dst.set(lin_dst,gd_src,lin_src);

					lin_src++;
					lin_dst++;
				}
				lin_src -= bx_src.getHigh(0) + 1;
				lin_dst -= bx_dst.getHigh(0) + 1;
				lin_src += gs_src.size_s(0);
				lin_dst += gs_dst.size_s(0);
			}
			lin_src -= (bx_src.getHigh(1) + 1)*gs_src.size_s(0);
			lin_dst -= (bx_dst.getHigh(1) + 1)*gs_dst.size_s(0);
			lin_src += gs_src.size_s(1);
			lin_dst += gs_dst.size_s(1);
		}
	}
};


/*! \brief This is a way to quickly copy a grid into another grid
 *
 *
 */
template<typename grid, typename ginfo>
struct copy_grid_fast<true,2,grid,ginfo>
{
	static void copy(ginfo & gs_src,
				   ginfo & gs_dst,
				   const Box<2,size_t> & bx_src,
				   const Box<2,size_t> & bx_dst,
				   const grid & gd_src,
				   grid & gd_dst,
				   grid_key_dx<2> (& cnt)[1] )
	{
		size_t lin_src = 0;
		size_t lin_dst = 0;


		lin_src += bx_src.getLow(1) * gs_src.size_s(0);
		lin_dst += bx_dst.getLow(1) * gs_dst.size_s(0);
		for (size_t j = bx_src.getLow(1) ; j <= bx_src.getHigh(1) ; j++)
		{
			lin_src += bx_src.getLow(0);
			lin_dst += bx_dst.getLow(0);
			for (size_t k = bx_src.getLow(0) ; k <= bx_src.getHigh(0) ; k++)
			{
				gd_dst.set(lin_dst,gd_src,lin_src);

				lin_src++;
				lin_dst++;
			}
			lin_src -= bx_src.getHigh(0) + 1;
			lin_dst -= bx_dst.getHigh(0) + 1;
			lin_src += gs_src.size_s(0);
			lin_dst += gs_dst.size_s(0);
		}

	}
};



/*! \brief This is a way to quickly copy a grid into another grid
 *
 *
 */
template<typename grid, typename ginfo>
struct copy_grid_fast<true,1,grid,ginfo>
{
	static void copy(ginfo & gs_src,
				   ginfo & gs_dst,
				   const Box<1,size_t> & bx_src,
				   const Box<1,size_t> & bx_dst,
				   const grid & gd_src,
				   grid & gd_dst,
				   grid_key_dx<1> (& cnt)[1] )
	{
		size_t lin_src = 0;
		size_t lin_dst = 0;

		lin_src += bx_src.getLow(0);
		lin_dst += bx_dst.getLow(0);
		for (size_t k = bx_src.getLow(0) ; k <= bx_src.getHigh(0) ; k++)
		{
			gd_dst.set(lin_dst,gd_src,lin_src);

			lin_src++;
			lin_dst++;
		}


	}
};

////// In case the property is not complex

template<unsigned int object_size>
void copy_grid_fast_longx_3(const Box<3,size_t> & bx_src,
							unsigned char * ptr_dst,
							unsigned char * ptr_src,
							striding<3> & sr)
{
	for (size_t i = bx_src.getLow(2) ; i <= bx_src.getHigh(2) ; i++)
	{
		for (size_t j = bx_src.getLow(1) ; j <= bx_src.getHigh(1) ; j++)
		{
			memcpy(ptr_dst,ptr_src,sr.n_cpy*object_size);

			ptr_dst += sr.striding_dst[0];
			ptr_src += sr.striding_src[0];
		}
		ptr_dst += sr.striding_dst[1] - sr.tot_y*sr.striding_dst[0];
		ptr_src += sr.striding_src[1] - sr.tot_y*sr.striding_src[0];
	}
}

template<unsigned int object_size, unsigned int n_cpy>
void copy_grid_fast_shortx_3(const Box<3,size_t> & bx_src,
		unsigned char * ptr_dst,
		unsigned char * ptr_src,
		striding<3> & sr)
{
	for (size_t i = bx_src.getLow(2) ; i <= bx_src.getHigh(2) ; i++)
	{
		for (size_t j = bx_src.getLow(1) ; j <= bx_src.getHigh(1) ; j++)
		{
			__builtin_memcpy(ptr_dst,ptr_src,n_cpy*object_size);

			ptr_dst += sr.striding_dst[0];
			ptr_src += sr.striding_src[0];
		}
		ptr_dst += sr.striding_dst[1] - sr.tot_y*sr.striding_dst[0];
		ptr_src += sr.striding_src[1] - sr.tot_y*sr.striding_src[0];
	}
}


////// In case the property is not complex

template<unsigned int object_size>
void copy_grid_fast_longx_2(const Box<2,size_t> & bx_src,
							unsigned char * ptr_dst,
							unsigned char * ptr_src,
							striding<2> & sr)
{
	for (size_t j = bx_src.getLow(1) ; j <= bx_src.getHigh(1) ; j++)
	{
		memcpy(ptr_dst,ptr_src,sr.n_cpy*object_size);

		ptr_dst += sr.striding_dst[0];
		ptr_src += sr.striding_src[0];
	}
}

template<unsigned int object_size, unsigned int n_cpy>
void copy_grid_fast_shortx_2(const Box<2,size_t> & bx_src,
		unsigned char * ptr_dst,
		unsigned char * ptr_src,
		striding<2> & sr)
{
	for (size_t j = bx_src.getLow(1) ; j <= bx_src.getHigh(1) ; j++)
	{
		__builtin_memcpy(ptr_dst,ptr_src,n_cpy*object_size);

		ptr_dst += sr.striding_dst[0];
		ptr_src += sr.striding_src[0];
	}
}

template<unsigned int dim>
struct copy_ndim_fast_selector
{
	template<unsigned int object_size>
	static void call(unsigned char * ptr_src,
			 unsigned char * ptr_dst,
			 striding<dim> & sr,
  const Box<dim,size_t> & bx_src)
	{
		std::cout << __FILE__ << ":" << __LINE__ << " error, wrong specialization. this case is not handled" << std::endl;
	}
};

template<>
struct copy_ndim_fast_selector<3>
{
	template<unsigned int object_size>
	static void call(unsigned char * ptr_src,
					 unsigned char * ptr_dst,
					 striding<3> & sr,
		   const Box<3,size_t> & bx_src)
	{
		switch (sr.n_cpy)
		{
		case 1:
				copy_grid_fast_shortx_3<object_size,1>(bx_src,ptr_dst,ptr_src,sr);
				break;
		case 2:
				copy_grid_fast_shortx_3<object_size,2>(bx_src,ptr_dst,ptr_src,sr);
				break;

		case 3:
				copy_grid_fast_shortx_3<object_size,3>(bx_src,ptr_dst,ptr_src,sr);
				break;

		case 4:
				copy_grid_fast_shortx_3<object_size,4>(bx_src,ptr_dst,ptr_src,sr);
				break;
		case 5:
				copy_grid_fast_shortx_3<object_size,5>(bx_src,ptr_dst,ptr_src,sr);
				break;
		case 6:
				copy_grid_fast_shortx_3<object_size,6>(bx_src,ptr_dst,ptr_src,sr);
				break;

		case 7:
				copy_grid_fast_shortx_3<object_size,7>(bx_src,ptr_dst,ptr_src,sr);
				break;

		case 8:
				copy_grid_fast_shortx_3<object_size,8>(bx_src,ptr_dst,ptr_src,sr);
				break;

		default:
				copy_grid_fast_longx_3<object_size>(bx_src,ptr_dst,ptr_src,sr);
		}
	}
};

template<>
struct copy_ndim_fast_selector<2>
{
	template<unsigned int object_size>
	static void call(unsigned char * ptr_src,
			 unsigned char * ptr_dst,
			 striding<2> & sr,
  const Box<2,size_t> & bx_src)
	{
		switch (sr.n_cpy)
		{
		case 1:
				copy_grid_fast_shortx_2<object_size,1>(bx_src,ptr_dst,ptr_src,sr);
				break;
		case 2:
				copy_grid_fast_shortx_2<object_size,2>(bx_src,ptr_dst,ptr_src,sr);
				break;

		case 3:
				copy_grid_fast_shortx_2<object_size,3>(bx_src,ptr_dst,ptr_src,sr);
				break;

		case 4:
				copy_grid_fast_shortx_2<object_size,4>(bx_src,ptr_dst,ptr_src,sr);
				break;
		case 5:
				copy_grid_fast_shortx_2<object_size,5>(bx_src,ptr_dst,ptr_src,sr);
				break;
		case 6:
				copy_grid_fast_shortx_2<object_size,6>(bx_src,ptr_dst,ptr_src,sr);
				break;

		case 7:
				copy_grid_fast_shortx_2<object_size,7>(bx_src,ptr_dst,ptr_src,sr);
				break;

		case 8:
				copy_grid_fast_shortx_2<object_size,8>(bx_src,ptr_dst,ptr_src,sr);
				break;

		default:
				copy_grid_fast_longx_2<object_size>(bx_src,ptr_dst,ptr_src,sr);
		}
	}
};


template<unsigned int dim_prp, unsigned int prp,typename grid_type>
struct get_pointer
{
	static unsigned char * get(const grid_type & gd, grid_key_dx<grid_type::dims> & key, unsigned int (& id)[dim_prp])
	{
		return (unsigned char *)(&gd.template get_unsafe<prp>(key));
	}
};

template<unsigned int prp,typename grid_type>
struct get_pointer<1,prp,grid_type>
{
	static unsigned char * get(const grid_type & gd, grid_key_dx<grid_type::dims> & key, unsigned int (& id)[1])
	{
		return (unsigned char *)(&gd.template get_unsafe<prp>(key)[id[0]]);
	}
};

template<unsigned int prp,typename grid_type>
struct get_pointer<2,prp,grid_type>
{
	static unsigned char * get(const grid_type & gd, grid_key_dx<grid_type::dims> & key, unsigned int (& id)[2])
	{
		return (unsigned char *)(&gd.template get_unsafe<prp>(key)[id[0]][id[1]]);
	}
};

template<unsigned int dim,unsigned int dim_prp,unsigned int prp, typename grid_type>
struct get_striding
{
	static striding<dim> get(const grid_type & gd_src,grid_type &gd_dst,
							   const Box<dim,size_t> & bx_src, unsigned int (& id)[dim_prp])
	{
		striding<dim> sr;

		grid_key_dx<dim> zero;
		zero.zero();
		grid_key_dx<dim> one = zero;
		one.set_d(1,1);

		unsigned char * ptr_final_src = get_pointer<dim_prp,prp,grid_type>::get(gd_src,one,id);
		unsigned char * ptr_start_src = get_pointer<dim_prp,prp,grid_type>::get(gd_src,zero,id);

		unsigned char * ptr_final_dst = get_pointer<dim_prp,prp,grid_type>::get(gd_dst,one,id);
		unsigned char * ptr_start_dst = get_pointer<dim_prp,prp,grid_type>::get(gd_dst,zero,id);

		sr.n_cpy = bx_src.getHigh(0) - bx_src.getLow(0) + 1;
		sr.tot_y = bx_src.getHigh(1) - bx_src.getLow(1) + 1;

		sr.striding_src[0] = ptr_final_src - ptr_start_src;
		sr.striding_dst[0] = ptr_final_dst - ptr_start_dst;

		if (dim == 3)
		{
			grid_key_dx<dim> one2 = zero;
			one2.set_d(2,1);

			unsigned char * ptr_final_src = get_pointer<dim_prp,prp,grid_type>::get(gd_src,one2,id);
			unsigned char * ptr_start_src = get_pointer<dim_prp,prp,grid_type>::get(gd_src,zero,id);

			unsigned char * ptr_final_dst = get_pointer<dim_prp,prp,grid_type>::get(gd_dst,one2,id);
			unsigned char * ptr_start_dst = get_pointer<dim_prp,prp,grid_type>::get(gd_dst,zero,id);

			sr.striding_src[1] = ptr_final_src - ptr_start_src;
			sr.striding_dst[1] = ptr_final_dst - ptr_start_dst;
		}

		return sr;
	};
};

template<unsigned int dim, typename T>
struct mp_funct_impl
{
	template<unsigned int prp, typename grid>
	static void process(const grid & gd_src,grid & gd_dst,
						const Box<dim,size_t> & bx_src, const Box<grid::dims,size_t> & bx_dst)
	{
		unsigned int id[0];
		striding<dim> sr = get_striding<dim,0,prp,grid>::get(gd_src,gd_dst,bx_src,id);

		unsigned char * ptr_src = (unsigned char *)(&gd_src.template get<prp>(bx_src.getKP1()));
		unsigned char * ptr_dst = (unsigned char *)(&gd_dst.template get<prp>(bx_dst.getKP1()));

		copy_ndim_fast_selector<dim>::template call<sizeof(T)>(ptr_src,ptr_dst,sr,bx_src);
	}
};

template<unsigned int dim, typename T, unsigned int N1>
struct mp_funct_impl<dim,T[N1]>
{
	template<unsigned int prp, typename grid>
	static void process(const grid & gd_src,grid & gd_dst,
						const Box<dim,size_t> & bx_src, const Box<grid::dims,size_t> & bx_dst)
	{
		unsigned int id[1];

		for (id[0] = 0 ; id[0] < N1 ; id[0]++)
		{
			striding<dim> sr = get_striding<dim,1,prp,grid>::get(gd_src,gd_dst,bx_src,id);

			unsigned char * ptr_src = (unsigned char *)(&gd_src.template get<prp>(bx_src.getKP1())[id[0]]);
			unsigned char * ptr_dst = (unsigned char *)(&gd_dst.template get<prp>(bx_dst.getKP1())[id[0]]);

			copy_ndim_fast_selector<dim>::template call<sizeof(T)>(ptr_src,ptr_dst,sr,bx_src);
		}
	}
};


template<unsigned int dim, typename T, unsigned int N1, unsigned int N2>
struct mp_funct_impl<dim,T[N1][N2]>
{
	template<unsigned int prp, typename grid>
	static void process(const grid & gd_src,grid & gd_dst,
						const Box<dim,size_t> & bx_src, const Box<grid::dims,size_t> & bx_dst)
	{
		unsigned int id[2];

		for (id[1] = 0 ; id[1] < N1 ; id[1]++)
		{
			for (id[0] = 0 ; id[0] < N2 ; id[0]++)
			{
				striding<dim> sr = get_striding<dim,2,prp,grid>::get(gd_src,gd_dst,bx_src,id);

				unsigned char * ptr_src = (unsigned char *)(&gd_src.template get<prp>(bx_src.getKP1())[id[1]][id[0]]);
				unsigned char * ptr_dst = (unsigned char *)(&gd_dst.template get<prp>(bx_dst.getKP1())[id[1]][id[0]]);

				copy_ndim_fast_selector<dim>::template call<sizeof(T)>(ptr_src,ptr_dst,sr,bx_src);
			}
		}
	}
};

template<typename grid, typename ginfo>
class mp_funct
{
	const grid & gd_src;
	grid & gd_dst;
	ginfo & gs_src;
	ginfo & gs_dst;

	const Box<grid::dims,size_t> & bx_src;
	const Box<grid::dims,size_t> & bx_dst;

	grid_key_dx<grid::dims> (& cnt)[1];

public:

	mp_funct(const grid & gd_src,grid & gd_dst,
			 ginfo & gs_src,ginfo & gs_dst,
			 const Box<grid::dims,size_t> & bx_src , const Box<grid::dims,size_t> & bx_dst,
			 grid_key_dx<grid::dims> (& cnt)[1])
	:gd_src(gd_src),gd_dst(gd_dst),gs_src(gs_src),gs_dst(gs_dst),bx_src(bx_src),bx_dst(bx_dst),cnt(cnt)
	{}

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename grid::value_type::type,boost::mpl::int_<T::value>>::type prp_type;

		mp_funct_impl<grid::dims,prp_type>::template process<T::value>(gd_src,gd_dst,bx_src,bx_dst);
	}
};

template<bool is_inte,unsigned int dim,typename grid, typename ginfo>
struct copy_grid_fast_layout_switch
{
	static void copy(const ginfo & gs_src,
			   const ginfo & gs_dst,
			   const Box<dim,size_t> & bx_src,
			   const Box<dim,size_t> & bx_dst,
			   const grid & gd_src,
			   grid & gd_dst,
			   grid_key_dx<dim> (& cnt)[1])
	{
		mp_funct_impl<dim,typename grid::value_type>::template process<0>(gd_src,gd_dst,bx_src,bx_dst);

//		copy_ndim_fast_selector<grid::dims>::template call<0>(gs_src,gs_dst,bx_src,bx_dst,
//							  	  	  	  	  	  	  	  	  gd_src,gd_dst,cnt);
	}
};

template<unsigned int dim, typename grid, typename ginfo>
struct copy_grid_fast_layout_switch<true,dim,grid,ginfo>
{
	static void copy(const ginfo & gs_src,
			   const ginfo & gs_dst,
			   const Box<dim,size_t> & bx_src,
			   const Box<dim,size_t> & bx_dst,
			   const grid & gd_src,
			   grid & gd_dst,
			   grid_key_dx<dim> (& cnt)[1])
	{
		mp_funct<grid,ginfo> mpf(gd_src,gd_dst,gs_src,gs_dst,bx_src,bx_dst,cnt);

		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,grid::value_type::max_prop>>(mpf);
	}
};

/*! \brief This is a way to quickly copy a grid into another grid
 *
 *
 */
template<typename grid, typename ginfo>
struct copy_grid_fast<false,3,grid,ginfo>
{
	static void copy(ginfo & gs_src,
				   ginfo & gs_dst,
				   const Box<3,size_t> & bx_src,
				   const Box<3,size_t> & bx_dst,
				   const grid & gd_src,
				   grid & gd_dst,
				   grid_key_dx<3> (& cnt)[1] )
	{
		copy_grid_fast_layout_switch<is_layout_inte<typename grid::layout_base_>::value,3,grid,ginfo>::copy(gs_src,gs_dst,bx_src,bx_dst,gd_src,gd_dst,cnt);
	}
};





/*! \brief This is a way to quickly copy a grid into another grid
 *
 *
 */
template<typename grid, typename ginfo>
struct copy_grid_fast<false,2,grid,ginfo>
{
	static void copy(ginfo & gs_src,
				   ginfo & gs_dst,
				   const Box<2,size_t> & bx_src,
				   const Box<2,size_t> & bx_dst,
				   const grid & gd_src,
				   grid & gd_dst,
				   grid_key_dx<2> (& cnt)[1] )
	{
		copy_grid_fast_layout_switch<is_layout_inte<typename grid::layout_base_>::value,2,grid,ginfo>::copy(gs_src,gs_dst,bx_src,bx_dst,gd_src,gd_dst,cnt);
	}
};

//////////////////// Pack grid fast


/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
 *
 * \tparam it type of iterator of the grid-structure
 * \tparam dtype type of the structure B
 * \tparam dim Dimensionality of the grid
 * \tparam properties to pack
 *
 */
template <bool is_complex,
		  unsigned int dim,
		  typename grid,
          typename encap_src,
		  typename encap_dst,
		  typename boost_vct,
		  typename it,
		  typename dtype,
		  int ... prp>
struct pack_with_iterator
{
	/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
	 *
	 * \param it Grid iterator
	 * \param obj object to pack
	 * \param dest where to pack
	 *
	 */
	static void pack(grid & gr, it & sub_it, dtype & dest)
	{
		size_t id = 0;

		// Packing the information
		while (sub_it.isNext())
		{
			// Copy only the selected properties
			object_si_d<encap_src,encap_dst,OBJ_ENCAP,prp...>(gr.get_o(sub_it.get()),dest.get(id));

			++id;
			++sub_it;
		}
	}
};




/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
 *
 * \tparam it type of iterator of the grid-structure
 * \tparam dtype type of the structure B
 * \tparam dim Dimensionality of the grid
 * \tparam properties to pack
 *
 */
template <typename grid,
          typename encap_src,
		  typename encap_dst,
		  typename boost_vct,
		  typename it,
		  typename dtype,
		  int ... prp>
struct pack_with_iterator<true,3,grid,encap_src,encap_dst,boost_vct,it,dtype,prp...>
{
	/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
	 *
	 * \param it Grid iterator
	 * \param obj object to pack
	 * \param dest where to pack
	 *
	 */
	static void pack(grid & gr, it & sub_it, dtype & dest)
	{
		size_t id = 0;

		size_t lin_src = 0;

		auto & gs_src = gr.getGrid();
		grid_key_dx<3> start = sub_it.getStart();
		grid_key_dx<3> stop = sub_it.getStop();

		lin_src += start.get(2) * gs_src.size_s(1);
		for (long int i = start.get(2) ; i <= stop.get(2) ; i++)
		{
			lin_src += start.get(1) * gs_src.size_s(0);
			for (long int j = start.get(1) ; j <= stop.get(1) ; j++)
			{
				lin_src += start.get(0);
				for (long int k = start.get(0) ; k <= stop.get(0) ; k++)
				{
					// Copy only the selected properties
					object_si_d<encap_src,encap_dst,OBJ_ENCAP,prp...>(gr.get_o(lin_src),dest.get(id));

					++id;
					++lin_src;
				}
				lin_src -= stop.get(0) + 1;
				lin_src += gs_src.size_s(0);
			}
			lin_src -= (stop.get(1) + 1)*gs_src.size_s(0);
			lin_src += gs_src.size_s(1);
		}
	}
};

/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
 *
 * \tparam it type of iterator of the grid-structure
 * \tparam dtype type of the structure B
 * \tparam dim Dimensionality of the grid
 * \tparam properties to pack
 *
 */
template <typename grid,
          typename encap_src,
		  typename encap_dst,
		  typename boost_vct,
		  typename it,
		  typename dtype,
		  int ... prp>
struct pack_with_iterator<true,2,grid,encap_src,encap_dst,boost_vct,it,dtype,prp...>
{
	/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
	 *
	 * \param it Grid iterator
	 * \param obj object to pack
	 * \param dest where to pack
	 *
	 */
	static void pack(grid & gr, it & sub_it, dtype & dest)
	{
		size_t id = 0;

		size_t lin_src = 0;

		auto & gs_src = gr.getGrid();
		grid_key_dx<2> start = sub_it.getStart();
		grid_key_dx<2> stop = sub_it.getStop();

		lin_src += start.get(1) * gs_src.size_s(0);
		for (long int j = start.get(1) ; j <= stop.get(1) ; j++)
		{
			lin_src += start.get(0);
			for (long int k = start.get(0) ; k <= stop.get(0) ; k++)
			{
				// Copy only the selected properties
				object_si_d<encap_src,encap_dst,OBJ_ENCAP,prp...>(gr.get_o(lin_src),dest.get(id));

				++id;
				++lin_src;
			}
			lin_src -= stop.get(0) + 1;
			lin_src += gs_src.size_s(0);
		}

	}
};

//////////// PACK ITERATOR LONG /////////////////////////

template<unsigned int dim, int obj_byte, typename git, typename grid>
struct pack_with_iterator_longx
{
	static void pack(grid & gr,
			  git & sub_it,
			  unsigned char * ptr_dest,
			  unsigned char * ptr,
			  size_t stride_x,
			  size_t n_cpy)
	{
		std::cout << __FILE__ << ":" << __LINE__ << " critical error, we shoukd never be here" << std::endl;
	}
};

template<int obj_byte, typename git, typename grid>
struct pack_with_iterator_longx<2,obj_byte,git,grid>
{
	static void pack(grid & gr,
			  git & sub_it,
			  unsigned char * ptr_dest,
			  unsigned char * ptr,
			  size_t stride_x,
			  size_t n_cpy)
	{
		size_t tot_y = sub_it.getStop().get(1) - sub_it.getStart().get(1) + 1;

		for (size_t i = 0 ; i < tot_y ; i++)
		{
			memcpy(ptr_dest,ptr,n_cpy * obj_byte);

			ptr += stride_x;
			ptr_dest += n_cpy * obj_byte;
		}
	}
};

template<int obj_byte, typename git, typename grid>
struct pack_with_iterator_longx<3,obj_byte,git,grid>
{
	static void pack(grid & gr,
			  git & sub_it,
			  unsigned char * ptr_dest,
			  unsigned char * ptr,
			  size_t stride_x,
			  size_t n_cpy)
	{
		size_t tot_y = sub_it.getStop().get(1) - sub_it.getStart().get(1) + 1;
		size_t tot_z = sub_it.getStop().get(2) - sub_it.getStart().get(2) + 1;

		grid_key_dx<3> zero;
		zero.zero();
		grid_key_dx<3> one = zero;
		one.set_d(2,1);

		unsigned char * ptr_final = (unsigned char *)&(gr.template get<0>(one));
		unsigned char * ptr_start = (unsigned char *)&(gr.template get<0>(zero));

		size_t stride_y = ptr_final - ptr_start;

		for (size_t i = 0 ; i < tot_z ; i++)
		{
			for (size_t i = 0 ; i < tot_y ; i++)
			{
				memcpy(ptr_dest,ptr,n_cpy * obj_byte);

				ptr += stride_x;
				ptr_dest += n_cpy * obj_byte;
			}
			ptr += stride_y - tot_y*stride_x;
		}
	}
};

/////////////////////////////////////////////////////////////////////

///////////////// PACK ITERATOR SHORT ///////////////////////////////

template<unsigned int dim, unsigned int n_cpy ,int obj_byte, typename git, typename grid>
struct pack_with_iterator_shortx
{
	static void pack(grid & gr,
			  git & sub_it,
			  unsigned char * ptr_dest,
			  unsigned char * ptr,
			  size_t stride_x)
	{
		std::cout << __FILE__ << ":" << __LINE__ << " critical error, we shoukd never be here" << std::endl;
	}
};

template<unsigned int n_cpy, int obj_byte, typename git, typename grid>
struct pack_with_iterator_shortx<2,n_cpy,obj_byte,git,grid>
{
	static void pack(grid & gr,
			  git & sub_it,
			  unsigned char * ptr_dest,
			  unsigned char * ptr,
			  size_t stride_x)
	{
		size_t tot_y = sub_it.getStop().get(1) - sub_it.getStart().get(1) + 1;

		for (size_t i = 0 ; i < tot_y ; i++)
		{
			__builtin_memcpy(ptr_dest,ptr,n_cpy * obj_byte);

			ptr += stride_x;
			ptr_dest += n_cpy * obj_byte;
		}
	}
};

template<unsigned int n_cpy, int obj_byte, typename git, typename grid>
struct pack_with_iterator_shortx<3,n_cpy,obj_byte,git,grid>
{
	static void pack(grid & gr,
			  git & sub_it,
			  unsigned char * ptr_dest,
			  unsigned char * ptr,
			  size_t stride_x)
	{
		size_t tot_y = sub_it.getStop().get(1) - sub_it.getStart().get(1) + 1;
		size_t tot_z = sub_it.getStop().get(2) - sub_it.getStart().get(2) + 1;

		grid_key_dx<3> zero;
		zero.zero();
		grid_key_dx<3> one = zero;
		one.set_d(2,1);

		unsigned char * ptr_final = (unsigned char *)(&gr.template get_unsafe<0>(one));
		unsigned char * ptr_start = (unsigned char *)(&gr.template get<0>(zero));

		size_t stride_y = ptr_final - ptr_start;

		for (size_t i = 0 ; i < tot_z ; i++)
		{
			for (size_t i = 0 ; i < tot_y ; i++)
			{
				__builtin_memcpy(ptr_dest,ptr,n_cpy * obj_byte);

				ptr += stride_x;
				ptr_dest += n_cpy * obj_byte;
			}
			ptr += stride_y - tot_y*stride_x;
		}
	}
};

///////////////////////////////////////////////////////////////////

/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
 *
 * \tparam it type of iterator of the grid-structure
 * \tparam dtype type of the structure B
 * \tparam dim Dimensionality of the grid
 * \tparam properties to pack
 *
 */
template <unsigned int dim,
		  typename grid,
          typename encap_src,
		  typename encap_dst,
		  typename boost_vct,
		  typename it,
		  typename dtype,
		  int ... prp>
struct pack_with_iterator<false,dim,grid,encap_src,encap_dst,boost_vct,it,dtype,prp...>
{
	/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
	 *
	 * \param it Grid iterator
	 * \param obj object to pack
	 * \param dest where to pack
	 *
	 */
	static void pack(grid & gr, it & sub_it, dtype & dest)
	{
		// We do not have an optimized version for dimension different from 3 and 2
		if (dim == 1 || dim > 3)
		{
			pack_with_iterator<true,dim,grid,encap_src,encap_dst,boost_vct,it,dtype,prp...>::pack(gr,sub_it,dest);
			return;
		}

		// Sending property object
		typedef object<typename object_creator<
										boost_vct,
										prp...>::type
					  > prp_object;


		// fast strided copy

		grid_key_dx<dim> zero;
		zero.zero();
		grid_key_dx<dim> one = zero;
		one.set_d(1,1);

		unsigned char * ptr_final = (unsigned char *)(&gr.template get_unsafe<0>(one));
		unsigned char * ptr_start = (unsigned char *)(&gr.template get<0>(zero));

		size_t stride = ptr_final - ptr_start;

		size_t n_tot = 0;
		for (size_t i = 1 ; i < dim ; i++)
		{n_tot += sub_it.getStop().get(i) - sub_it.getStart().get(i) + 1;}

		size_t n_cpy = sub_it.getStop().get(0) - sub_it.getStart().get(0) + 1;
		unsigned char * ptr = (unsigned char *)&(gr.template get<first_variadic<prp...>::type::value>(sub_it.getStart()));
		unsigned char * ptr_dest = (unsigned char *)dest.getPointer();

		switch (n_cpy)
		{
		case 1:
			pack_with_iterator_shortx<dim,1,sizeof(prp_object),it,grid>::pack(gr,
															sub_it,
															ptr_dest,
															ptr,
															stride);
			break;
		case 2:
			pack_with_iterator_shortx<dim,2,sizeof(prp_object),it,grid>::pack(gr,
															sub_it,
															ptr_dest,
															ptr,
															stride);
			break;
		case 3:
			pack_with_iterator_shortx<dim,3,sizeof(prp_object),it,grid>::pack(gr,
															sub_it,
															ptr_dest,
															ptr,
															stride);
			break;
		case 4:
			pack_with_iterator_shortx<dim,4,sizeof(prp_object),it,grid>::pack(gr,
															sub_it,
															ptr_dest,
															ptr,
															stride);
			break;
		case 5:
			pack_with_iterator_shortx<dim,5,sizeof(prp_object),it,grid>::pack(gr,
															sub_it,
															ptr_dest,
															ptr,
															stride);
			break;
		case 6:
			pack_with_iterator_shortx<dim,6,sizeof(prp_object),it,grid>::pack(gr,
															sub_it,
															ptr_dest,
															ptr,
															stride);
			break;
		case 7:
			pack_with_iterator_shortx<dim,7,sizeof(prp_object),it,grid>::pack(gr,
															sub_it,
															ptr_dest,
															ptr,
															stride);
			break;
		case 8:
			pack_with_iterator_shortx<dim,8,sizeof(prp_object),it,grid>::pack(gr,
															sub_it,
															ptr_dest,
															ptr,
															stride);
			break;
		default:
			pack_with_iterator_longx<dim,sizeof(prp_object),it,grid>::pack(gr,
															sub_it,
															ptr_dest,
															ptr,
															stride,
															n_cpy);

		}

	}
};

/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
 *
 * \tparam it type of iterator of the grid-structure
 * \tparam dtype type of the structure B
 * \tparam dim Dimensionality of the grid
 * \tparam properties to pack
 *
 */
template <typename grid,
          typename encap_src,
		  typename encap_dst,
		  typename boost_vct,
		  typename it,
		  typename dtype,
		  int ... prp>
struct pack_with_iterator<true,1,grid,encap_src,encap_dst,boost_vct,it,dtype,prp...>
{
	/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
	 *
	 * \param it Grid iterator
	 * \param obj object to pack
	 * \param dest where to pack
	 *
	 */
	static void pack(grid & gr, it & sub_it, dtype & dest)
	{
		size_t id = 0;

		size_t lin_src = 0;

		auto & gs_src = gr.getGrid();
		grid_key_dx<1> start = sub_it.getStart();
		grid_key_dx<1> stop = sub_it.getStop();


		lin_src += start.get(0);
		for (long int k = start.get(0) ; k <= stop.get(0) ; k++)
		{
			// Copy only the selected properties
			object_si_d<encap_src,encap_dst,OBJ_ENCAP,prp...>(gr.get_o(lin_src),dest.get(id));

			++id;
			++lin_src;
		}
		lin_src -= stop.get(0) + 1;
		lin_src += gs_src.size_s(0);

	}
};

//////////////////////////// Unpack grid fast ////////////////////////////

/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
 *
 * \tparam it type of iterator of the grid-structure
 * \tparam dtype type of the structure B
 * \tparam dim Dimensionality of the grid
 * \tparam properties to pack
 *
 */
template <unsigned int dim,
		  typename grid,
          typename encap_src,
		  typename encap_dst,
		  typename boost_vct,
		  typename it,
		  typename stype,
		  int ... prp>
struct unpack_with_iterator
{
	/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
	 *
	 * \param it Grid iterator
	 * \param obj object to pack
	 * \param dest where to pack
	 *
	 */
	static void unpack(grid & gr, it & sub_it, stype & src)
	{
		size_t id = 0;

		// unpacking the information
		while (sub_it.isNext())
		{

			// Copy only the selected properties
			object_s_di<encap_src,encap_dst,OBJ_ENCAP,prp...>(src.get(id),gr.get_o(sub_it.get()));

			++id;
			++sub_it;
		}
	}
};


/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
 *
 * \tparam it type of iterator of the grid-structure
 * \tparam dtype type of the structure B
 * \tparam dim Dimensionality of the grid
 * \tparam properties to pack
 *
 */
template <typename grid,
          typename encap_src,
		  typename encap_dst,
		  typename boost_vct,
		  typename it,
		  typename stype,
		  int ... prp>
struct unpack_with_iterator<3,grid,
							encap_src,
							encap_dst,
							boost_vct,
							it,
							stype,
							prp ...>
{
	/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
	 *
	 * \param it Grid iterator
	 * \param obj object to pack
	 * \param dest where to pack
	 *
	 */
	static void unpack(grid & gr, it & sub_it, stype & src)
	{
		size_t id = 0;

		size_t lin_dst = 0;

		auto & gs_dst = gr.getGrid();
		grid_key_dx<3> start = sub_it.getStart();
		grid_key_dx<3> stop = sub_it.getStop();

		// unpacking the information

		lin_dst += start.get(2) * gs_dst.size_s(1);
		for (long int i = start.get(2) ; i <= stop.get(2) ; i++)
		{
			lin_dst += start.get(1) * gs_dst.size_s(0);
			for (long int j = start.get(1) ; j <= stop.get(1) ; j++)
			{
				lin_dst += start.get(0);
				for (long int k = start.get(0) ; k <= stop.get(0) ; k++)
				{
					// Copy only the selected properties
					object_s_di<encap_src,encap_dst,OBJ_ENCAP,prp...>(src.get(id),gr.get_o(lin_dst));

					++id;
					++lin_dst;
				}
				lin_dst -= stop.get(0) + 1;
				lin_dst += gs_dst.size_s(0);
			}
			lin_dst -= (stop.get(1) + 1)*gs_dst.size_s(0);
			lin_dst += gs_dst.size_s(1);
		}
	}
};

#endif /* OPENFPM_DATA_SRC_GRID_COPY_GRID_FAST_HPP_ */
