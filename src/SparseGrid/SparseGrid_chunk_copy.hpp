/*
 * SparseGrid_chunk_copy.hpp
 *
 *  Created on: Mar 18, 2020
 *      Author: i-bird
 */

#ifndef SPARSEGRID_CHUNK_COPY_HPP_
#define SPARSEGRID_CHUNK_COPY_HPP_

//! Linearize a set of index
template<typename vmpl, typename a> __device__ __host__ inline size_t Lin_vmpl(a v)
{
	return v*vmpl_reduce_prod_stop<vmpl,(int)vmpl::size::value - 2>::type::value;
}

/*! \brief linearize an arbitrary set of index
 *
 * linearize an arbitrary set of index
 *
 */
template<typename vmpl, typename a, typename ...lT>
__device__ __host__ inline size_t Lin_vmpl(a v,lT...t)
{
		return v*vmpl_reduce_prod_stop<vmpl,(int)vmpl::size::value - sizeof...(t) - 2>::type::value + Lin_vmpl<vmpl>(t...);
}




//! Copy XY surface in 3D
template<int layout_type, int prop,int stencil_size, typename chunking,bool is_cross>
struct copy_xy_3
{
	template<unsigned int i_src, unsigned int i_dest, unsigned int N1, typename T, typename headerType, typename chunkType>
	inline static void copy(T ptr[N1], unsigned char mask[N1], headerType & h , const chunkType & chunk)
	{
		int s = stencil_size + stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		for (int v = 0 ; v < stencil_size ; v++)
		{
		for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
		{
			for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
			{
				const int id = Lin_vmpl<typename chunking::type>(k,j,i_src+v);

				ptr[s] = chunk.template get<prop>()[id];
				mask[s] = exist_sub(h,id);

				s++;
			}
			s += 2*stencil_size;
		}

		s+= 2*stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
		}
	}
};


//! Copy XY surface in 3D
template<int layout_type, int prop , typename chunking,bool is_cross>
struct copy_xy_3<layout_type,prop,1,chunking,is_cross>
{
	template<unsigned int i_src, unsigned int i_dest, unsigned int N1, typename T, typename headerType,  typename chunkType>
	inline static void copy(T ptr[N1],  unsigned char mask[N1],const headerType & h,  const chunkType & chunk)
	{
		int s = 1 + 1*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);

		for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
		{
			for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
			{
				const int id = Lin_vmpl<typename chunking::type>(k,j,i_src);

				ptr[s] = chunk.template get<prop>()[id];
				mask[s] = exist_sub(h,id);

				s++;
			}
			s += 2*1;
		}
	}
};

//! Copy XZ surface in 3D
template<int layout_type, int prop, int stencil_size, typename chunking,bool is_cross>
struct copy_xz_3
{
	template<unsigned int j_src, unsigned int j_dest, unsigned int N1, typename T, typename headerType , typename chunkType>
	inline static void copy(T ptr[N1], unsigned char mask[N1], headerType & h , const chunkType & chunk)
	{
		int s = stencil_size + j_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		for (int v = 0 ; v < stencil_size ; v++)
		{
		for (int i = 0 ; i < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; i++)
		{
			for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
			{
				const int id = Lin_vmpl<typename chunking::type>(k,j_src+v,i);

				ptr[s] = chunk.template get<prop>()[id];
				mask[s] = exist_sub(h,id);

				s++;
			}

			s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value + 2*stencil_size) - boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value;
		}

		s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value + 2*stencil_size) - boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);
		}
	}
};

//! Copy XZ surface in 3D
template<int layout_type, int prop, typename chunking,bool is_cross>
struct copy_xz_3<layout_type,prop,1,chunking,is_cross>
{
	template<unsigned int j_src, unsigned int j_dest, unsigned int N1, typename T, typename headerType , typename chunkType>
	inline static void copy(T ptr[N1], unsigned char mask[N1], headerType & h , const chunkType & chunk)
	{
		int s = 1 + j_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1) +
				1*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);

		for (int i = 0 ; i < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; i++)
		{
			for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
			{
				const int id = Lin_vmpl<typename chunking::type>(k,j_src,i);

				ptr[s] = chunk.template get<prop>()[id];
				mask[s] = exist_sub(h,id);

				s++;
			}

			s += 1 * (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value + 2*1) - boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value;
		}
	}
};

//! Copy YZ surface in 3D
template<int layout_type, int prop, int stencil_size, typename chunking,bool is_cross>
struct copy_yz_3
{
	template<unsigned int k_src, unsigned int k_dest, unsigned int N1, typename T, typename headerType , typename chunkType>
	inline static void copy(T ptr[N1], unsigned char mask[N1], headerType & h , const chunkType & chunk)
	{
		int s = k_dest + stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		for (int v = 0 ; v < stencil_size ; v++)
		{
		for (int i = 0 ; i < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; i++)
		{
			for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
			{
				const int id = Lin_vmpl<typename chunking::type>(k_src+v,j,i);

				ptr[s] = chunk.template get<prop>()[id];
				mask[s] = exist_sub(h,id);

				s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
			}

			s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value + 2*stencil_size) - boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
		}

		s += 1 - boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		}
	}
};

//! Copy YZ surface in 3D
template<int layout_type, int prop, typename chunking,bool is_cross>
struct copy_yz_3<layout_type,prop,1,chunking,is_cross>
{
	template<unsigned int k_src, unsigned int k_dest, unsigned int N1, typename T, typename headerType, typename chunkType>
	inline static void copy(T ptr[N1], unsigned char mask[N1], headerType & h , const chunkType & chunk)
	{
		int s = k_dest + 1*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1) +
				1*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);

		for (int i = 0 ; i < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; i++)
		{
			for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
			{
				const int id = Lin_vmpl<typename chunking::type>(k_src,j,i);

				ptr[s] = chunk.template get<prop>()[id];
				exist_sub(h,id);

				s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1);
			}

			s += 1 * (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value + 2*1) - boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1);
		}
	}
};

//! Copy x edge in 3D
template<int layout_type, int prop, int stencil_size, typename chunking,bool is_cross>
struct copy_x_3
{
	template<unsigned int i_src, unsigned int i_dest,unsigned int j_src, unsigned int j_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{
		int s = stencil_size + j_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		for (int v1 = 0 ; v1 < stencil_size ; v1++)
		{
			for (int v2 = 0 ; v2 < stencil_size ; v2++)
			{
				for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
				{
					ptr[s] = chunk.template get<prop>()[Lin_vmpl<typename chunking::type>(k,j_src+v2,i_src+v1)];

					s++;
				}

				s+= 2*stencil_size;
			}

			s+= (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size) - stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
		}
	}
};

//! Copy x edge in 3D
template<int layout_type, int prop, typename chunking,bool is_cross>
struct copy_x_3<layout_type,prop,1,chunking,is_cross>
{
	template<unsigned int i_src, unsigned int i_dest,unsigned int j_src, unsigned int j_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{
		int s = 1 + j_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);


		for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
		{
			ptr[s] = chunk.template get<prop>()[Lin_vmpl<typename chunking::type>(k,j_src,i_src)];

			s++;
		}
	}
};

//! Copy x edge in 3D
template<int layout_type, int prop, int stencil_size, typename chunking>
struct copy_x_3<layout_type,prop,stencil_size,chunking,true>
{
	template<unsigned int i_src, unsigned int i_dest,unsigned int j_src, unsigned int j_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{}
};

//! Copy y edge in 3D
template<int layout_type, int prop, int stencil_size, typename chunking,bool is_cross>
struct copy_y_3
{
	template<unsigned int i_src, unsigned int i_dest,unsigned int k_src, unsigned int k_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{
		int s = k_dest + stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		for (int v1 = 0 ; v1 < stencil_size ; v1++)
		{
			for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
			{
				for (int v2 = 0 ; v2 < stencil_size ; v2++)
				{
					ptr[s] = chunk.template get<prop>()[Lin_vmpl<typename chunking::type>(k_src+v2,j,i_src+v1)];
					s+= 1;
				}

				s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) - stencil_size;
			}

			s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size) - boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
		}
	}
};

//! Copy y edge in 3D
template<int layout_type, int prop, typename chunking,bool is_cross>
struct copy_y_3<layout_type,prop,1,chunking,is_cross>
{
	template<unsigned int i_src, unsigned int i_dest,unsigned int k_src, unsigned int k_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{
		int s = k_dest + 1*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);


		for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; j++)
		{
			ptr[s] = chunk.template get<prop>()[Lin_vmpl<typename chunking::type>(k_src,j,i_src)];

			s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1);
		}
	}
};

//! Copy y edge in 3D
template<int layout_type, int prop, int stencil_size, typename chunking>
struct copy_y_3<layout_type,prop,stencil_size,chunking,true>
{
	template<unsigned int i_src, unsigned int i_dest,unsigned int k_src, unsigned int k_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{}
};


//! Copy z edge in 3D
template<int layout_type, int prop, int stencil_size, typename chunking,bool is_cross>
struct copy_z_3
{
	template<unsigned int j_src, unsigned int j_dest,unsigned int k_src, unsigned int k_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{
		int s = k_dest + j_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);


		for (int i = 0 ; i < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; i++)
		{
			for (int v1 = 0 ; v1 < stencil_size ; v1++)
			{
					for (int v2 = 0 ; v2 < stencil_size ; v2++)
					{
						ptr[s] = chunk.template get<prop>()[Lin_vmpl<typename chunking::type>(k_src+v2,j_src+v1,i)];

						s++;
					}

					s+= 2*stencil_size;
			}

			s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);
		}
	}
};

//! Copy z edge in 3D
template<int layout_type, int prop, typename chunking,bool is_cross>
struct copy_z_3<layout_type,prop,1,chunking,is_cross>
{
	template<unsigned int j_src, unsigned int j_dest,unsigned int k_src, unsigned int k_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{
		int s = k_dest + j_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1) +
				1*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);


		for (int i = 0 ; i < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; i++)
		{
			ptr[s] = chunk.template get<prop>()[Lin_vmpl<typename chunking::type>(k_src,j_src,i)];

			s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);
		}
	}
};

//! Copy z edge in 3D
template<int layout_type, int prop, int stencil_size, typename chunking>
struct copy_z_3<layout_type,prop,stencil_size,chunking,true>
{
	template<unsigned int j_src, unsigned int j_dest,unsigned int k_src, unsigned int k_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{}
};

//! Copy point in 3D
template<int layout_type, int prop, int stencil_size, typename chunking,bool is_cross>
struct copy_corner_3
{
	template<unsigned int i_src, unsigned int i_dest,unsigned int j_src, unsigned int j_dest,unsigned int k_src, unsigned int k_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{
		int s = k_dest + j_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);


		for (int i = 0 ; i < stencil_size ; i++)
		{
			for (int j = 0 ; j < stencil_size ; j++)
			{
				for (int k = 0 ; k < stencil_size ; k++)
				{
					ptr[s] = chunk.template get<prop>()[Lin_vmpl<typename chunking::type>(k_src+k,j_src+j,i_src+i)];

					s++;
				}

				s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) - stencil_size;
			}

			s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size) - stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
		}
	}
};

//! Copy point in 3D
template<int layout_type, int prop, typename chunking,bool is_cross>
struct copy_corner_3<layout_type,prop,1,chunking,is_cross>
{
	template<unsigned int i_src, unsigned int i_dest,unsigned int j_src, unsigned int j_dest,unsigned int k_src, unsigned int k_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{
		int s = k_dest + j_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);


		ptr[s] = chunk.template get<prop>()[Lin_vmpl<typename chunking::type>(k_src,j_src,i_src)];
	}
};

//! Copy point in 3D
template<int layout_type, int prop, int stencil_size, typename chunking>
struct copy_corner_3<layout_type,prop,stencil_size,chunking,true>
{
	template<unsigned int i_src, unsigned int i_dest,unsigned int j_src, unsigned int j_dest,unsigned int k_src, unsigned int k_dest , unsigned int N1, typename T, typename chunkType>
	inline static void copy(T ptr[N1], const chunkType & chunk)
	{}
};

#endif /* SPARSEGRID_CHUNK_COPY_HPP_ */
