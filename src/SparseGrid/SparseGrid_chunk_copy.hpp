/*
 * SparseGrid_chunk_copy.hpp
 *
 *  Created on: Mar 18, 2020
 *      Author: i-bird
 */

#ifndef SPARSEGRID_CHUNK_COPY_HPP_
#define SPARSEGRID_CHUNK_COPY_HPP_

#include <Vc/Vc>
#include <immintrin.h>

/*! \brief Check if the point in the chunk exist
 *
 * \param h header
 * \param sub_id index of the sub-domain
 *
 * \return true if exist, false if does not
 *
 */
template<typename headerType>
inline bool exist_sub(headerType & h, int sub_id)
{
	size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

	return h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check;
}

template<unsigned int v>
struct exist_sub_v_impl
{
	typedef unsigned char type;

	template<typename headerType>
	static inline void exist(headerType & h, int sub_id, unsigned char * pmask)
	{
		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		pmask[0] = h.mask[sub_id >> BIT_SHIFT_SIZE_T] & mask_check;
	}
};

template<>
struct exist_sub_v_impl<2>
{
	typedef unsigned short type;

	template<typename headerType>
	static inline void exist(headerType & h, int sub_id, unsigned char * pmask)
	{
		size_t m = h.mask[sub_id >> BIT_SHIFT_SIZE_T];

		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		pmask[0] = m & mask_check;
		pmask[1] = m & mask_check << 1;
	}
};

template<>
struct exist_sub_v_impl<4>
{
	typedef unsigned int type;

	template<typename headerType>
	static inline void exist(headerType & h, int sub_id, unsigned char * pmask)
	{
		size_t m = h.mask[sub_id >> BIT_SHIFT_SIZE_T];

		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		pmask[0] = (m & mask_check) != 0;
		pmask[1] = (m & (mask_check << 1)) != 0;
		pmask[2] = (m & (mask_check << 2)) != 0;
		pmask[3] = (m & (mask_check << 3)) != 0;
	}
};

template<>
struct exist_sub_v_impl<8>
{
	typedef unsigned long int type;

	template<typename headerType>
	static inline void exist(headerType & h, int sub_id, unsigned char * pmask)
	{
		size_t m = h.mask[sub_id >> BIT_SHIFT_SIZE_T];

		size_t mask_check = (size_t)1 << (sub_id & ((1 << BIT_SHIFT_SIZE_T) - 1));

		pmask[0] = m & mask_check;
		pmask[1] = m & mask_check << 1;
		pmask[2] = m & mask_check << 2;
		pmask[3] = m & mask_check << 3;
		pmask[4] = m & mask_check << 4;
		pmask[5] = m & mask_check << 5;
		pmask[6] = m & mask_check << 6;
		pmask[7] = m & mask_check << 7;
	}
};

/*! \brief Check if the point in the chunk exist (Vectorial form)
 *
 * \param h header
 * \param sub_id index of the sub-domain
 *
 * \return true if exist, false if does not
 *
 */
template<unsigned int v, typename headerType>
inline void exist_sub_v(headerType & h, int sub_id, unsigned char * pmask)
{
	exist_sub_v_impl<v>::exist(h,sub_id,pmask);
}

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


//! Copy block in 3D
template<int layout_type, int prop, int stencil_size ,typename chunking,bool is_cross>
struct copy_xyz
{
	template<unsigned int N1, typename T, typename headerType, typename chunkType>
	inline static void copy(T ptr[N1], unsigned char mask[N1], headerType & h , const chunkType & chunk)
	{
		int s = stencil_size + stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		int s2 = 0;

		for (int v = 0 ; v < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; v++)
		{
			for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
			{
				for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
				{
					ptr[s] = chunk[s2];
					mask[s] = exist_sub(h,s2);

					s++;
					s2++;
				}
				s += 2*stencil_size;
			}
		s+= 2*stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
		}
	}

	template<unsigned int N1, typename T, typename chunkType>
	inline static void store(T ptr[N1] , chunkType & chunk)
	{
		int s2 = 0;

		for (int v = 0 ; v < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; v++)
		{
			for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
			{
				for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
				{
					chunk[s2] = ptr[s2];

					s2++;
				}
			}
		}
	}
};

template<unsigned int i>
struct multi_mask
{};

template<>
struct multi_mask<1>
{
	typedef unsigned char type;
};

template<>
struct multi_mask<2>
{
	typedef unsigned short int type;
};

template<>
struct multi_mask<4>
{
	typedef unsigned int type;
};

template<>
struct multi_mask<8>
{
	typedef unsigned long int type;
};

template<>
struct multi_mask<16>
{
	typedef unsigned long int type;
};

//! Copy block in 3D vectorized
template<int prop, int stencil_size ,typename chunking,bool is_cross>
struct copy_xyz<1,prop,stencil_size,chunking,is_cross>
{
	template<unsigned int N1, typename T, typename headerType, typename chunkType>
	inline static void copy(T ptr[N1], unsigned char mask[N1], headerType & h , const chunkType & chunk)
	{
		int s = stencil_size + stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		typedef boost::mpl::int_<boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value * sizeof(chunk[0]) /
				                                Vc::float_v::Size / sizeof(float)> n_it_lead;

		int s2 = 0;

		for (int v = 0 ; v < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; v++)
		{
			for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
			{
				for (int k = 0 ; k < n_it_lead::value ; k+=4)
				{
					exist_sub_v<Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))>(h,s2,&mask[s]);
					exist_sub_v<Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))>(h,s2+Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float)),&mask[s+Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))]);
					exist_sub_v<Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))>(h,s2+2*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float)),&mask[s+2*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))]);
					exist_sub_v<Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))>(h,s2+3*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float)),&mask[s+3*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))]);

					Vc::float_v tmp = Vc::float_v((float *)&chunk[s2],Vc::Aligned);
					Vc::float_v tmp2 = Vc::float_v((float *)&chunk[s2+Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Aligned);
					Vc::float_v tmp3 = Vc::float_v((float *)&chunk[s2+2*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Aligned);
					Vc::float_v tmp4 = Vc::float_v((float *)&chunk[s2+3*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Aligned);
					tmp.store((float *)&ptr[s],Vc::Unaligned);
					tmp2.store((float *)&ptr[s+Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Unaligned);
					tmp3.store((float *)&ptr[s+2*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Unaligned);
					tmp4.store((float *)&ptr[s+3*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Unaligned);

					s += 4*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float));
					s2 += 4*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float));
				}
				s += 2*stencil_size ;
			}
		s+= 2*stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
		}
	}

	template<unsigned int N1, typename T, typename chunkType>
	inline static void store(T ptr[N1] , chunkType & chunk)
	{
		typedef boost::mpl::int_<boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value * sizeof(chunk[0]) /
				                                Vc::float_v::Size / sizeof(float)> n_it_lead;

		int s2 = 0;

		for (int v = 0 ; v < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; v++)
		{
			for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
			{
				for (int k = 0 ; k < n_it_lead::value ; k += 4)
				{
					Vc::float_v tmp = Vc::float_v((float *)&ptr[s2],Vc::Aligned);
					Vc::float_v tmp1 = Vc::float_v((float *)&ptr[s2+1*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Aligned);
					Vc::float_v tmp2 = Vc::float_v((float *)&ptr[s2+2*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Aligned);
					Vc::float_v tmp3 = Vc::float_v((float *)&ptr[s2+3*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Aligned);

					tmp.store((float *)&chunk[s2],Vc::Aligned);
					tmp1.store((float *)&chunk[s2+1*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Aligned);
					tmp2.store((float *)&chunk[s2+2*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Aligned);
					tmp3.store((float *)&chunk[s2+3*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float))],Vc::Aligned);

					s2 += 4*Vc::float_v::Size / (sizeof(chunk[0])/sizeof(float));
				}
			}
		}
	}
};

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

	template<unsigned int i_dest, unsigned int N1>
	inline static void mask_null(unsigned char mask[N1])
	{
		int s = stencil_size + stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		for (int v = 0 ; v < stencil_size ; v++)
		{
		for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
		{
			for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
			{
				mask[s] = 0;

				s++;
			}
			s += 2*stencil_size;
		}

		s+= 2*stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
		}
	}
};

//! Copy XY surface in 3D
template<int prop,int stencil_size, typename chunking,bool is_cross>
struct copy_xy_3<1,prop,stencil_size,chunking,is_cross>
{
	template<unsigned int i_src, unsigned int i_dest, unsigned int N1, typename T, typename headerType, typename chunkType>
	inline static void copy(T ptr[N1], unsigned char mask[N1], headerType & h , const chunkType & chunk)
	{
		int s = stencil_size + stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		int s2 = i_src*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value);

		typedef boost::mpl::int_<boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value * sizeof(chunk.template get<prop>()[0]) /
				                                Vc::float_v::Size / sizeof(float)> n_it_lead;

		for (int v = 0 ; v < stencil_size ; v++)
		{
		for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
		{
			for (int k = 0 ; k < n_it_lead::value ; k+=4)
			{
				exist_sub_v<Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))>(h,s2,&mask[s]);
				exist_sub_v<Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))>(h,s2+1*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float)),&mask[s+1*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))]);
				exist_sub_v<Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))>(h,s2+2*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float)),&mask[s+2*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))]);
				exist_sub_v<Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))>(h,s2+3*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float)),&mask[s+3*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))]);


				Vc::float_v tmp = Vc::float_v((float *)&chunk.template get<prop>()[s2],Vc::Unaligned);
				Vc::float_v tmp1 = Vc::float_v((float *)&chunk.template get<prop>()[s2+1*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))],Vc::Unaligned);
				Vc::float_v tmp2 = Vc::float_v((float *)&chunk.template get<prop>()[s2+2*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))],Vc::Unaligned);
				Vc::float_v tmp3 = Vc::float_v((float *)&chunk.template get<prop>()[s2+3*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))],Vc::Unaligned);
				tmp.store((float *)&ptr[s],Vc::Unaligned);
				tmp1.store((float *)&ptr[s+1*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))],Vc::Unaligned);
				tmp2.store((float *)&ptr[s+2*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))],Vc::Unaligned);
				tmp3.store((float *)&ptr[s+3*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float))],Vc::Unaligned);

				s += 4*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float));
				s2 += 4*Vc::float_v::Size / (sizeof(chunk.template get<prop>()[0])/sizeof(float));
			}
			s += 2*stencil_size;

		}

		s+= 2*stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
		}
	}

	template<unsigned int i_dest, unsigned int N1>
	inline static void mask_null(unsigned char mask[N1])
	{
		int s = stencil_size + stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size) +
				i_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*stencil_size);

		for (int v = 0 ; v < stencil_size ; v++)
		{
		for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
		{
			for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
			{
				mask[s] = 0;

				s++;
			}
			s += 2*stencil_size;
		}

		s+= 2*stencil_size*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*stencil_size);
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

	template<unsigned int j_dest, unsigned int N1>
	inline static void mask_null(unsigned char mask[N1])
	{
		int s = 1 + j_dest*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1) +
				1*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);

		for (int i = 0 ; i < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; i++)
		{
			for (int k = 0 ; k < boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value ; k++)
			{
				mask[s] = 0;

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
				mask[s] = exist_sub(h,id);

				s += (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1);
			}

			s += 1 * (boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value + 2*1) - boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1);
		}
	}

	template<unsigned int k_dest, unsigned int N1>
	inline static void mask_null(unsigned char mask[N1])
	{
		int s = k_dest + 1*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1) +
				1*(boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type::value+2*1)*(boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value+2*1);

		for (int i = 0 ; i < boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type::value ; i++)
		{
			for (int j = 0 ; j < boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type::value ; j++)
			{
				mask[s] = 0;

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
