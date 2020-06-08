/*
 * SparseGrid_chunk_copy.hpp
 *
 *  Created on: Mar 18, 2020
 *      Author: i-bird
 */

#ifndef SPARSEGRID_CHUNK_COPY_HPP_
#define SPARSEGRID_CHUNK_COPY_HPP_

#ifndef __NVCC__
#include <Vc/Vc>
#endif
#include "util/mathutil.hpp"


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
	return h.mask[sub_id];
}

template<unsigned int v>
struct exist_sub_v_impl
{
	typedef unsigned char type;

	template<typename headerType>
	static inline void exist(headerType & h, int sub_id, unsigned char * pmask)
	{
		pmask[0] = h.mask[sub_id];
	}
};

template<>
struct exist_sub_v_impl<2>
{
	typedef unsigned short type;

	template<typename headerType>
	static inline void exist(headerType & h, int sub_id, unsigned char * pmask)
	{
		pmask[0] = h.mask[sub_id];
		pmask[1] = h.mask[sub_id+1];
	}
};

template<>
struct exist_sub_v_impl<4>
{
	typedef unsigned int type;

	template<typename headerType>
	static inline void exist(headerType & h, int sub_id, unsigned char * pmask)
	{
		pmask[0] = h.mask[sub_id];
		pmask[1] = h.mask[sub_id+1];
		pmask[2] = h.mask[sub_id+2];
		pmask[3] = h.mask[sub_id+3];
	}
};

template<>
struct exist_sub_v_impl<8>
{
	typedef unsigned long int type;

	template<typename headerType>
	static inline void exist(headerType & h, int sub_id, unsigned char * pmask)
	{
		pmask[0] = h.mask[sub_id];
		pmask[1] = h.mask[sub_id+1];
		pmask[2] = h.mask[sub_id+2];
		pmask[3] = h.mask[sub_id+3];
		pmask[4] = h.mask[sub_id+4];
		pmask[5] = h.mask[sub_id+5];
		pmask[6] = h.mask[sub_id+6];
		pmask[7] = h.mask[sub_id+7];
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

#ifndef __NVCC__

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

#endif

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

#ifndef __NVCC__

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

#endif

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

///////////////////////////////////////////////////////////////////////// Chunk missalignment mapping

template<unsigned int dim>
struct key_int
{
	unsigned char i;
	grid_key_dx<dim,long int> k;
};

/*! to move chunks from one Sparse grid to another SparseGrid chunk can be miss-aligned this function create a map between
 *  the missaligned source and the destination chunks
 *
 * For example in an 8x8 chunk in which the chunk in the source sparse grid are shifted by two from the destination sparse-grid
 * you will get a map like this
 *
/verbatim

   2 2 3 3 3 3 3 3
   2 2 3 3 3 3 3 3
   2 2 3 3 3 3 3 3
   2 2 3 3 3 3 3 3
   0 0 1 1 1 1 1 1
   0 0 1 1 1 1 1 1

/endverbatim
 *
 *  and a vk filled with {-1,-1},0 {0,-1},1 {-1,0},2 {0,0},3
 *
 * in case we are comping a slice some of the index in vk can be filtered out. For example if we are copyng a slice in x = 0
 * the vector will only have:
 *
 * {-1,-1},0 {-1,0},2
 *
 */
template<unsigned int dim, unsigned int N, typename chunking>
void construct_chunk_missalign_map(unsigned char miss_al_map[N],
		                           short int mp_off[N],
									 const Box<dim,long int> & b_src,
									 const Box<dim,long int> & b_dst,
									 openfpm::vector<key_int<dim>> & vk)
{
	typedef typename vmpl_create_constant<dim,1>::type one_vmpl;

	get_block_sizes<dim,0,typename chunking::type,one_vmpl> gbs;
	boost::mpl::for_each_ref< boost::mpl::range_c<int,0,dim> >(gbs);

	grid_key_dx<dim> k_src = b_src.getKP1();
	grid_key_dx<dim> k_dst = b_dst.getKP1();

	grid_key_dx<dim,long int> kl1;
	grid_key_dx<dim,long int> kl2;

	// shift the key
	key_shift<dim,chunking>::shift(k_src,kl1);

	// shift the key
	key_shift<dim,chunking>::shift(k_dst,kl2);

	grid_key_dx<dim,long int> key_sub;

	for (int i = 0 ; i < dim ; i++)
	{
		key_sub.set_d(i,kl2.get(i) - kl1.get(i) );
	}

	int c_id[openfpm::math::pow(2,dim)];
	size_t sz[dim];

	for (int i  = 0 ; i < dim ; i++)
	{sz[i] = 3;}

//	grid_sm<dim,void> g2(sz);
//	grid_key_dx_iterator<dim> it2(g2);

/*	grid_key_dx<dim> one;
	one.one();

	while (it2.isNext())
	{
		auto k1 = it2.get();

		g2.LinId(k1);

		key_int<dim> kk;

		kk.i = g2.LinId(k1);
		kk.k = k1 - one;

		++it2;
	}*/

	grid_sm<dim,void> gcnk(gbs.sz_block);
	grid_key_dx_iterator<dim> it(gcnk);
	while (it.isNext())
	{
		auto k = it.get();

		int off = 0;
		int bid = 0;
		int stride = 1;
		int stride_off = 1;
		for (int i = 0 ; i < dim ; i++)
		{
			if ((long int)k.get(i) + key_sub.get(i) >= (long int)gbs.sz_block[i])
			{
				bid += 2*stride;

			}
			else if ((long int)k.get(i) + key_sub.get(i) < 0)
			{
				bid += 0*stride;

			}
			else
			{
				bid += 1*stride;
			}
			off += (openfpm::math::positive_modulo(key_sub.get(i) + k.get(i),gbs.sz_block[i]))*stride_off;
			stride *= 3;
			stride_off *= gbs.sz_block[i];
		}

		miss_al_map[gcnk.LinId(k)] = bid;
		mp_off[gcnk.LinId(k)] = off;

		++it;
	}

	// Filtering: sometime in particular for ghost copy we have a slice to copy, in this case some of the index in the map
	// are never touched so we can eliminate elements in the vk

	Box<dim,long int> slice;

	// calculate slice
	for (int i = 0 ; i < dim ; i++)
	{
		if (kl1.get(i) + (b_src.getHigh(i) - b_src.getLow(i) + 1) > gbs.sz_block[i])
		{
			slice.setLow(i,0);
			slice.setHigh(i,gbs.sz_block[i]-1);
		}
		else
		{
			slice.setLow(i,kl1.get(i));
			slice.setHigh(i,kl1.get(i) + (b_src.getHigh(i) - b_src.getLow(i) + 1));
		}
	}

	key_int<dim> l;
	l.i = 0;
	l.k.zero();
	vk.add(l);

	for (int i  = 0 ; i < dim ; i++)
	{
		if (key_sub.get(i) >= 0)
		{
			size_t bord = gbs.sz_block[i] - key_sub.get(i);

			if (slice.getLow(i) < bord)
			{
				// set the component
				for (int j = 0 ; j < vk.size() ; j++)
				{
					vk.get(j).k.set_d(i,0);
				}
				if (bord <=  slice.getHigh(i) )
				{
					int n_dup = vk.size();

					// add the other parts with zero
					for (int j = 0 ; j < n_dup ; j++)
					{
						key_int<dim> tmp;
						tmp.k = vk.get(j).k;

						tmp.k.set_d(i,1);

						vk.add(tmp);
					}
				}
			}
			else
			{
				// set the component
				for (int j = 0 ; j < vk.size() ; j++)
				{
					vk.get(j).k.set_d(i,0);
				}
			}
		}
		else
		{
			size_t bord = -key_sub.get(i);

			if (slice.getLow(i) < bord)
			{
				// set the component
				for (int j = 0 ; j < vk.size() ; j++)
				{
					vk.get(j).k.set_d(i,-1);
				}
				if (bord <=  slice.getHigh(i) )
				{
					int n_dup = vk.size();

					// add the other parts with zero
					for (int j = 0 ; j < n_dup ; j++)
					{
						key_int<dim> tmp;
						tmp.k = vk.get(j).k;

						tmp.k.set_d(i,0);

						vk.add(tmp);
					}
				}
			}
			else
			{
				// set the component
				for (int j = 0 ; j < vk.size() ; j++)
				{
					vk.get(j).k.set_d(i,0);
				}
			}
		}

/*		size_t bord = (gbs.sz_block[i] - key_sub.get(i));

		if (slice.getLow(i) < (gbs.sz_block[i] - key_sub.get(i)))
		{
			// set the component
			for (int j = 0 ; j < vk.size() ; j++)
			{
				vk.get(j).k.set_d(i,0);
			}
			if ((gbs.sz_block[i] - key_sub.get(i)) <=  slice.getHigh(i) )
			{
				int n_dup = vk.size();

				// add the other parts with zero
				for (int j = 0 ; j < n_dup ; j++)
				{
					key_int<dim> tmp;
					tmp.k = vk.get(j).k;

					tmp.k.set_d(i,0);

					vk.add(tmp);
				}
			}
		}
		else
		{
			// set the component
			for (int j = 0 ; j < vk.size() ; j++)
			{
				vk.get(j).k.set_d(i,0);
			}
		}*/
	}
}


template<typename SparseGridType>
void copy_remove_to_impl(const SparseGridType & grid_src,
							   SparseGridType & grid_dst,
        		  	  	  const Box<SparseGridType::dims,long int> & b_src,
        		  	  	  const Box<SparseGridType::dims,long int> & b_dst)
{
	typedef typename vmpl_create_constant<SparseGridType::dims,1>::type one_vmpl;

	get_block_sizes<SparseGridType::dims,0,typename SparseGridType::chunking_type::type,one_vmpl> gbs;
	boost::mpl::for_each_ref< boost::mpl::range_c<int,0,SparseGridType::dims> >(gbs);

	typedef typename vmpl_reduce_prod<typename SparseGridType::chunking_type::type>::type sizeBlock;

	unsigned char miss_al_map[sizeBlock::value];
	short int mp_off[sizeBlock::value];

	openfpm::vector<key_int<SparseGridType::dims>> vk;

	construct_chunk_missalign_map<SparseGridType::dims,
	                              sizeBlock::value,
	                              typename SparseGridType::chunking_type>
	                             (miss_al_map,mp_off,
								  b_src,b_dst,vk);

	// Now we intersect every chunk source

	Box<SparseGridType::dims,long int> b_inte;
	Box<SparseGridType::dims,long int> b_c;

	auto & data_src = grid_src.private_get_data();
	auto & header_src = grid_src.private_get_header();

	auto & data_dst = grid_src.private_get_data();
	auto & header_dst = grid_src.private_get_header();

	grid_sm<SparseGridType::dims,void> gb;

	int chunk_pos[openfpm::math::pow(2,SparseGridType::dims)];

	for (int i = 0 ; i < header_src.size() ; i++)
	{
		for (int j = 0 ; j < SparseGridType::dims ; j++)
		{
			b_c.setLow(j,header_src.get(i).pos.get(j));
			b_c.setHigh(j,header_src.get(i).pos.get(j) + gbs.sz_block[j] - 1);
		}

		bool inte = b_src.Intersect(b_c,b_inte);

		if (inte == true)
		{
			// Prepare destination chunks

			for (int s = 0 ; s < vk.size() ; s++)
			{
				chunk_pos[vk.get(s).i] = grid_dst.getChunkCreate();
			}

			grid_key_dx_iterator_sub<SparseGridType::dims> it(gb,b_inte.getKP1(),b_inte.getKP2());

			auto & block_src = data_src.get(i);

			while (it.isNext())
			{
				auto id = gb.LinId(it.get());

				int dest_bid = chunk_pos[miss_al_map[id]];

				data_dst.get(dest_bid)[mp_off[id]] = block_src[i];

				++it;
			}
		}
	}
}

#endif /* SPARSEGRID_CHUNK_COPY_HPP_ */
