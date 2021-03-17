/*
 * SparseGrid_conv_opt.hpp
 *
 *  Created on: Jul 19, 2020
 *      Author: i-bird
 */

#ifndef SPARSEGRID_CONV_OPT_HPP_
#define SPARSEGRID_CONV_OPT_HPP_

template<unsigned int l>
union data_il
{
};

template<>
union data_il<8>
{
	typedef long int type;

	unsigned char uc[8];
	long int i;
};

template<>
union data_il<4>
{
	typedef int type;

	unsigned char uc[4];
	int i;
};

template<>
union data_il<2>
{
	typedef short int type;

	unsigned char uc[2];
	short int i;
};

template<>
union data_il<1>
{
	typedef char type;

	unsigned char uc[4];
	char i;
};


template<unsigned int dim, unsigned int sz>
struct ids_crs
{
	long int sumdm[dim];
	long int sumdp[dim];

	long int s2;
	bool mask_row[sz];
	int k;
};



template<unsigned int dim>
struct conv_impl
{
	template<unsigned int prop_src, unsigned int prop_dst, unsigned int stencil_size , unsigned int N, typename SparseGridType, typename lambda_f, typename ... ArgsT >
	void conv(int (& stencil)[N][3], grid_key_dx<3> & start, grid_key_dx<3> & stop, SparseGridType & grid , lambda_f func, ArgsT ... args)
	{
#ifndef __NVCC__
		std::cout << __FILE__ << ":" << __LINE__ << " error conv operation not implemented for this dimension " << std::endl;
#else
		std::cout << __FILE__ << ":" << __LINE__ << " error conv is unsupported when compiled on NVCC " << std::endl;
#endif
	}

	template<bool findNN, unsigned int prop_src, unsigned int prop_dst, unsigned int stencil_size, typename SparseGridType, typename lambda_f, typename ... ArgsT >
	static void conv_cross(grid_key_dx<3> & start, grid_key_dx<3> & stop, SparseGridType & grid , lambda_f func, ArgsT ... args)
	{
#ifndef __NVCC__
		std::cout << __FILE__ << ":" << __LINE__ << " error conv_cross operation not implemented for this dimension " << std::endl;
#else
		std::cout << __FILE__ << ":" << __LINE__ << " error conv_cross is unsupported when compiled on NVCC " << std::endl;
#endif
	}

	template<bool findNN, typename NNType, unsigned int prop_src1, unsigned int prop_src2,
			 unsigned int prop_dst1, unsigned int prop_dst2,
			 unsigned int stencil_size , unsigned int N,
			 typename SparseGridType, typename lambda_f, typename ... ArgsT >
	static void conv2(int (& stencil)[N][3], grid_key_dx<3> & start, grid_key_dx<3> & stop, SparseGridType & grid , lambda_f func, ArgsT ... args)
	{
#ifndef __NVCC__
		std::cout << __FILE__ << ":" << __LINE__ << " error conv2 operation not implemented for this dimension " << std::endl;
#else
		std::cout << __FILE__ << ":" << __LINE__ << " error conv2 is unsupported when compiled on NVCC " << std::endl;
#endif
	}

	template<bool findNN, unsigned int prop_src1, unsigned int prop_src2, unsigned int prop_dst1, unsigned int prop_dst2, unsigned int stencil_size, typename SparseGridType, typename lambda_f, typename ... ArgsT >
	static void conv_cross2(grid_key_dx<3> & start, grid_key_dx<3> & stop, SparseGridType & grid , lambda_f func, ArgsT ... args)
	{
#ifndef __NVCC__
		std::cout << __FILE__ << ":" << __LINE__ << " error conv_cross2 operation not implemented for this dimension " << std::endl;
#else
		std::cout << __FILE__ << ":" << __LINE__ << " error conv_cross2 is unsupported when compiled on NVCC " << std::endl;
#endif
	}
};

#if !defined(__NVCC__) || defined(CUDA_ON_CPU)


template<unsigned int dir,int p, unsigned int prop_src1,typename chunk_type, typename vect_type, typename ids_type>
void load_crs(vect_type & cs1, chunk_type & chunk, ids_type & ids)
{
	if (dir == 0 && p < 0)
	{
		Vc::Vector<typename vect_type::EntryType> cmd1(&chunk.template get<prop_src1>()[ids.s2]);

		cs1 = cmd1;
		cs1 = cs1.shifted(-1);
		cs1[0] = chunk.template get<prop_src1>()[ids.sumdm[dir]];
	}
	else if (dir == 0 && p > 0)
	{
		Vc::Vector<typename vect_type::EntryType> cmd1(&chunk.template get<prop_src1>()[ids.s2]);

		cs1 = cmd1;
		cs1 = cs1.shifted(1);
		cs1[Vc::Vector<typename vect_type::EntryType>::Size - 1] = chunk.template get<prop_src1>()[ids.sumdp[0]];
	}
	else if (p < 0)
	{
		cs1.load(&chunk.template get<prop_src1>()[ids.sumdm[dir]],Vc::Aligned);
	}
	else if (p > 0)
	{
		cs1.load(&chunk.template get<prop_src1>()[ids.sumdp[dir]],Vc::Aligned);
	}
	else
	{
		Vc::Vector<typename vect_type::EntryType> cmd1(&chunk.template get<prop_src1>()[ids.s2]);

		cs1 = cmd1;
	}
}

template<unsigned int prop_dst1,typename chunk_type, typename vect_type, typename ids_type>
void store_crs(chunk_type & chunk, vect_type & res, ids_type & ids)
{
	Vc::Mask<typename vect_type::EntryType> m(&ids.mask_row[ids.k]);

	res.store(&chunk.template get<prop_dst1>()[ids.s2],m,Vc::Aligned);
}

template<unsigned int prop_dst1,unsigned int comp, typename chunk_type, typename vect_type, typename ids_type>
void store_crs_v(chunk_type & chunk, vect_type & res, ids_type & ids)
{
	Vc::Mask<typename vect_type::EntryType> m(&ids.mask_row[ids.k]);

	res.store(&chunk.template get<prop_dst1>()[comp][ids.s2],m,Vc::Aligned);
}

template<unsigned int dir,int p, unsigned int comp, unsigned int prop_src1,typename chunk_type, typename vect_type, typename ids_type>
void load_crs_v(vect_type & cs1, chunk_type & chunk,  ids_type & ids)
{
	if (dir == 0 && p < 0)
	{
		Vc::Vector<typename vect_type::EntryType> cmd1(&chunk.template get<prop_src1>()[comp][ids.s2]);

		cs1 = cmd1;
		cs1 = cs1.shifted(-1);
		cs1[0] = chunk.template get<prop_src1>()[comp][ids.sumdm[dir]];
	}
	else if (dir == 0 && p > 0)
	{
		Vc::Vector<typename vect_type::EntryType> cmd1(&chunk.template get<prop_src1>()[comp][ids.s2]);

		cs1 = cmd1;
		cs1 = cs1.shifted(1);
		cs1[Vc::Vector<typename vect_type::EntryType>::Size - 1] = chunk.template get<prop_src1>()[comp][ids.sumdp[dir]];
	}
	else if (p < 0)
	{
		cs1.load(&chunk.template get<prop_src1>()[comp][ids.sumdm[dir]],Vc::Aligned);
	}
	else if (p > 0)
	{
		cs1.load(&chunk.template get<prop_src1>()[comp][ids.sumdp[dir]],Vc::Aligned);
	}
	else
	{
		Vc::Vector<typename vect_type::EntryType> cmd1(&chunk.template get<prop_src1>()[comp][ids.s2]);

		cs1 = cmd1;
	}
}

struct cross_stencil_v
{
	Vc::double_v xm;
	Vc::double_v xp;
	Vc::double_v ym;
	Vc::double_v yp;
	Vc::double_v zm;
	Vc::double_v zp;
};

template<>
struct conv_impl<3>
{
	template<bool findNN, typename NNtype, unsigned int prop_src, unsigned int prop_dst, unsigned int stencil_size , unsigned int N, typename SparseGridType, typename lambda_f, typename ... ArgsT >
	static void conv(int (& stencil)[N][3], grid_key_dx<3> & start, grid_key_dx<3> & stop, SparseGridType & grid , lambda_f func, ArgsT ... args)
	{
		auto it = grid.template getBlockIterator<stencil_size>(start,stop);

		typedef typename boost::mpl::at<typename SparseGridType::value_type::type, boost::mpl::int_<prop_src>>::type prop_type;

		unsigned char mask[decltype(it)::sizeBlockBord];
		unsigned char mask_sum[decltype(it)::sizeBlockBord];
		unsigned char mask_unused[decltype(it)::sizeBlock];
		__attribute__ ((aligned (32))) prop_type block_bord_src[decltype(it)::sizeBlockBord];
		__attribute__ ((aligned (32))) prop_type block_bord_dst[decltype(it)::sizeBlock];

		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<0>>::type sz0;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<1>>::type sz1;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<2>>::type sz2;

		while (it.isNext())
		{
			it.template loadBlockBorder<prop_src,NNtype,findNN>(block_bord_src,mask);

			if (it.start_b(2) != stencil_size || it.start_b(1) != stencil_size || it.start_b(0) != stencil_size ||
			    it.stop_b(2) != sz2::value+stencil_size || it.stop_b(1) != sz1::value+stencil_size || it.stop_b(0) != sz0::value+stencil_size)
			{
				auto & header_mask = grid.private_get_header_mask();
				auto & header_inf = grid.private_get_header_inf();

				loadBlock_impl<prop_dst,0,3,typename decltype(it)::vector_blocks_exts_type, typename decltype(it)::vector_ext_type>::template loadBlock<decltype(it)::sizeBlock>(block_bord_dst,grid,it.getChunkId(),mask_unused);
			}

			// Sum the mask
			for (int k = it.start_b(2) ; k < it.stop_b(2) ; k++)
			{
				for (int j = it.start_b(1) ; j < it.stop_b(1) ; j++)
				{
					int cc = it.LinB(it.start_b(0),j,k);
					int c[N];

					for (int s = 0 ; s < N ; s++)
					{
						c[s] = it.LinB(it.start_b(0)+stencil[s][0],j+stencil[s][1],k+stencil[s][2]);
					}

					for (int i = it.start_b(0) ; i < it.stop_b(0) ; i += sizeof(size_t))
					{
						size_t cmd = *(size_t *)&mask[cc];

						if (cmd != 0)
						{
							size_t xm[N];

							for (int s = 0 ; s < N ; s++)
							{
								xm[s] = *(size_t *)&mask[c[s]];
							}

							size_t sum = 0;
							for (int s = 0 ; s < N ; s++)
							{
								sum += xm[s];
							}

							*(size_t *)&mask_sum[cc] = sum;
						}

						cc += sizeof(size_t);
						for (int s = 0 ; s < N ; s++)
						{
							c[s] += sizeof(size_t);
						}
					}
				}
			}

			for (int k = it.start_b(2) ; k < it.stop_b(2) ; k++)
			{
				for (int j = it.start_b(1) ; j < it.stop_b(1) ; j++)
				{
					int cc = it.LinB(it.start_b(0),j,k);
					int c[N];

					int cd = it.LinB_off(it.start_b(0),j,k);

					for (int s = 0 ; s < N ; s++)
					{
						c[s] = it.LinB(it.start_b(0)+stencil[s][0],j+stencil[s][1],k+stencil[s][2]);
					}

					for (int i = it.start_b(0) ; i < it.stop_b(0) ; i += Vc::Vector<prop_type>::Size)
					{
						Vc::Mask<prop_type> cmp;

						for (int s = 0 ; s < Vc::Vector<prop_type>::Size ; s++)
						{
							cmp[s] = (mask[cc+s] == true && i+s < it.stop_b(0));
						}

						// we do only if exist the point
						if (Vc::none_of(cmp) == false)
						{
							Vc::Mask<prop_type> surround;

							Vc::Vector<prop_type> xs[N+1];

							xs[0] = Vc::Vector<prop_type>(&block_bord_src[cc],Vc::Unaligned);

							for (int s = 1 ; s < N+1 ; s++)
							{
								xs[s] = Vc::Vector<prop_type>(&block_bord_src[c[s-1]],Vc::Unaligned);
							}

							auto res = func(xs, &mask_sum[cc], args ...);

							res.store(&block_bord_dst[cd],cmp,Vc::Aligned);
						}

						cc += Vc::Vector<prop_type>::Size;
						for (int s = 0 ; s < N ; s++)
						{
							c[s] += Vc::Vector<prop_type>::Size;
						}
						cd += Vc::Vector<prop_type>::Size;
					}
				}
			}

			it.template storeBlock<prop_dst>(block_bord_dst);

			++it;
		}
	}

	template<bool findNN, unsigned int prop_src, unsigned int prop_dst, unsigned int stencil_size, typename SparseGridType, typename lambda_f, typename ... ArgsT >
	static void conv_cross(grid_key_dx<3> & start, grid_key_dx<3> & stop, SparseGridType & grid , lambda_f func, ArgsT ... args)
	{
		auto it = grid.template getBlockIterator<1>(start,stop);

		auto & datas = grid.private_get_data();
		auto & headers = grid.private_get_header_mask();

		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<0>>::type sz0;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<1>>::type sz1;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<2>>::type sz2;

		typedef typename SparseGridType::chunking_type chunking;

		typedef typename boost::mpl::at<typename SparseGridType::value_type::type, boost::mpl::int_<prop_src>>::type prop_type;

		while (it.isNext())
		{
			// Load
			long int offset_jump[6];

			size_t cid = it.getChunkId();

			auto chunk = datas.get(cid);
			auto & mask = headers.get(cid);

			bool exist;
			grid_key_dx<3> p = grid.getChunkPos(cid) + grid_key_dx<3>({-1,0,0});
			long int r = grid.getChunk(p,exist);
			offset_jump[0] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({1,0,0});
			r = grid.getChunk(p,exist);
			offset_jump[1] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,-1,0});
			r = grid.getChunk(p,exist);
			offset_jump[2] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,1,0});
			r = grid.getChunk(p,exist);
			offset_jump[3] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,0,-1});
			r = grid.getChunk(p,exist);
			offset_jump[4] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,0,1});
			r = grid.getChunk(p,exist);
			offset_jump[5] = (r-cid)*decltype(it)::sizeBlock;

			// Load offset jumps

			// construct a row mask

			long int s2 = 0;

			typedef typename boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type sz;
			typedef typename boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type sy;
			typedef typename boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type sx;


			bool mask_row[sx::value];

			for (int k = 0 ; k < sx::value ; k++)
			{
				mask_row[k] = (k >= it.start(0) && k < it.stop(0))?true:false;
			}

			for (int v = it.start(2) ; v < it.stop(2) ; v++)
			{
				for (int j = it.start(1) ; j < it.stop(1) ; j++)
				{
					s2 = it.Lin(0,j,v);
					for (int k = 0 ; k < sx::value ; k += Vc::Vector<prop_type>::Size)
					{
						// we do only id exist the point
						if (*(int *)&mask.mask[s2] == 0) {s2 += Vc::Vector<prop_type>::Size; continue;}

						data_il<Vc::Vector<prop_type>::Size> mxm;
						data_il<Vc::Vector<prop_type>::Size> mxp;
						data_il<Vc::Vector<prop_type>::Size> mym;
						data_il<Vc::Vector<prop_type>::Size> myp;
						data_il<Vc::Vector<prop_type>::Size> mzm;
						data_il<Vc::Vector<prop_type>::Size> mzp;

						cross_stencil_v cs;

						Vc::Vector<prop_type> cmd(&chunk.template get<prop_src>()[s2]);

						// Load x-1
						long int sumxm = s2-1;
						sumxm += (k==0)?offset_jump[0] + sx::value:0;

						// Load x+1
						long int sumxp = s2+Vc::Vector<prop_type>::Size;
						sumxp += (k+Vc::Vector<prop_type>::Size == sx::value)?offset_jump[1] - sx::value:0;

						long int sumym = (j == 0)?offset_jump[2] + (sy::value-1)*sx::value:-sx::value;
						sumym += s2;
						long int sumyp = (j == sy::value-1)?offset_jump[3] - (sy::value - 1)*sx::value:sx::value;
						sumyp += s2;
						long int sumzm = (v == 0)?offset_jump[4] + (sz::value-1)*sx::value*sy::value:-sx::value*sy::value;
						sumzm += s2;
						long int sumzp = (v == sz::value-1)?offset_jump[5] - (sz::value - 1)*sx::value*sy::value:sx::value*sy::value;
						sumzp += s2;

						if (Vc::Vector<prop_type>::Size == 2)
						{
							mxm.i = *(short int *)&mask.mask[s2];
							mxm.i = mxm.i << 8;
							mxm.i |= (short int)mask.mask[sumxm];

							mxp.i = *(short int *)&mask.mask[s2];
							mxp.i = mxp.i >> 8;
							mxp.i |= ((short int)mask.mask[sumxp]) << (Vc::Vector<prop_type>::Size - 1)*8;

							mym.i = *(short int *)&mask.mask[sumym];
							myp.i = *(short int *)&mask.mask[sumyp];

							mzm.i = *(short int *)&mask.mask[sumzm];
							mzp.i = *(short int *)&mask.mask[sumzp];
						}
						else if (Vc::Vector<prop_type>::Size == 4)
						{
							mxm.i = *(int *)&mask.mask[s2];
							mxm.i = mxm.i << 8;
							mxm.i |= (int)mask.mask[sumxm];

							mxp.i = *(int *)&mask.mask[s2];
							mxp.i = mxp.i >> 8;
							mxp.i |= ((int)mask.mask[sumxp]) << (Vc::Vector<prop_type>::Size - 1)*8;

							mym.i = *(int *)&mask.mask[sumym];
							myp.i = *(int *)&mask.mask[sumyp];

							mzm.i = *(int *)&mask.mask[sumzm];
							mzp.i = *(int *)&mask.mask[sumzp];
						}
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " UNSUPPORTED" << std::endl;
						}

						cs.xm = cmd;
						cs.xm = cs.xm.shifted(-1);
						cs.xm[0] = chunk.template get<prop_src>()[sumxm];


						cs.xp = cmd;
						cs.xp = cs.xp.shifted(1);
						cs.xp[Vc::Vector<prop_type>::Size - 1] = chunk.template get<prop_src>()[sumxp];

						// Load y and z direction

						cs.ym.load(&chunk.template get<prop_src>()[sumym],Vc::Aligned);
						cs.yp.load(&chunk.template get<prop_src>()[sumyp],Vc::Aligned);
						cs.zm.load(&chunk.template get<prop_src>()[sumzm],Vc::Aligned);
						cs.zp.load(&chunk.template get<prop_src>()[sumzp],Vc::Aligned);

						// Calculate

						data_il<Vc::Vector<prop_type>::Size> tot_m;
						tot_m.i = mxm.i + mxp.i + mym.i + myp.i + mzm.i + mzp.i;

						Vc::Vector<prop_type> res = func(cmd,cs,tot_m.uc,args ... );

						Vc::Mask<prop_type> m(&mask_row[k]);

						res.store(&chunk.template get<prop_dst>()[s2],m,Vc::Aligned);

						s2 += Vc::Vector<prop_type>::Size;
					}
				}
			}

			++it;
		}
	}


	template<bool findNN, typename NNType, unsigned int prop_src1, unsigned int prop_src2,
			 unsigned int prop_dst1, unsigned int prop_dst2,
			 unsigned int stencil_size , unsigned int N,
			 typename SparseGridType, typename lambda_f, typename ... ArgsT >
	static void conv2(int (& stencil)[N][3], grid_key_dx<3> & start, grid_key_dx<3> & stop, SparseGridType & grid , lambda_f func, ArgsT ... args)
	{
		auto it = grid.template getBlockIterator<stencil_size>(start,stop);

		typedef typename boost::mpl::at<typename SparseGridType::value_type::type, boost::mpl::int_<prop_src1>>::type prop_type;

		unsigned char mask[decltype(it)::sizeBlockBord];
		unsigned char mask_sum[decltype(it)::sizeBlockBord];
		__attribute__ ((aligned (64))) prop_type block_bord_src1[decltype(it)::sizeBlockBord];
		__attribute__ ((aligned (64))) prop_type block_bord_dst1[decltype(it)::sizeBlock+16];
		__attribute__ ((aligned (64))) prop_type block_bord_src2[decltype(it)::sizeBlockBord];
		__attribute__ ((aligned (64))) prop_type block_bord_dst2[decltype(it)::sizeBlock+16];

		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<0>>::type sz0;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<1>>::type sz1;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<2>>::type sz2;

		while (it.isNext())
		{
			it.template loadBlockBorder<prop_src1,NNType,findNN>(block_bord_src1,mask);
			it.template loadBlockBorder<prop_src2,NNType,findNN>(block_bord_src2,mask);

			// Sum the mask
			for (int k = it.start_b(2) ; k < it.stop_b(2) ; k++)
			{
				for (int j = it.start_b(1) ; j < it.stop_b(1) ; j++)
				{
					int cc = it.LinB(it.start_b(0),j,k);
					int c[N];

					for (int s = 0 ; s < N ; s++)
					{
						c[s] = it.LinB(it.start_b(0)+stencil[s][0],j+stencil[s][1],k+stencil[s][2]);
					}

					for (int i = it.start_b(0) ; i < it.stop_b(0) ; i += sizeof(size_t))
					{
						size_t cmd = *(size_t *)&mask[cc];

						if (cmd == 0) {continue;}


						size_t xm[N];

						for (int s = 0 ; s < N ; s++)
						{
							xm[s] = *(size_t *)&mask[c[s]];
						}

						size_t sum = 0;
						for (int s = 0 ; s < N ; s++)
						{
							sum += xm[s];
						}

						*(size_t *)&mask_sum[cc] = sum;

						cc += sizeof(size_t);
						for (int s = 0 ; s < N ; s++)
						{
							c[s] += sizeof(size_t);
						}
					}
				}
			}

			for (int k = it.start_b(2) ; k < it.stop_b(2) ; k++)
			{
				for (int j = it.start_b(1) ; j < it.stop_b(1) ; j++)
				{
					int cc = it.LinB(it.start_b(0),j,k);
					int c[N];

					int cd = it.LinB_off(it.start_b(0),j,k);

					for (int s = 0 ; s < N ; s++)
					{
						c[s] = it.LinB(it.start_b(0)+stencil[s][0],j+stencil[s][1],k+stencil[s][2]);
					}

					for (int i = it.start_b(0) ; i < it.stop_b(0) ; i += Vc::Vector<prop_type>::Size)
					{
						Vc::Mask<prop_type> cmp;

						for (int s = 0 ; s < Vc::Vector<prop_type>::Size ; s++)
						{
							cmp[s] = (mask[cc+s] == true);
						}

						// we do only id exist the point
						if (Vc::none_of(cmp) == true) {continue;}

						Vc::Mask<prop_type> surround;

						Vc::Vector<prop_type> xs1[N+1];
						Vc::Vector<prop_type> xs2[N+1];

						xs1[0] = Vc::Vector<prop_type>(&block_bord_src1[cc],Vc::Unaligned);
						xs2[0] = Vc::Vector<prop_type>(&block_bord_src2[cc],Vc::Unaligned);

						for (int s = 1 ; s < N+1 ; s++)
						{
							xs1[s] = Vc::Vector<prop_type>(&block_bord_src1[c[s-1]],Vc::Unaligned);
							xs2[s] = Vc::Vector<prop_type>(&block_bord_src2[c[s-1]],Vc::Unaligned);
						}

						Vc::Vector<prop_type> vo1;
						Vc::Vector<prop_type> vo2;

						func(vo1, vo2, xs1, xs2, &mask_sum[cc], args ...);

						vo1.store(&block_bord_dst1[cd],cmp,Vc::Unaligned);
						vo2.store(&block_bord_dst2[cd],cmp,Vc::Unaligned);

						cc += Vc::Vector<prop_type>::Size;
						for (int s = 0 ; s < N ; s++)
						{
							c[s] += Vc::Vector<prop_type>::Size;
						}
						cd += Vc::Vector<prop_type>::Size;
					}
				}
			}

			it.template storeBlock<prop_dst1>(block_bord_dst1);
			it.template storeBlock<prop_dst2>(block_bord_dst2);

			++it;
		}
	}

	template<bool findNN, unsigned int prop_src1, unsigned int prop_src2, unsigned int prop_dst1, unsigned int prop_dst2, unsigned int stencil_size, typename SparseGridType, typename lambda_f, typename ... ArgsT >
	static void conv_cross2(grid_key_dx<3> & start, grid_key_dx<3> & stop, SparseGridType & grid , lambda_f func, ArgsT ... args)
	{
		auto it = grid.template getBlockIterator<stencil_size>(start,stop);

		auto & datas = grid.private_get_data();
		auto & headers = grid.private_get_header_mask();

		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<0>>::type sz0;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<1>>::type sz1;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<2>>::type sz2;

		typedef typename SparseGridType::chunking_type chunking;

		typedef typename boost::mpl::at<typename SparseGridType::value_type::type, boost::mpl::int_<prop_src1>>::type prop_type;

		while (it.isNext())
		{
			// Load
			long int offset_jump[6];

			size_t cid = it.getChunkId();

			auto chunk = datas.get(cid);
			auto & mask = headers.get(cid);

			bool exist;
			grid_key_dx<3> p = grid.getChunkPos(cid) + grid_key_dx<3>({-1,0,0});
			long int r = grid.getChunk(p,exist);
			offset_jump[0] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({1,0,0});
			r = grid.getChunk(p,exist);
			offset_jump[1] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,-1,0});
			r = grid.getChunk(p,exist);
			offset_jump[2] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,1,0});
			r = grid.getChunk(p,exist);
			offset_jump[3] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,0,-1});
			r = grid.getChunk(p,exist);
			offset_jump[4] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,0,1});
			r = grid.getChunk(p,exist);
			offset_jump[5] = (r-cid)*decltype(it)::sizeBlock;

			// Load offset jumps

			// construct a row mask

			long int s2 = 0;

			typedef typename boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type sz;
			typedef typename boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type sy;
			typedef typename boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type sx;


			bool mask_row[sx::value];

			for (int k = 0 ; k < sx::value ; k++)
			{
				mask_row[k] = (k >= it.start(0) && k < it.stop(0))?true:false;
			}

			for (int v = it.start(2) ; v < it.stop(2) ; v++)
			{
				for (int j = it.start(1) ; j < it.stop(1) ; j++)
				{
					s2 = it.Lin(0,j,v);
					for (int k = 0 ; k < sx::value ; k += Vc::Vector<prop_type>::Size)
					{
						// we do only id exist the point
						if (*(int *)&mask.mask[s2] == 0) {s2 += Vc::Vector<prop_type>::Size; continue;}

						data_il<4> mxm;
						data_il<4> mxp;
						data_il<4> mym;
						data_il<4> myp;
						data_il<4> mzm;
						data_il<4> mzp;

						cross_stencil_v cs1;
						cross_stencil_v cs2;

						Vc::Vector<prop_type> cmd1(&chunk.template get<prop_src1>()[s2]);
						Vc::Vector<prop_type> cmd2(&chunk.template get<prop_src2>()[s2]);

						// Load x-1
						long int sumxm = s2-1;
						sumxm += (k==0)?offset_jump[0] + sx::value:0;

						// Load x+1
						long int sumxp = s2+Vc::Vector<prop_type>::Size;
						sumxp += (k+Vc::Vector<prop_type>::Size == sx::value)?offset_jump[1] - sx::value:0;

						long int sumym = (j == 0)?offset_jump[2] + (sy::value-1)*sx::value:-sx::value;
						sumym += s2;
						long int sumyp = (j == sy::value-1)?offset_jump[3] - (sy::value - 1)*sx::value:sx::value;
						sumyp += s2;
						long int sumzm = (v == 0)?offset_jump[4] + (sz::value-1)*sx::value*sy::value:-sx::value*sy::value;
						sumzm += s2;
						long int sumzp = (v == sz::value-1)?offset_jump[5] - (sz::value - 1)*sx::value*sy::value:sx::value*sy::value;
						sumzp += s2;

						if (Vc::Vector<prop_type>::Size == 2)
						{
							mxm.i = *(short int *)&mask.mask[s2];
							mxm.i = mxm.i << 8;
							mxm.i |= (short int)mask.mask[sumxm];

							mxp.i = *(short int *)&mask.mask[s2];
							mxp.i = mxp.i >> 8;
							mxp.i |= ((short int)mask.mask[sumxp]) << (Vc::Vector<prop_type>::Size - 1)*8;

							mym.i = *(short int *)&mask.mask[sumym];
							myp.i = *(short int *)&mask.mask[sumyp];

							mzm.i = *(short int *)&mask.mask[sumzm];
							mzp.i = *(short int *)&mask.mask[sumzp];
						}
						else if (Vc::Vector<prop_type>::Size == 4)
						{
							mxm.i = *(int *)&mask.mask[s2];
							mxm.i = mxm.i << 8;
							mxm.i |= (int)mask.mask[sumxm];

							mxp.i = *(int *)&mask.mask[s2];
							mxp.i = mxp.i >> 8;
							mxp.i |= ((int)mask.mask[sumxp]) << (Vc::Vector<prop_type>::Size - 1)*8;

							mym.i = *(int *)&mask.mask[sumym];
							myp.i = *(int *)&mask.mask[sumyp];

							mzm.i = *(int *)&mask.mask[sumzm];
							mzp.i = *(int *)&mask.mask[sumzp];
						}
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " UNSUPPORTED" << std::endl;
						}

						cs1.xm = cmd1;
						cs1.xm = cs1.xm.shifted(-1);
						cs1.xm[0] = chunk.template get<prop_src1>()[sumxm];

						cs2.xm = cmd2;
						cs2.xm = cs2.xm.shifted(-1);
						cs2.xm[0] = chunk.template get<prop_src2>()[sumxm];

						cs1.xp = cmd1;
						cs1.xp = cs1.xp.shifted(1);
						cs1.xp[Vc::Vector<prop_type>::Size - 1] = chunk.template get<prop_src1>()[sumxp];

						cs2.xp = cmd2;
						cs2.xp = cs2.xp.shifted(1);
						cs2.xp[Vc::Vector<prop_type>::Size - 1] = chunk.template get<prop_src2>()[sumxp];

						// Load y and z direction

						cs1.ym.load(&chunk.template get<prop_src1>()[sumym],Vc::Aligned);
						cs1.yp.load(&chunk.template get<prop_src1>()[sumyp],Vc::Aligned);
						cs1.zm.load(&chunk.template get<prop_src1>()[sumzm],Vc::Aligned);
						cs1.zp.load(&chunk.template get<prop_src1>()[sumzp],Vc::Aligned);

						cs2.ym.load(&chunk.template get<prop_src2>()[sumym],Vc::Aligned);
						cs2.yp.load(&chunk.template get<prop_src2>()[sumyp],Vc::Aligned);
						cs2.zm.load(&chunk.template get<prop_src2>()[sumzm],Vc::Aligned);
						cs2.zp.load(&chunk.template get<prop_src2>()[sumzp],Vc::Aligned);

						// Calculate

						data_il<4> tot_m;
						tot_m.i = mxm.i + mxp.i + mym.i + myp.i + mzm.i + mzp.i;

						Vc::Vector<prop_type> res1;
						Vc::Vector<prop_type> res2;

						func(res1,res2,cmd1,cmd2,cs1,cs2,tot_m.uc,args ... );

						Vc::Mask<prop_type> m(&mask_row[k]);

						res1.store(&chunk.template get<prop_dst1>()[s2],m,Vc::Aligned);
						res2.store(&chunk.template get<prop_dst2>()[s2],m,Vc::Aligned);

						s2 += Vc::Vector<prop_type>::Size;
					}
				}
			}

			++it;
		}
	}

	template<bool findNN, unsigned int stencil_size, typename prop_type, typename SparseGridType, typename lambda_f, typename ... ArgsT >
	static void conv_cross_ids(grid_key_dx<3> & start, grid_key_dx<3> & stop, SparseGridType & grid , lambda_f func, ArgsT ... args)
	{
		auto it = grid.template getBlockIterator<stencil_size>(start,stop);

		auto & datas = grid.private_get_data();
		auto & headers = grid.private_get_header_mask();

		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<0>>::type sz0;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<1>>::type sz1;
		typedef typename boost::mpl::at<typename decltype(it)::stop_border_vmpl,boost::mpl::int_<2>>::type sz2;

		typedef typename SparseGridType::chunking_type chunking;

		while (it.isNext())
		{
			// Load
			long int offset_jump[6];

			size_t cid = it.getChunkId();

			auto chunk = datas.get(cid);
			auto & mask = headers.get(cid);

			bool exist;
			grid_key_dx<3> p = grid.getChunkPos(cid) + grid_key_dx<3>({-1,0,0});
			long int r = grid.getChunk(p,exist);
			offset_jump[0] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({1,0,0});
			r = grid.getChunk(p,exist);
			offset_jump[1] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,-1,0});
			r = grid.getChunk(p,exist);
			offset_jump[2] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,1,0});
			r = grid.getChunk(p,exist);
			offset_jump[3] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,0,-1});
			r = grid.getChunk(p,exist);
			offset_jump[4] = (r-cid)*decltype(it)::sizeBlock;

			p = grid.getChunkPos(cid) + grid_key_dx<3>({0,0,1});
			r = grid.getChunk(p,exist);
			offset_jump[5] = (r-cid)*decltype(it)::sizeBlock;

			// Load offset jumps

			// construct a row mask

			long int s2 = 0;

			typedef typename boost::mpl::at<typename chunking::type,boost::mpl::int_<2>>::type sz;
			typedef typename boost::mpl::at<typename chunking::type,boost::mpl::int_<1>>::type sy;
			typedef typename boost::mpl::at<typename chunking::type,boost::mpl::int_<0>>::type sx;

			ids_crs<3,sx::value> ids;

			for (int k = 0 ; k < sx::value ; k++)
			{
				ids.mask_row[k] = (k >= it.start(0) && k < it.stop(0))?true:false;
			}

			for (int v = it.start(2) ; v < it.stop(2) ; v++)
			{
				for (int j = it.start(1) ; j < it.stop(1) ; j++)
				{
					s2 = it.Lin(0,j,v);
					for (int k = 0 ; k < sx::value ; k += Vc::Vector<prop_type>::Size)
					{
						// we do only id exist the point
						if (*(int *)&mask.mask[s2] == 0) {s2 += Vc::Vector<prop_type>::Size; continue;}

						data_il<4> mxm;
						data_il<4> mxp;
						data_il<4> mym;
						data_il<4> myp;
						data_il<4> mzm;
						data_il<4> mzp;

						ids.k = k;

						// Load x-1
						ids.sumdm[0] = s2-1;
						ids.sumdm[0] += (k==0)?offset_jump[0] + sx::value:0;

						// Load x+1
						ids.sumdp[0] = s2+Vc::Vector<prop_type>::Size;
						ids.sumdp[0] += (k+Vc::Vector<prop_type>::Size == sx::value)?offset_jump[1] - sx::value:0;

						ids.sumdm[1] = (j == 0)?offset_jump[2] + (sy::value-1)*sx::value:-sx::value;
						ids.sumdm[1] += s2;
						ids.sumdp[1] = (j == sy::value-1)?offset_jump[3] - (sy::value - 1)*sx::value:sx::value;
						ids.sumdp[1] += s2;
						ids.sumdm[2] = (v == 0)?offset_jump[4] + (sz::value-1)*sx::value*sy::value:-sx::value*sy::value;
						ids.sumdm[2] += s2;
						ids.sumdp[2] = (v == sz::value-1)?offset_jump[5] - (sz::value - 1)*sx::value*sy::value:sx::value*sy::value;
						ids.sumdp[2] += s2;

						ids.s2 = s2;

                        if (Vc::Vector<prop_type>::Size == 2)
                        {
                            mxm.i = *(short int *)&mask.mask[s2];
                            mxm.i = mxm.i << 8;
                            mxm.i |= (short int)mask.mask[ids.sumdm[0]];

                            mxp.i = *(short int *)&mask.mask[s2];
                            mxp.i = mxp.i >> 8;
                            mxp.i |= ((short int)mask.mask[ids.sumdp[0]]) << (Vc::Vector<prop_type>::Size - 1)*8;

                            mym.i = *(short int *)&mask.mask[ids.sumdm[1]];
                            myp.i = *(short int *)&mask.mask[ids.sumdp[1]];

                            mzm.i = *(short int *)&mask.mask[ids.sumdm[2]];
                            mzp.i = *(short int *)&mask.mask[ids.sumdp[2]];
                        }
                        else if (Vc::Vector<prop_type>::Size == 4)
                        {
                            mxm.i = *(int *)&mask.mask[s2];
                            mxm.i = mxm.i << 8;
                            mxm.i |= (int)mask.mask[ids.sumdm[0]];

                            mxp.i = *(int *)&mask.mask[s2];
                            mxp.i = mxp.i >> 8;
                            mxp.i |= ((int)mask.mask[ids.sumdp[0]]) << (Vc::Vector<prop_type>::Size - 1)*8;

                        	mym.i = *(int *)&mask.mask[ids.sumdm[1]];
                            myp.i = *(int *)&mask.mask[ids.sumdp[1]];

                        	mzm.i = *(int *)&mask.mask[ids.sumdm[2]];
                            mzp.i = *(int *)&mask.mask[ids.sumdp[2]];
                        }
                        else
                        {
                            std::cout << __FILE__ << ":" << __LINE__ << " UNSUPPORTED" << std::endl;
                        }

						// Calculate

						data_il<4> tot_m;
						tot_m.i = mxm.i + mxp.i + mym.i + myp.i + mzm.i + mzp.i;

						func(chunk,ids,tot_m.uc,args ... );

						s2 += Vc::Vector<prop_type>::Size;
					}
				}
			}

			++it;
		}
	}

};

#endif


#endif /* SPARSEGRID_CONV_OPT_HPP_ */
