#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <Vector/map_vector.hpp>
#include <Vector/cuda/map_vector_sparse_cuda_kernels.cuh>
#include "util/cuda/moderngpu/kernel_scan.hxx"
#include "util/cuda/moderngpu/kernel_merge.hxx"

BOOST_AUTO_TEST_SUITE( vector_sparse_cuda_kernels )

BOOST_AUTO_TEST_CASE(test_maps_create)
{
	openfpm::vector_gpu<aggregate<int>> merge_indexes;

	openfpm::vector_gpu<aggregate<int,int,int,int>> p_ids;
	openfpm::vector_gpu<aggregate<int,int,int,int>> s_ids;

	openfpm::vector_gpu<aggregate<int>> copy_old;
	openfpm::vector_gpu<aggregate<int>> outputMap;
	openfpm::vector_gpu<aggregate<int>> segments_oldData;

	// old data has 9 points
	// new data has 7 segments

	merge_indexes.resize(16);
	merge_indexes.template get<0>(0) = 9;      // new
	merge_indexes.template get<0>(1) = 10;     // new
	merge_indexes.template get<0>(2) = 0;      // old id:0
	merge_indexes.template get<0>(3) = 1;      // old id:1
	merge_indexes.template get<0>(4) = 2;      // old id:2
	merge_indexes.template get<0>(5) = 3;      // old id:3
	merge_indexes.template get<0>(6) = 4;      // old id:4
	merge_indexes.template get<0>(7) = 5;      // old id:5
	merge_indexes.template get<0>(8) = 11;     // new
	merge_indexes.template get<0>(9) = 6;      // old id:6
	merge_indexes.template get<0>(10) = 12;     // new
	merge_indexes.template get<0>(11) = 7;     // old id:7
	merge_indexes.template get<0>(12) = 8;     // old id: 8
	merge_indexes.template get<0>(13) = 13;    // new
	merge_indexes.template get<0>(14) = 14;    // new
	merge_indexes.template get<0>(15) = 15;    // new

	merge_indexes.template hostToDevice<0>();

	s_ids.resize(16+1);
	p_ids.resize(16+1);

	// fill the last of compute predicates
	p_ids.template get<0>(p_ids.size()-1) = 0;
	p_ids.template get<1>(p_ids.size()-1) = 0;
	p_ids.template get<2>(p_ids.size()-1) = 0;
	p_ids.template get<3>(p_ids.size()-1) = 0;

	p_ids.template hostToDevice<0,1,2,3>(p_ids.size()-1,p_ids.size()-1);

	auto ite = merge_indexes.getGPUIterator();

	CUDA_LAUNCH(compute_predicate,ite,merge_indexes.toKernel(),9,p_ids.toKernel());

	mgpu::standard_context_t context(false);
	mgpu::scan((int *)p_ids.template getDeviceBuffer<0>(),
				s_ids.size(),
	            (int *)s_ids.template getDeviceBuffer<0>(),
                context);

	mgpu::scan((int *)p_ids.template getDeviceBuffer<1>(),
				s_ids.size(),
	            (int *)s_ids.template getDeviceBuffer<1>(),
                context);

	mgpu::scan((int *)p_ids.template getDeviceBuffer<2>(),
				s_ids.size(),
	            (int *)s_ids.template getDeviceBuffer<2>(),
                context);

	mgpu::scan((int *)p_ids.template getDeviceBuffer<3>(),
				s_ids.size(),
	            (int *)s_ids.template getDeviceBuffer<3>(),
                context);

	s_ids.template deviceToHost<0,1,2,3>();
	p_ids.template deviceToHost<0,1,2,3>();

	size_t copy_old_size = s_ids.template get<3>(s_ids.size()-1) + p_ids.template get<3>(p_ids.size()-1);
	size_t seg_old_size = s_ids.template get<1>(s_ids.size()-1) + p_ids.template get<1>(p_ids.size()-1);
	size_t out_map_size = s_ids.template get<1>(s_ids.size()-1) + p_ids.template get<1>(p_ids.size()-1);

	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(0),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(0),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(0),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(0),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(1),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(1),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(1),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(1),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(2),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(2),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(2),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(2),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(3),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(3),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(3),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(3),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(4),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(4),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(4),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(4),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(5),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(5),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(5),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(5),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(6),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(6),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(6),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(6),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(7),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(7),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(7),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(7),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(8),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(8),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(8),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(8),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(9),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(9),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(9),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(9),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(10),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(10),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(10),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(10),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(11),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(11),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(11),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(11),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(12),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(12),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(12),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(12),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(13),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(13),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(13),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(13),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(14),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(14),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(14),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(14),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<0>(15),0);
	BOOST_REQUIRE_EQUAL(p_ids.template get<1>(15),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<2>(15),1);
	BOOST_REQUIRE_EQUAL(p_ids.template get<3>(15),0);

	segments_oldData.resize(seg_old_size);
	outputMap.resize(out_map_size);
	copy_old.resize(copy_old_size);

	CUDA_LAUNCH(maps_create,ite,s_ids.toKernel(),p_ids.toKernel(),segments_oldData.toKernel(),outputMap.toKernel(),copy_old.toKernel());

	segments_oldData.template deviceToHost<0>();
	outputMap.template deviceToHost<0>();
	copy_old.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(seg_old_size,7);
	BOOST_REQUIRE_EQUAL(out_map_size,7);
	BOOST_REQUIRE_EQUAL(copy_old_size,6);

	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(0),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(0),0);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(1),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(1),1);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(2),5);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(2),7);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(3),6);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(3),8);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(4),8);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(4),10);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(5),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(5),11);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(6),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(6),12);

	BOOST_REQUIRE_EQUAL(copy_old.template get<0>(0),2);
	BOOST_REQUIRE_EQUAL(copy_old.template get<0>(1),3);
	BOOST_REQUIRE_EQUAL(copy_old.template get<0>(2),4);
	BOOST_REQUIRE_EQUAL(copy_old.template get<0>(3),5);
	BOOST_REQUIRE_EQUAL(copy_old.template get<0>(4),6);
	BOOST_REQUIRE_EQUAL(copy_old.template get<0>(5),9);
}

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_kernels_use )
{
	openfpm::vector_gpu<aggregate<int>> block_insert;
	openfpm::vector_gpu<aggregate<int>> block_n;
	openfpm::vector_gpu<aggregate<int>> block_n_scan;
	openfpm::vector_gpu<aggregate<int>> output_list_0;
	openfpm::vector_gpu<aggregate<int>> output_list_1;

	openfpm::vector_gpu<aggregate<int,float,double>> block_data;
	openfpm::vector_gpu<aggregate<int,float,double>> cont_data;

	int nblock = 100;
	int nslot = 512;
	block_insert.resize(nblock*nslot);
	block_data.resize(nblock*nslot);
	block_n.resize(nblock+1);
	block_n_scan.resize(nblock+1);

	// fill block insert of some data
	for (int i = 0 ; i < nblock ; i++)
	{
		block_n.template get<0>(i) = ((float)rand() / RAND_MAX) * 511;
		for (int j = 0 ; j < block_n.template get<0>(i) ; j++)
		{
			block_insert.template get<0>(i*nslot + j) = ((float)rand() / RAND_MAX) * 511;

			block_data.template get<0>(i*nslot + j) = ((float)rand() / RAND_MAX) * 511;
			block_data.template get<1>(i*nslot + j) = ((float)rand() / RAND_MAX) * 511;
			block_data.template get<2>(i*nslot + j) = ((float)rand() / RAND_MAX) * 511;
		}
	}

	block_n.template get<0>(block_n.size()-1) = 0;
	block_insert.template hostToDevice<0>();
	block_n.template hostToDevice<0>();
	block_data.template hostToDevice<0,1,2>();

	mgpu::standard_context_t context(false);
	mgpu::scan((int *)block_n.template getDeviceBuffer<0>(), block_n.size(), (int *)block_n_scan.template getDeviceBuffer<0>() , context);

	block_n_scan.template deviceToHost<0>(block_n_scan.size()-1,block_n_scan.size()-1);
	size_t n_ele = block_n_scan.template get<0>(block_n.size()-1);
	cont_data.resize(n_ele);

	output_list_0.resize(n_ele);
	output_list_1.resize(n_ele);

	dim3 wthr;
	dim3 thr;

	wthr.x = nblock;
	wthr.y = 1;
	wthr.z = 1;
	thr.x = 64;
	thr.y = 1;
	thr.z = 1;

	construct_insert_list<<<wthr,thr>>>(block_insert.toKernel(),block_n.toKernel(),
										block_n_scan.toKernel(),output_list_0.toKernel(),output_list_1.toKernel(),
										block_data.toKernel(),cont_data.toKernel(),
										nslot);

	output_list_0.template deviceToHost<0>();
	block_n_scan.template deviceToHost<0>();
	cont_data.template deviceToHost<0,1,2>();

	// we check

	bool check = true;

	// fill block insert of some data
	for (int i = 0 ; i < nblock ; i++)
	{
		for (int j = 0 ; j < block_n.template get<0>(i) ; j++)
		{
			check &= output_list_0.template get<0>(block_n_scan.template get<0>(i) + j) == block_insert.template get<0>(i*nslot + j);

			check &= cont_data.template get<0>(block_n_scan.template get<0>(i) + j) == block_data.template get<0>(i*nslot + j);
			check &= cont_data.template get<1>(block_n_scan.template get<0>(i) + j) == block_data.template get<1>(i*nslot + j);
			check &= cont_data.template get<2>(block_n_scan.template get<0>(i) + j) == block_data.template get<2>(i*nslot + j);
		}
	}

	BOOST_REQUIRE_EQUAL(check,true);
}

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_kernels_merge_use )
{
	openfpm::vector_gpu<aggregate<int,int>> vct_index_old;
	openfpm::vector_gpu<aggregate<int,int>> vct_index;
	openfpm::vector_gpu<aggregate<int,int>> vct_add_index;

	openfpm::vector_gpu<aggregate<int,float,double>> vct_data;
	openfpm::vector_gpu<aggregate<int,float,double>> vct_add_data;

	vct_index_old.resize(1024);

	for (size_t i = 0 ; i < vct_index_old.size() ; i++)
	{
		vct_index_old.template get<0>(i) = 17*(float)rand()/RAND_MAX + i * 17;
		vct_index_old.template get<1>(i) = i;
	}

	vct_add_index.resize(100);

	for (size_t i = 0 ; i < vct_add_index.size() ; i++)
	{
		vct_add_index.template get<0>(i) = 17*(float)rand()/RAND_MAX + i * 17;
		vct_add_index.template get<1>(i) = i+vct_index_old.size();
	}

	// Now we merge

	vct_index.resize(vct_add_index.size() + vct_index_old.size());

	mgpu::standard_context_t ctx(false);

	// host to device
	vct_index_old.template hostToDevice<0,1>();
	vct_add_index.template hostToDevice<0,1>();

	mgpu::merge((int *)vct_index_old.template getDeviceBuffer<0>(),(int *)vct_index_old.template getDeviceBuffer<1>(),vct_index_old.size(),
			    (int *)vct_add_index.template getDeviceBuffer<0>(),(int *)vct_add_index.template getDeviceBuffer<1>(),vct_add_index.size(),
			    (int *)vct_index.template getDeviceBuffer<0>(),(int *)vct_index.template getDeviceBuffer<1>(),mgpu::less_t<int>(),ctx);

	vct_index.template deviceToHost<0,1>();

	// Check

	bool match = true;

	size_t i;
	for ( i = 0 ; i < vct_index.size() - 1 ; i++)
	{
		match &= vct_index.template get<0>(i) <= vct_index.template get<0>(i+1);

		int a = vct_index.template get<1>(i);

		if (a >= vct_index_old.size())
		{
			a -= vct_index_old.size();
			match &= vct_add_index.template get<0>(a) == vct_index.template get<0>(i);
		}
		else
		{
			match &= vct_index_old.template get<0>(a) == vct_index.template get<0>(i);
		}
	}

	int a = vct_index.template get<1>(i);

	if (a >= vct_index_old.size())
	{
		a -= vct_index_old.size();
		match &= vct_add_index.template get<0>(a) == vct_index.template get<0>(i);
	}
	else
	{
		match &= vct_index_old.template get<0>(a) == vct_index.template get<0>(i);
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_kernels_solve_conflicts_use )
{
	openfpm::vector_gpu<aggregate<int,int>> vct_index_old;
	openfpm::vector_gpu<aggregate<int>> vct_index;
	openfpm::vector_gpu<aggregate<int>> merge_indexes;
	openfpm::vector_gpu<aggregate<int,int>> vct_add_index;
	openfpm::vector_gpu<aggregate<int>> vct_index_out;

	openfpm::vector_gpu<aggregate<int,float,double>> vct_data_old;
	openfpm::vector_gpu<aggregate<int,float,double>> vct_add_data;

	openfpm::vector_gpu<aggregate<int,float,double>> vct_data_out;
	openfpm::vector_gpu<aggregate<int,int,int>> vct_tot_out;

	vct_index_old.resize(1024);
	vct_data_old.resize(vct_index_old.size());

	for (size_t i = 0 ; i < vct_index_old.size() ; i++)
	{
		vct_index_old.template get<0>(i) = 17*(float)rand()/RAND_MAX + i * 17;
		vct_index_old.template get<1>(i) = i;

		vct_data_old.template get<0>(i) = 128*(float)rand()/RAND_MAX;
		vct_data_old.template get<1>(i) = 128*(float)rand()/RAND_MAX;
		vct_data_old.template get<2>(i) = 128*(float)rand()/RAND_MAX;
	}

	vct_add_index.resize(100);
	vct_add_data.resize(100);

	for (size_t i = 0 ; i < vct_add_index.size() ; i++)
	{
		vct_add_index.template get<0>(i) = 17*(float)rand()/RAND_MAX + i * 17;
		vct_add_index.template get<1>(i) = i+vct_index_old.size();

		vct_add_data.template get<0>(i) = 128*(float)rand()/RAND_MAX;
		vct_add_data.template get<1>(i) = 128*(float)rand()/RAND_MAX;
		vct_add_data.template get<2>(i) = 128*(float)rand()/RAND_MAX;
	}

	// Now we merge

	vct_index.resize(vct_add_index.size() + vct_index_old.size());
	merge_indexes.resize(vct_index.size());

	mgpu::standard_context_t ctx(false);

	// host to device
	vct_index_old.template hostToDevice<0,1>();
	vct_add_index.template hostToDevice<0,1>();
	vct_data_old.template hostToDevice<0,1,2>();
	vct_add_data.template hostToDevice<0,1,2>();

	mgpu::merge((int *)vct_index_old.template getDeviceBuffer<0>(),(int *)vct_index_old.template getDeviceBuffer<1>(),vct_index_old.size(),
			    (int *)vct_add_index.template getDeviceBuffer<0>(),(int *)vct_add_index.template getDeviceBuffer<1>(),vct_add_index.size(),
			    (int *)vct_index.template getDeviceBuffer<0>(),(int *)merge_indexes.template getDeviceBuffer<0>(),mgpu::less_t<int>(),ctx);

	constexpr int bdim = 128;

	auto ite = vct_index.getGPUIterator(bdim);

	vct_index_out.resize(vct_index.size());
	vct_data_out.resize(vct_index.size());
	vct_tot_out.resize(ite.wthr.x);

	solve_conflicts<decltype(vct_index.toKernel()),decltype(vct_data_old.toKernel()),decltype(vct_tot_out.toKernel()),bdim,sadd_<0>,smin_<1>,smax_<2>><<<ite.wthr,ite.thr>>>(vct_index.toKernel(),vct_data_old.toKernel(),
						  merge_indexes.toKernel(),vct_add_data.toKernel(),
						  vct_index_out.toKernel(),vct_data_out.toKernel(),
						  vct_tot_out.toKernel(),vct_data_old.size());

	vct_index.template deviceToHost<0>();
	vct_data_out.template deviceToHost<0,1,2>();
	vct_tot_out.template deviceToHost<0>();
	vct_index_out.template deviceToHost<0>();
	merge_indexes.template deviceToHost<0>();

	size_t cnt = 0;

	bool match = true;

	while (cnt * bdim < vct_index.size())
	{
		size_t scan_block = 0;
		for (size_t i = 0 ; i < bdim - 1 && (cnt * bdim + i + 1) < vct_index.size() ; i++)
		{
			if (vct_index.template get<0>(cnt * bdim + i) == vct_index.template get<0>(cnt * bdim + i + 1))
			{
				// they match enable reduction pattern

				match &= vct_index_out.template get<0>(cnt * bdim + scan_block) == vct_index.template get<0>(cnt * bdim + i);

				int index1 = merge_indexes.template get<0>(cnt * bdim + i);
				int index2 = merge_indexes.template get<0>(cnt * bdim + i + 1) - vct_index_old.size();

				match &= vct_data_out.template get<0>(cnt * bdim + scan_block) == vct_data_old.template get<0>(index1) + vct_add_data.template get<0>(index2);

				match &= vct_data_out.template get<1>(cnt * bdim + scan_block) == smin_<1>::red(vct_data_old.template get<1>(index1),vct_add_data.template get<1>(index2));
				match &= vct_data_out.template get<2>(cnt * bdim + scan_block) == smax_<2>::red(vct_data_old.template get<2>(index1),vct_add_data.template get<2>(index2));

				i++;
			}
			else
			{
				match &= vct_index_out.template get<0>(cnt * bdim + scan_block) == vct_index.template get<0>(cnt * bdim + i);

				int index1 = merge_indexes.template get<0>(cnt * bdim + i);
				if (index1 < vct_index_old.size())
				{
					match &= vct_data_out.template get<0>(cnt * bdim + scan_block) == vct_data_old.template get<0>(index1);
					match &= vct_data_out.template get<1>(cnt * bdim + scan_block) == vct_data_old.template get<1>(index1);
					match &= vct_data_out.template get<2>(cnt * bdim + scan_block) == vct_data_old.template get<2>(index1);
				}
				else
				{
					index1 -= vct_index_old.size();
					match &= vct_data_out.template get<0>(cnt * bdim + scan_block) == vct_add_data.template get<0>(index1);
					match &= vct_data_out.template get<1>(cnt * bdim + scan_block) == vct_add_data.template get<1>(index1);
					match &= vct_data_out.template get<2>(cnt * bdim + scan_block) == vct_add_data.template get<2>(index1);
				}
			}
			scan_block++;
		}

		++cnt;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_kernels_realign_use )
{
	openfpm::vector_gpu<aggregate<int>> vct_index;
	openfpm::vector_gpu<aggregate<int,float,double>> vct_data;
	openfpm::vector_gpu<aggregate<int>> vct_index_out;
	openfpm::vector_gpu<aggregate<int,float,double>> vct_data_out;

	openfpm::vector_gpu<aggregate<int,int,int>> vct_tot_out;

	vct_index.resize(1024);
	vct_data.resize(vct_index.size());
	vct_tot_out.resize(vct_index.size() / 128);

	for (size_t i = 0 ; i < vct_index.size() ; i++)
	{
		vct_index.template get<0>(i) = 17*(float)rand()/RAND_MAX + i * 17;

		vct_data.template get<0>(i) = 128*(float)rand()/RAND_MAX;
		vct_data.template get<1>(i) = 128*(float)rand()/RAND_MAX;
		vct_data.template get<2>(i) = 128*(float)rand()/RAND_MAX;
	}

	for (size_t i = 0 ; i < vct_tot_out.size() ; i++)
	{
		vct_tot_out.template get<0>(i) = 128*(float)rand()/RAND_MAX;
		vct_tot_out.template get<2>(i) = 1;
	}

	vct_index.template hostToDevice<0>();
	vct_data.template hostToDevice<0,1,2>();
	vct_tot_out.template hostToDevice<0,2>();

	mgpu::standard_context_t ctx(false);
	mgpu::scan((int *)vct_tot_out.getDeviceBuffer<0>(),vct_tot_out.size(),(int *)vct_tot_out.getDeviceBuffer<1>(),ctx);

	vct_tot_out.deviceToHost<0,1>();
	vct_index_out.resize(vct_index.size());
	vct_data_out.resize(vct_index.size());

	int nblock = vct_tot_out.size();
	realign<<<nblock,128>>>(vct_index.toKernel(),vct_data.toKernel(),vct_index_out.toKernel(),vct_data_out.toKernel(),vct_tot_out.toKernel());

	vct_index_out.deviceToHost<0>();
	vct_data_out.deviceToHost<0,1,2>();

	int base = vct_tot_out.template get<1>(vct_tot_out.size()-1);
	int last = vct_tot_out.template get<0>(vct_tot_out.size()-1);

	vct_index_out.resize(base+last);

	size_t pr = 0;
	bool match = true;
	for (size_t i = 0 ; i < vct_tot_out.size() ; i++)
	{
		int tot = vct_tot_out.template get<0>(i);
		for (size_t j = 0 ; j < tot ; j++)
		{
			match &= vct_index_out.template get<0>(pr) == vct_index.template get<0>(i*128 + j);
			match &= vct_data_out.template get<0>(pr) == vct_data.template get<0>(i*128 + j);
			match &= vct_data_out.template get<1>(pr) == vct_data.template get<1>(i*128 + j);
			match &= vct_data_out.template get<2>(pr) == vct_data.template get<2>(i*128 + j);

			++pr;
		}
	}

	BOOST_REQUIRE_EQUAL(vct_index_out.template get<0>(vct_index_out.size()-1), vct_index.template get<0>(last + vct_index.size() - 128 -1 ) );
	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_SUITE_END()
