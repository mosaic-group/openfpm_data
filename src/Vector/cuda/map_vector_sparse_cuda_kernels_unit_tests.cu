#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <Vector/map_vector.hpp>
#include <Vector/cuda/map_vector_sparse_cuda_kernels.cuh>
#include "util/cuda/scan_ofp.cuh"
#include "util/cuda/merge_ofp.cuh"

BOOST_AUTO_TEST_SUITE( vector_sparse_cuda_kernels )

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_kernels_use )
{
	openfpm::vector_gpu<aggregate<int>> block_insert;
	openfpm::vector_gpu<aggregate<int>> block_n;
	openfpm::vector_gpu<aggregate<int>> block_n_scan;
	openfpm::vector_gpu<aggregate<int>> output_list_0;
	openfpm::vector_gpu<aggregate<int>> output_list_1;

	int nblock = 100;
	int nslot = 512;
	block_insert.resize(nblock*nslot);
	block_n.resize(nblock+1);
	block_n_scan.resize(nblock+1);

	// fill block insert of some data
	for (int i = 0 ; i < nblock ; i++)
	{
		block_n.template get<0>(i) = ((float)rand() / (float)RAND_MAX) * 511;
		for (int j = 0 ; j < block_n.template get<0>(i) ; j++)
		{
			block_insert.template get<0>(i*nslot + j) = ((float)rand() / (float)RAND_MAX) * 511;
		}
	}

	block_n.template get<0>(block_n.size()-1) = 0;
	block_insert.template hostToDevice<0>();
	block_n.template hostToDevice<0>();

	mgpu::ofp_context_t context;
	openfpm::scan((int *)block_n.template getDeviceBuffer<0>(), block_n.size(), (int *)block_n_scan.template getDeviceBuffer<0>() , context);

	block_n_scan.template deviceToHost<0>(block_n_scan.size()-1,block_n_scan.size()-1);
	size_t n_ele = block_n_scan.template get<0>(block_n.size()-1);

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

	CUDA_LAUNCH_DIM3(construct_insert_list_key_only,wthr,thr,block_insert.toKernel(),block_n.toKernel(),
										block_n_scan.toKernel(),output_list_0.toKernel(),output_list_1.toKernel(),
										nslot);

	output_list_0.template deviceToHost<0>();
	block_n_scan.template deviceToHost<0>();

	// we check

	bool check = true;

	// fill block insert of some data
	for (int i = 0 ; i < nblock ; i++)
	{
		for (int j = 0 ; j < block_n.template get<0>(i) ; j++)
		{
			check &= output_list_0.template get<0>(block_n_scan.template get<0>(i) + j) == block_insert.template get<0>(i*nslot + j);
		}
	}

	BOOST_REQUIRE_EQUAL(check,true);
}

BOOST_AUTO_TEST_CASE( vector_sparse_cuda_kernels_use_small_pool )
{
	openfpm::vector_gpu<aggregate<int>> block_insert;
	openfpm::vector_gpu<aggregate<int>> block_n;
	openfpm::vector_gpu<aggregate<int>> block_n_scan;
	openfpm::vector_gpu<aggregate<int>> output_list_0;
	openfpm::vector_gpu<aggregate<int>> output_list_1;

	int nblock = 100;
	int nslot = 17;
	block_insert.resize(nblock*nslot);
	block_n.resize(nblock+1);
	block_n_scan.resize(nblock+1);

	// fill block insert of some data
	for (int i = 0 ; i < nblock ; i++)
	{
		block_n.template get<0>(i) = ((float)rand() / (float)RAND_MAX) * 16;
		for (int j = 0 ; j < block_n.template get<0>(i) ; j++)
		{
			block_insert.template get<0>(i*nslot + j) = ((float)rand() / (float)RAND_MAX) * 511;
		}
	}

	block_n.template get<0>(block_n.size()-1) = 0;
	block_insert.template hostToDevice<0>();
	block_n.template hostToDevice<0>();

	mgpu::ofp_context_t context;
	openfpm::scan((int *)block_n.template getDeviceBuffer<0>(), block_n.size(), (int *)block_n_scan.template getDeviceBuffer<0>() , context);

	block_n_scan.template deviceToHost<0>(block_n_scan.size()-1,block_n_scan.size()-1);
	size_t n_ele = block_n_scan.template get<0>(block_n.size()-1);

	output_list_0.resize(n_ele);
	output_list_1.resize(n_ele);

	auto ite = block_insert.getGPUIterator();

	CUDA_LAUNCH(construct_insert_list_key_only_small_pool,ite,block_insert.toKernel(),block_n.toKernel(),
										block_n_scan.toKernel(),output_list_0.toKernel(),output_list_1.toKernel(),
										nslot);

	output_list_0.template deviceToHost<0>();
	block_n_scan.template deviceToHost<0>();

	// we check

	bool check = true;

	// fill block insert of some data
	for (int i = 0 ; i < nblock ; i++)
	{
		for (int j = 0 ; j < block_n.template get<0>(i) ; j++)
		{
			check &= output_list_0.template get<0>(block_n_scan.template get<0>(i) + j) == block_insert.template get<0>(i*nslot + j);
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
		vct_index_old.template get<0>(i) = 17*(float)rand()/(float)RAND_MAX + i * 17;
		vct_index_old.template get<1>(i) = i;
	}

	vct_add_index.resize(100);

	for (size_t i = 0 ; i < vct_add_index.size() ; i++)
	{
		vct_add_index.template get<0>(i) = 17*(float)rand()/(float)RAND_MAX + i * 17;
		vct_add_index.template get<1>(i) = i+vct_index_old.size();
	}

	// Now we merge

	vct_index.resize(vct_add_index.size() + vct_index_old.size());

	mgpu::ofp_context_t ctx;

	// host to device
	vct_index_old.template hostToDevice<0,1>();
	vct_add_index.template hostToDevice<0,1>();

	openfpm::merge((int *)vct_index_old.template getDeviceBuffer<0>(),(int *)vct_index_old.template getDeviceBuffer<1>(),vct_index_old.size(),
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
		vct_index_old.template get<0>(i) = 17*(float)rand()/(float)RAND_MAX + i * 17;
		vct_index_old.template get<1>(i) = i;

		vct_data_old.template get<0>(i) = 128*(float)rand()/(float)RAND_MAX;
		vct_data_old.template get<1>(i) = 128*(float)rand()/(float)RAND_MAX;
		vct_data_old.template get<2>(i) = 128*(float)rand()/(float)RAND_MAX;
	}

	vct_add_index.resize(100);
	vct_add_data.resize(100);

	for (size_t i = 0 ; i < vct_add_index.size() ; i++)
	{
		vct_add_index.template get<0>(i) = 17*(float)rand()/(float)RAND_MAX + i * 17;
		vct_add_index.template get<1>(i) = i+vct_index_old.size();

		vct_add_data.template get<0>(i) = 128*(float)rand()/(float)RAND_MAX;
		vct_add_data.template get<1>(i) = 128*(float)rand()/(float)RAND_MAX;
		vct_add_data.template get<2>(i) = 128*(float)rand()/(float)RAND_MAX;
	}

	// Now we merge

	vct_index.resize(vct_add_index.size() + vct_index_old.size());
	merge_indexes.resize(vct_index.size());

	mgpu::ofp_context_t ctx;

	// host to device
	vct_index_old.template hostToDevice<0,1>();
	vct_add_index.template hostToDevice<0,1>();
	vct_data_old.template hostToDevice<0,1,2>();
	vct_add_data.template hostToDevice<0,1,2>();

	openfpm::merge((int *)vct_index_old.template getDeviceBuffer<0>(),(int *)vct_index_old.template getDeviceBuffer<1>(),vct_index_old.size(),
			    (int *)vct_add_index.template getDeviceBuffer<0>(),(int *)vct_add_index.template getDeviceBuffer<1>(),vct_add_index.size(),
			    (int *)vct_index.template getDeviceBuffer<0>(),(int *)merge_indexes.template getDeviceBuffer<0>(),mgpu::less_t<int>(),ctx);

	constexpr int bdim = 128;

	auto ite = vct_index.getGPUIterator(bdim);

	vct_index_out.resize(vct_index.size());
	vct_data_out.resize(vct_index.size());
	vct_tot_out.resize(ite.wthr.x);

	CUDA_LAUNCH_DIM3((solve_conflicts<decltype(vct_index.toKernel()),decltype(vct_data_old.toKernel()),decltype(vct_tot_out.toKernel()),bdim,sadd_<0>,smin_<1>,smax_<2>>),ite.wthr,ite.thr,vct_index.toKernel(),vct_data_old.toKernel(),
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
		vct_index.template get<0>(i) = 17*(float)rand()/(float)RAND_MAX + i * 17;

		vct_data.template get<0>(i) = 128*(float)rand()/(float)RAND_MAX;
		vct_data.template get<1>(i) = 128*(float)rand()/(float)RAND_MAX;
		vct_data.template get<2>(i) = 128*(float)rand()/(float)RAND_MAX;
	}

	for (size_t i = 0 ; i < vct_tot_out.size() ; i++)
	{
		vct_tot_out.template get<0>(i) = 128*(float)rand()/(float)RAND_MAX;
		vct_tot_out.template get<2>(i) = 1;
	}

	vct_index.template hostToDevice<0>();
	vct_data.template hostToDevice<0,1,2>();
	vct_tot_out.template hostToDevice<0,2>();

	mgpu::ofp_context_t ctx;
	openfpm::scan((int *)vct_tot_out.getDeviceBuffer<0>(),vct_tot_out.size(),(int *)vct_tot_out.getDeviceBuffer<1>(),ctx);

	vct_tot_out.deviceToHost<0,1>();
	vct_index_out.resize(vct_index.size());
	vct_data_out.resize(vct_index.size());

	int nblock = vct_tot_out.size();
	CUDA_LAUNCH_DIM3(realign,nblock,128,vct_index.toKernel(),vct_data.toKernel(),vct_index_out.toKernel(),vct_data_out.toKernel(),vct_tot_out.toKernel());

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

			if (match == false)
			{
				std::cout << 0 << " " << pr << " " << i*128 + j << "  " << vct_index_out.template get<0>(pr) << "  " << vct_index.template get<0>(i*128 + j) << std::endl;
				std::cout << 1 << " " << pr << " " << i*128 + j << "  " << vct_data_out.template get<0>(pr) << "  " << vct_data.template get<0>(i*128 + j) << std::endl;
				std::cout << 2 << " " << pr << " " << i*128 + j << "  " << vct_data_out.template get<1>(pr) << "  " << vct_data.template get<1>(i*128 + j) << std::endl;
				std::cout << 3 << " " << pr << " " << i*128 + j << "  " << vct_data_out.template get<2>(pr) << "  " << vct_data.template get<2>(i*128 + j) << std::endl;
			}

			++pr;
		}
	}

	BOOST_REQUIRE_EQUAL(vct_index_out.template get<0>(vct_index_out.size()-1), vct_index.template get<0>(last + vct_index.size() - 128 -1 ) );
	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_SUITE_END()
