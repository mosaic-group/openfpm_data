//
// Created by tommaso on 30/05/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/BlockMapGpu.hpp"
#include "SparseGridGpu/BlockMapGpu_ker.cuh"
#include "SparseGridGpu/BlockMapGpu_kernels.cuh"
#include "SparseGridGpu/DataBlock.cuh"
#include "Vector/cuda/map_vector_sparse_cuda_kernels.cuh"

BOOST_AUTO_TEST_SUITE(BlockMapGpu_kernels_tests)

BOOST_AUTO_TEST_CASE(testSegreduce_total)
{
		typedef float ScalarT;
		typedef DataBlock<unsigned char, 64> MaskBlockT;
		typedef DataBlock<ScalarT, 64> BlockT;
		openfpm::vector_gpu<aggregate<int>> segments;
		openfpm::vector_gpu<aggregate<int>> segment_dataMap;
		openfpm::vector_gpu<aggregate<int>> outputMap;
		openfpm::vector_gpu<aggregate<int>> segments_oldData;

		// segment reduction for the new added data
		segments.resize(8);
		segments.template get<0>(0) = 0;
		segments.template get<0>(1) = 4;
		segments.template get<0>(2) = 5;
		segments.template get<0>(3) = 7;
		segments.template get<0>(4) = 8;
		segments.template get<0>(5) = 11;
		segments.template get<0>(6) = 17;
		segments.template get<0>(7) = 18; // Id of first non-existent data

		// segment old data indicate if new data has also old data
		// to be reduced with
		segments_oldData.resize(8);
		segments_oldData.template get<0>(0) = -1;
		segments_oldData.template get<0>(1) = 2;
		segments_oldData.template get<0>(2) = -1;
		segments_oldData.template get<0>(3) = 5;
		segments_oldData.template get<0>(4) = 7;
		segments_oldData.template get<0>(5) = -1;
		segments_oldData.template get<0>(6) = -1;
		segments_oldData.template get<0>(7) = -1; // Id of first non-existent data

		// for each added index we have a chunk in the vct_add_data vector
		// undortunately is not

		segment_dataMap.resize(19);
		segment_dataMap.template get<0>(0) = 10;
		segment_dataMap.template get<0>(1) = 1;
		segment_dataMap.template get<0>(2) = 50;
		segment_dataMap.template get<0>(3) = 11;
		segment_dataMap.template get<0>(4) = 13;
		segment_dataMap.template get<0>(5) = 87;
		segment_dataMap.template get<0>(6) = 54;
		segment_dataMap.template get<0>(7) = 33;
		segment_dataMap.template get<0>(8) = 22;
		segment_dataMap.template get<0>(9) = 17;
		segment_dataMap.template get<0>(10) = 40;
		segment_dataMap.template get<0>(11) = 32;
		segment_dataMap.template get<0>(12) = 80;
		segment_dataMap.template get<0>(13) = 52;
		segment_dataMap.template get<0>(14) = 21;
		segment_dataMap.template get<0>(15) = 76;
		segment_dataMap.template get<0>(16) = 65;
		segment_dataMap.template get<0>(17) = 54;
		segment_dataMap.template get<0>(18) = 3;

		outputMap.resize(7);
		outputMap.template get<0>(0) = 9;
		outputMap.template get<0>(1) = 11;
		outputMap.template get<0>(2) = 13;
		outputMap.template get<0>(3) = 34;
		outputMap.template get<0>(4) = 23;
		outputMap.template get<0>(5) = 90;
		outputMap.template get<0>(6) = 21;

		segments.template hostToDevice<0>();
		segment_dataMap.hostToDevice<0>();
		segments_oldData.hostToDevice<0>();
		outputMap.hostToDevice<0>();

		const unsigned int BITMASK = 0, BLOCK = 1;
		BlockT block;
		MaskBlockT mask;
		BlockT block_old;
		MaskBlockT mask_old;
		for (int i = 0; i < 32; ++i)
		{
			block[i] = i + 1;
			mask[i] = 1;
			block_old[i] = 100 + i + 1;
			mask_old[i] = 1;
		}
		for (int i = 32; i < 64; ++i)
		{
			block[i] = 666;
			mask[i] = 0;
			block_old[i] = 666;
			mask_old[i] = 0;
		}
		block_old[33] = 665;
		mask_old[33] = 1;

		openfpm::vector_gpu<aggregate<MaskBlockT, BlockT>> data_old;
		openfpm::vector_gpu<aggregate<MaskBlockT, BlockT>> data_new;
		data_new.resize(100);
		data_old.resize(100);
		for (int i = 0; i < 100; ++i)
		{
			data_new.template get<BITMASK>(i) = mask;
			data_new.template get<BLOCK>(i) = block;
			if (i < data_old.size())
			{
				data_old.template get<BITMASK>(i) = mask_old;
				data_old.template get<BLOCK>(i) = block_old;
			}
		}

		data_new.template hostToDevice<BITMASK, BLOCK>();
		data_old.template hostToDevice<BITMASK, BLOCK>();

		// Allocate output buffer
		openfpm::vector_gpu<aggregate<MaskBlockT, BlockT>> outputData;
		outputData.resize(100);

		CUDA_LAUNCH_DIM3((BlockMapGpuKernels::segreduce_total<BLOCK, 0, BITMASK, 2, mgpu::plus_t<ScalarT>>),segments.size()-1, 2*BlockT::size,
		data_new.toKernel(),
		data_old.toKernel(),
		segments.toKernel(),
		segment_dataMap.toKernel(),
		segments_oldData.toKernel(),
		outputMap.toKernel(),
		outputData.toKernel());

		// Segreduce on mask
		CUDA_LAUNCH_DIM3((BlockMapGpuKernels::segreduce_total<BITMASK, 0, BITMASK, 2, mgpu::maximum_t<unsigned char>>),segments.size()-1, 2*BlockT::size,
		data_new.toKernel(),
		data_old.toKernel(),
		segments.toKernel(),
		segment_dataMap.toKernel(),
		segments_oldData.toKernel(),
		outputMap.toKernel(),
		outputData.toKernel());

		// Check

        outputData.template deviceToHost<BITMASK, BLOCK>();

        std::cout << "------------------------------------------------------------" << std::endl;

		for (int j = 0; j < outputMap.size(); ++j)
		{
			size_t j_ = outputMap.template get<0>(j);

			std::cout << "Index: " << j_ << std::endl;

			BlockT outBlock = outputData.template get<BLOCK>(j_);
			MaskBlockT outMask = outputData.template get<BITMASK>(j_);

			for (int i = 0; i < BlockT::size; ++i)
			{
				std::cout << outBlock[i] << ":" << (int)outMask[i] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "------------------------------------------------------------" << std::endl;

        // Check

        for (int j = 0 ; j < outputMap.size() ; j++)
        {
        	int out_id = outputMap.template get<0>(j);
        	int seg_mult = segments.template get<0>(j+1) - segments.template get<0>(j);

			BlockT outBlock = outputData.template get<BLOCK>(out_id);
			MaskBlockT outMask = outputData.template get<BITMASK>(out_id);

        	if (segments_oldData.template get<0>(j) != -1)
        	{
    			for (int i = 0; i < BlockT::size / 2; ++i)
    			{
    				BOOST_REQUIRE_EQUAL(outMask[i],1);
    				BOOST_REQUIRE_EQUAL(outBlock[i],100 + (i+1)*(seg_mult+1));
    			}
				BOOST_REQUIRE_EQUAL(outMask[33],1);
				BOOST_REQUIRE_EQUAL(outBlock[33],665);
        	}
        	else
        	{
    			for (int i = 0; i < BlockT::size / 2; ++i)
    			{
    				BOOST_REQUIRE_EQUAL(outMask[i],1);
    				BOOST_REQUIRE_EQUAL(outBlock[i],(i+1)*seg_mult);
    			}
        	}
        }
}


BOOST_AUTO_TEST_SUITE_END() // SparseGridGpu_kernels_tests

BOOST_AUTO_TEST_SUITE(BlockMapGpu_functors_tests)


BOOST_AUTO_TEST_CASE(test_maps_create)
{
	openfpm::vector_gpu<aggregate<int>> merge_indexes;
	openfpm::vector_gpu<aggregate<int>> merge_keys;

	openfpm::vector_gpu<aggregate<int,int,int,int,int>> p_ids;
	openfpm::vector_gpu<aggregate<int,int,int,int,int>> s_ids;

	openfpm::vector_gpu<aggregate<int>> copy_old_src;
	openfpm::vector_gpu<aggregate<int>> copy_old_dst;
	openfpm::vector_gpu<aggregate<int>> outputMap;
	openfpm::vector_gpu<aggregate<int>> segments_oldData;

	merge_keys.resize(16);
	merge_keys.template get<0>(0) = 22;      // new
	merge_keys.template get<0>(1) = 23;     // new
	merge_keys.template get<0>(2) = 33;      // old id:0
	merge_keys.template get<0>(3) = 34;      // old id:1
	merge_keys.template get<0>(4) = 36;      // old id:2
	merge_keys.template get<0>(5) = 37;      // old id:3
	merge_keys.template get<0>(6) = 43;      // old id:4
	merge_keys.template get<0>(7) = 45;      // old id:5
	merge_keys.template get<0>(8) = 46;     // new
	merge_keys.template get<0>(9) = 56;      // old id:6
	merge_keys.template get<0>(10) = 56;     // new
	merge_keys.template get<0>(11) = 60;     // old id:7
	merge_keys.template get<0>(12) = 61;     // old id: 8
	merge_keys.template get<0>(13) = 63;    // new
	merge_keys.template get<0>(14) = 64;    // new
	merge_keys.template get<0>(15) = 65;    // new

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

	merge_keys.template hostToDevice<0>();
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

	CUDA_LAUNCH(BlockMapGpuKernels::compute_predicate,ite,merge_keys.toKernel(),merge_indexes.toKernel(),9,p_ids.toKernel());

	mgpu::ofp_context_t context;
	openfpm::scan((int *)p_ids.template getDeviceBuffer<0>(),
				s_ids.size(),
	            (int *)s_ids.template getDeviceBuffer<0>(),
                context);

	openfpm::scan((int *)p_ids.template getDeviceBuffer<1>(),
				s_ids.size(),
	            (int *)s_ids.template getDeviceBuffer<1>(),
                context);

	openfpm::scan((int *)p_ids.template getDeviceBuffer<2>(),
				s_ids.size(),
	            (int *)s_ids.template getDeviceBuffer<2>(),
                context);

	openfpm::scan((int *)p_ids.template getDeviceBuffer<3>(),
				s_ids.size(),
	            (int *)s_ids.template getDeviceBuffer<3>(),
                context);

	openfpm::scan((int *)p_ids.template getDeviceBuffer<4>(),
				s_ids.size(),
	            (int *)s_ids.template getDeviceBuffer<4>(),
                context);

	s_ids.template deviceToHost<0,1,2,3,4>();
	p_ids.template deviceToHost<0,1,2,3,4>();

	size_t copy_old_size = s_ids.template get<3>(s_ids.size()-1) + p_ids.template get<3>(p_ids.size()-1);
	size_t seg_old_size = s_ids.template get<1>(s_ids.size()-1) + p_ids.template get<1>(p_ids.size()-1);
	size_t out_map_size = s_ids.template get<1>(s_ids.size()-1) + p_ids.template get<1>(p_ids.size()-1);

	segments_oldData.resize(seg_old_size);
	outputMap.resize(out_map_size);
	copy_old_src.resize(copy_old_size);
	copy_old_dst.resize(copy_old_size);

	CUDA_LAUNCH(BlockMapGpuKernels::maps_create,ite,s_ids.toKernel(),p_ids.toKernel(),segments_oldData.toKernel(),outputMap.toKernel(),copy_old_dst.toKernel(),copy_old_src.toKernel());

	segments_oldData.template deviceToHost<0>();
	outputMap.template deviceToHost<0>();
	copy_old_dst.template deviceToHost<0>();
	copy_old_src.template deviceToHost<0>();

	BOOST_REQUIRE_EQUAL(seg_old_size,7);
	BOOST_REQUIRE_EQUAL(out_map_size,7);
	BOOST_REQUIRE_EQUAL(copy_old_size,8);

	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(0),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(0),0);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(1),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(1),1);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(2),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(2),8);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(3),6);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(3),9);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(4),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(4),12);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(5),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(5),13);
	BOOST_REQUIRE_EQUAL(segments_oldData.template get<0>(6),-1);
	BOOST_REQUIRE_EQUAL(outputMap.template get<0>(6),14);

	BOOST_REQUIRE_EQUAL(copy_old_dst.template get<0>(0),2);
	BOOST_REQUIRE_EQUAL(copy_old_dst.template get<0>(1),3);
	BOOST_REQUIRE_EQUAL(copy_old_dst.template get<0>(2),4);
	BOOST_REQUIRE_EQUAL(copy_old_dst.template get<0>(3),5);
	BOOST_REQUIRE_EQUAL(copy_old_dst.template get<0>(4),6);
	BOOST_REQUIRE_EQUAL(copy_old_dst.template get<0>(5),7);
	BOOST_REQUIRE_EQUAL(copy_old_dst.template get<0>(6),10);
	BOOST_REQUIRE_EQUAL(copy_old_dst.template get<0>(7),11);

	BOOST_REQUIRE_EQUAL(copy_old_src.template get<0>(0),0);
	BOOST_REQUIRE_EQUAL(copy_old_src.template get<0>(1),1);
	BOOST_REQUIRE_EQUAL(copy_old_src.template get<0>(2),2);
	BOOST_REQUIRE_EQUAL(copy_old_src.template get<0>(3),3);
	BOOST_REQUIRE_EQUAL(copy_old_src.template get<0>(4),4);
	BOOST_REQUIRE_EQUAL(copy_old_src.template get<0>(5),5);
	BOOST_REQUIRE_EQUAL(copy_old_src.template get<0>(6),7);
	BOOST_REQUIRE_EQUAL(copy_old_src.template get<0>(7),8);
}

BOOST_AUTO_TEST_CASE (testSolve_conflicts)
{
	typedef float ScalarT;
	typedef DataBlock<ScalarT, 64> BlockT;
	typedef DataBlock<unsigned char, 64> MaskBlockT;
	const unsigned int p=0, pMask=1, pInd=0;
	openfpm::vector_gpu<aggregate<unsigned int>> keys, mergeIndices, tmpIndices, keysOut, trivial_map;
	openfpm::vector_gpu<aggregate<unsigned int,unsigned int>> segments_new;
	openfpm::vector_gpu<aggregate<BlockT, MaskBlockT>> dataOld, dataNew, tmpData, dataOut;
	mgpu::ofp_context_t ctx;

	// Keys
	keys.resize(14);
	keys.template get<pInd>(0) = 0;
	keys.template get<pInd>(1) = 1;
	keys.template get<pInd>(2) = 2;
	keys.template get<pInd>(3) = 4;
	keys.template get<pInd>(4) = 4;
	keys.template get<pInd>(5) = 5;
	keys.template get<pInd>(6) = 6;
	keys.template get<pInd>(7) = 7;
	keys.template get<pInd>(8) = 8;
	keys.template get<pInd>(9) = 9;
	keys.template get<pInd>(10) = 10;
	keys.template get<pInd>(11) = 11;
	keys.template get<pInd>(12) = 13;
	keys.template get<pInd>(13) = 13;

	// Merge Indices
	mergeIndices.resize(14);
	mergeIndices.template get<pInd>(0) = 5;  // 0
	mergeIndices.template get<pInd>(1) = 0;
	mergeIndices.template get<pInd>(2) = 6;  // 1
	mergeIndices.template get<pInd>(3) = 1;
	mergeIndices.template get<pInd>(4) = 7;  // 2
	mergeIndices.template get<pInd>(5) = 8;  // 3
	mergeIndices.template get<pInd>(6) = 9;  // 4
	mergeIndices.template get<pInd>(7) = 10; // 5
	mergeIndices.template get<pInd>(8) = 2;
	mergeIndices.template get<pInd>(9) = 11; // 6
	mergeIndices.template get<pInd>(10) = 3;
	mergeIndices.template get<pInd>(11) = 12; // 7
	mergeIndices.template get<pInd>(12) = 13; // 8
	mergeIndices.template get<pInd>(13) = 14; // 9

	// segments new
	segments_new.resize(10);
	segments_new.template get<0>(0) = 0;   // 5
	segments_new.template get<0>(1) = 1;   // 6
	segments_new.template get<0>(2) = 2;   // 7
	segments_new.template get<0>(3) = 3;   // 8
	segments_new.template get<0>(4) = 4;   // 9
	segments_new.template get<0>(5) = 5;   // 10
	segments_new.template get<0>(6) = 6;   // 11
	segments_new.template get<0>(7) = 7;   // 12
	segments_new.template get<0>(8) = 8;   // 13,14
	segments_new.template get<0>(9) = 10;

	// Fill the data
	{ // We want i to live in a confined namespace
		int i = 0;
		dataOld.resize(5);
		for (; i < dataOld.size(); ++i)
		{
			for (int j = 0; j < BlockT::size; ++j)
			{
				dataOld.template get<p>(i)[j] = i;
				dataOld.template get<pMask>(i)[j] = 1;
			};
			dataOld.template get<p>(i)[0] = 1;
		}
		dataNew.resize(10);
		for (; i < dataOld.size() + dataNew.size(); ++i)
		{
			int ii = i - dataOld.size();
			for (int j = 0; j < BlockT::size; ++j)
			{
				dataNew.template get<p>(ii)[j] = i;
				dataNew.template get<pMask>(ii)[j] = 1;
			};
			dataNew.template get<p>(ii)[0] = 1;
		}
	}

    //// Create a trivial map

    trivial_map.resize(dataNew.size()+1);

    for (size_t i = 0 ; i < trivial_map.size(); i++)
    {
    	trivial_map.template get<0>(i) = i;
    }

    trivial_map.hostToDevice<0>();

	// Copy to device
	keys.hostToDevice<pInd>();
	mergeIndices.hostToDevice<pInd>();
	dataOld.hostToDevice<p, pMask>();
	dataNew.hostToDevice<p, pMask>();
	segments_new.hostToDevice<0>();

	BlockMapGpuFunctors::BlockFunctor<128> obj;

	// Now perform the compaction

	obj.solve_conflicts<0,
			decltype(keys),
			decltype(segments_new),
			decltype(dataOld),
			sadd_<p>
			>(
					keys, mergeIndices, segments_new, trivial_map,
					dataOld, dataNew,
					keysOut, dataOut,
					ctx
					);

	// Now retrieve the dataDst vector
	keysOut.deviceToHost<pInd>();
	dataOut.deviceToHost<p, pMask>();

	// Debug output
	for (int i=0; i<dataOut.size(); ++i)
	{
		std::cout
			<< "dataOut[" << i << "][0] = " << dataOut.template get<p>(i)[0]
			<< ", dataOut[" << i << "][1] = " << dataOut.template get<p>(i)[1]
			<< std::endl;
	}

	// Validation
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(0)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(1)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(2)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(3)[0], 2);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(4)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(5)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(6)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(7)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(8)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(9)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(10)[0], 1);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(11)[0], 2);

	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(0)[1], 5);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(1)[1], 0);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(2)[1], 6);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(3)[1], 8);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(4)[1], 8);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(5)[1], 9);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(6)[1], 10);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(7)[1], 2);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(8)[1], 11);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(9)[1], 3);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(10)[1], 12);
	BOOST_REQUIRE_EQUAL(dataOut.template get<p>(11)[1], 27);
}

BOOST_AUTO_TEST_SUITE_END() //SparseGridGpu_functors_tests

