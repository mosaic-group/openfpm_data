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

BOOST_AUTO_TEST_CASE (testSegreduce)
{
		typedef float ScalarT;
		typedef DataBlock<unsigned char, 64> MaskBlockT;
		typedef DataBlock<ScalarT, 64> BlockT;
		openfpm::vector_gpu<aggregate<int>> segments;
		segments.resize(8);
		segments.template get<0>(0) = 0;
		segments.template get<0>(1) = 4;
		segments.template get<0>(2) = 5;
		segments.template get<0>(3) = 7;
		segments.template get<0>(4) = 8;
		segments.template get<0>(5) = 11;
		segments.template get<0>(6) = 17;
		segments.template get<0>(7) = 18; // Id of first non-existent data

		segments.template hostToDevice<0>();

		const unsigned int BITMASK = 0, BLOCK = 1;
		BlockT block;
		MaskBlockT mask;
		for (int i = 0; i < 32; ++i)
		{
			block[i] = i + 1;
			mask[i] = 1;
		}
		for (int i = 32; i < 64; ++i)
		{
			block[i] = 666;
			mask[i] = 0;
		}

		openfpm::vector_gpu<aggregate<MaskBlockT, BlockT>> data;
		data.resize(18);
		for (int i = 0; i < 18; ++i)
		{
			data.template get<BITMASK>(i) = mask;
			data.template get<BLOCK>(i) = block;
		}

		data.template hostToDevice<BITMASK, BLOCK>();

		// Allocate output buffer
		openfpm::vector_gpu<aggregate<MaskBlockT, BlockT>> outputData;
		outputData.resize(segments.size()-1);

		BlockMapGpuKernels::segreduce<BLOCK, 0, BITMASK, 2, mgpu::plus_t<ScalarT>> <<< outputData.size(), 2*BlockT::size >>> (
		data.toKernel(),
		segments.toKernel(),
		data.toKernel(),
		outputData.toKernel()
		);

		// Segreduce on mask
		BlockMapGpuKernels::segreduce<BITMASK, 0, BITMASK, 2, mgpu::maximum_t<unsigned char>> <<< outputData.size(), 2*BlockT::size >>> (
		data.toKernel(),
		segments.toKernel(),
		data.toKernel(),
		outputData.toKernel()
		);

		outputData.template deviceToHost<BITMASK, BLOCK>();

		for (int j = 0; j < outputData.size(); ++j)
		{
			BlockT outBlock = outputData.template get<BLOCK>(j);
			MaskBlockT outMask = outputData.template get<BITMASK>(j);

			for (int i = 0; i < BlockT::size; ++i)
			{
				std::cout << outBlock[i] << ":" << (int)outMask[i] << " ";
			}
			std::cout << std::endl;
		}
}

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

		segment_dataMap.resize(20);
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


		for (int j = 0; j < outputData.size(); ++j)
		{
			BlockT outBlock = outputData.template get<BLOCK>(j);
			MaskBlockT outMask = outputData.template get<BITMASK>(j);

			for (int i = 0; i < BlockT::size; ++i)
			{
				std::cout << outBlock[i] << ":" << (int)outMask[i] << " ";
			}
			std::cout << std::endl;
		}

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

BOOST_AUTO_TEST_CASE (testSegreduce_vector)
{
		typedef float ScalarT;
		typedef DataBlock<unsigned char, 64> MaskBlockT;
		typedef DataBlock<ScalarT, 64> BlockT;
		openfpm::vector_gpu<aggregate<int>> segments;
		segments.resize(8);
		segments.template get<0>(0) = 0;
		segments.template get<0>(1) = 4;
		segments.template get<0>(2) = 5;
		segments.template get<0>(3) = 7;
		segments.template get<0>(4) = 8;
		segments.template get<0>(5) = 11;
		segments.template get<0>(6) = 17;
		segments.template get<0>(7) = 18; // Id of first non-existent data

		segments.template hostToDevice<0>();

		const unsigned int BITMASK = 0, BLOCK = 1;
		BlockT block[3];
		MaskBlockT mask;
		for (int i = 0; i < 32; ++i)
		{
			block[0][i] = i + 1;
			block[1][i] = i + 1;
			block[2][i] = i + 1;
			mask[i] = 1;
		}
		for (int i = 32; i < 64; ++i)
		{
			block[0][i] = 666;
			block[1][i] = 666;
			block[2][i] = 666;
			mask[i] = 0;
		}

		openfpm::vector_gpu<aggregate<MaskBlockT, BlockT[3]>> data;
		data.resize(18);
		for (int i = 0; i < 18; ++i)
		{
			data.template get<BITMASK>(i) = mask;
			for (int k=0; k<3; k++)
			{
				data.template get<BLOCK>(i)[k] = block[k];
			}
		}

		data.template hostToDevice<BITMASK, BLOCK>();

		// Allocate output buffer
		openfpm::vector_gpu<aggregate<MaskBlockT, BlockT[3]>> outputData;
		outputData.resize(segments.size()-1);

		BlockMapGpuKernels::segreduce<BLOCK, 0, BITMASK, 2, mgpu::plus_t<ScalarT>> <<< outputData.size(), 2*BlockT::size >>> (
		data.toKernel(),
		segments.toKernel(),
		data.toKernel(),
		outputData.toKernel()
		);

		// Segreduce on mask
		BlockMapGpuKernels::segreduce<BITMASK, 0, BITMASK, 2, mgpu::maximum_t<unsigned char>> <<< outputData.size(), 2*BlockT::size >>> (
		data.toKernel(),
		segments.toKernel(),
		data.toKernel(),
		outputData.toKernel()
		);

		outputData.template deviceToHost<BITMASK, BLOCK>();

		for (int j = 0; j < outputData.size(); ++j)
		{
			BlockT outBlock[3];
			for (int k = 0 ; k < 3 ;k++)
			{outBlock[k] = outputData.template get<BLOCK>(j)[k];}
			MaskBlockT outMask = outputData.template get<BITMASK>(j);

			for (int i = 0; i < BlockT::size; ++i)
			{
				std::cout << "(" << outBlock[0][i] << "," << outBlock[1][i] << "," << outBlock[2][i] << ")" << ":" << (int)outMask[i] << " ";
			}
			std::cout << std::endl;
		}
}

//    BOOST_AUTO_TEST_CASE (testSolveConflicts)
//    {
//        typedef float ScalarT;
//        typedef DataBlock<unsigned char, 64> MaskBlockT;
//        typedef DataBlock<ScalarT, 64> BlockT;
//        openfpm::vector_gpu<aggregate<int>> segments;
//        segments.resize(8);
//        segments.template get<0>(0) = 0;
//        segments.template get<0>(1) = 4;
//        segments.template get<0>(2) = 5;
//        segments.template get<0>(3) = 7;
//        segments.template get<0>(4) = 8;
//        segments.template get<0>(5) = 11;
//        segments.template get<0>(6) = 17;
//        segments.template get<0>(7) = 18; // Id of first non-existent data
//
//        segments.template hostToDevice<0>();
//
//        const unsigned int BITMASK = 0, BLOCK = 1;
//        BlockT block;
//        MaskBlockT mask;
//        for (int i = 0; i < 32; ++i)
//        {
//            block[i] = i + 1;
//            mask[i] = 1;
//        }
//        for (int i = 32; i < 64; ++i)
//        {
//            block[i] = 666;
//            mask[i] = 0;
//        }
//
//        openfpm::vector_gpu<aggregate<MaskBlockT, BlockT>> data;
//        data.resize(18);
//        for (int i = 0; i < 18; ++i)
//        {
//            data.template get<BITMASK>(i) = mask;
//            data.template get<BLOCK>(i) = block;
//        }
//
//        data.template hostToDevice<BITMASK, BLOCK>();
//
//        // Allocate output buffer
//        openfpm::vector_gpu<aggregate<MaskBlockT, BlockT>> outputData;
//        outputData.resize(segments.size()-1);
//
//        BlockMapGpuKernels::segreduce<BLOCK, 0, BITMASK, 2, mgpu::plus_t<ScalarT>> <<< outputData.size(), 2*BlockT::size >>> (
//                data.toKernel(),
//                        segments.toKernel(),
//                        data.toKernel(),
//                        outputData.toKernel()
//        );
//
//        // Segreduce on mask
//        BlockMapGpuKernels::segreduce<BITMASK, 0, BITMASK, 2, mgpu::maximum_t<unsigned char>> <<< outputData.size(), 2*BlockT::size >>> (
//                data.toKernel(),
//                        segments.toKernel(),
//                        data.toKernel(),
//                        outputData.toKernel()
//        );
//
//        outputData.template deviceToHost<BITMASK, BLOCK>();
//
//        for (int j = 0; j < outputData.size(); ++j)
//        {
//            BlockT outBlock = outputData.template get<BLOCK>(j);
//            MaskBlockT outMask = outputData.template get<BITMASK>(j);
//
//            for (int i = 0; i < BlockT::size; ++i)
//            {
//                std::cout << outBlock[i] << ":" << (int)outMask[i] << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
//
//    BOOST_AUTO_TEST_CASE (testSolveConflicts_vector)
//    {
//        typedef float ScalarT;
//        typedef DataBlock<unsigned char, 64> MaskBlockT;
//        typedef DataBlock<ScalarT, 64> BlockT;
//        openfpm::vector_gpu<aggregate<int>> segments;
//        segments.resize(8);
//        segments.template get<0>(0) = 0;
//        segments.template get<0>(1) = 4;
//        segments.template get<0>(2) = 5;
//        segments.template get<0>(3) = 7;
//        segments.template get<0>(4) = 8;
//        segments.template get<0>(5) = 11;
//        segments.template get<0>(6) = 17;
//        segments.template get<0>(7) = 18; // Id of first non-existent data
//
//        segments.template hostToDevice<0>();
//
//        const unsigned int BITMASK = 0, BLOCK = 1;
//        BlockT block[3];
//        MaskBlockT mask;
//        for (int i = 0; i < 32; ++i)
//        {
//            block[0][i] = i + 1;
//            block[1][i] = i + 1;
//            block[2][i] = i + 1;
//            mask[i] = 1;
//        }
//        for (int i = 32; i < 64; ++i)
//        {
//            block[0][i] = 666;
//            block[1][i] = 666;
//            block[2][i] = 666;
//            mask[i] = 0;
//        }
//
//        openfpm::vector_gpu<aggregate<MaskBlockT, BlockT[3]>> data;
//        data.resize(18);
//        for (int i = 0; i < 18; ++i)
//        {
//            data.template get<BITMASK>(i) = mask;
//            for (int k=0; k<3; k++)
//            {
//                data.template get<BLOCK>(i)[k] = block[k];
//            }
//        }
//
//        data.template hostToDevice<BITMASK, BLOCK>();
//
//        // Allocate output buffer
//        openfpm::vector_gpu<aggregate<MaskBlockT, BlockT[3]>> outputData;
//        outputData.resize(segments.size()-1);
//
//        BlockMapGpuKernels::segreduce<BLOCK, 0, BITMASK, 2, mgpu::plus_t<ScalarT>> <<< outputData.size(), 2*BlockT::size >>> (
//                data.toKernel(),
//                        segments.toKernel(),
//                        data.toKernel(),
//                        outputData.toKernel()
//        );
//
//        // Segreduce on mask
//        BlockMapGpuKernels::segreduce<BITMASK, 0, BITMASK, 2, mgpu::maximum_t<unsigned char>> <<< outputData.size(), 2*BlockT::size >>> (
//                data.toKernel(),
//                        segments.toKernel(),
//                        data.toKernel(),
//                        outputData.toKernel()
//        );
//
//        outputData.template deviceToHost<BITMASK, BLOCK>();
//
//        for (int j = 0; j < outputData.size(); ++j)
//        {
//            BlockT outBlock[3];
//            for (int k = 0 ; k < 3 ;k++)
//            {outBlock[k] = outputData.template get<BLOCK>(j)[k];}
//            MaskBlockT outMask = outputData.template get<BITMASK>(j);
//
//            for (int i = 0; i < BlockT::size; ++i)
//            {
//                std::cout << "(" << outBlock[0][i] << "," << outBlock[1][i] << "," << outBlock[2][i] << ")" << ":" << (int)outMask[i] << " ";
//            }
//            std::cout << std::endl;
//        }
//    }

BOOST_AUTO_TEST_SUITE_END() // SparseGridGpu_kernels_tests

BOOST_AUTO_TEST_SUITE(BlockMapGpu_functors_tests)
BOOST_AUTO_TEST_CASE (testCompact)
{
typedef float ScalarT;
typedef DataBlock<ScalarT, 64> BlockT;
const unsigned int numPools = 5;
const unsigned int poolSize = 10;
openfpm::vector_gpu<aggregate<unsigned int>> starts;
openfpm::vector_gpu<aggregate<BlockT>> dataSrc;
openfpm::vector_gpu<aggregate<BlockT>> dataDst;

const unsigned int numEl = numPools * poolSize;
starts.resize(numPools+1);
dataSrc.resize(numEl);

// Now fill the starts...
unsigned int curStart = 0;
for (int i=0; i <= numPools; ++i)
{
	curStart += i; // Each pool contains i+1 elements
	starts.template get<0>(i) = curStart;
//            curStart += 2; // Each pool contains 2 elements
//            curStart += 5; // Each pool contains RHS elements
}
dataDst.resize(starts.template get<0>(numPools));

//        // Debug
//        std::cout << "dataDst.size() = " << dataDst.size() << std::endl;

// Fill the data
for (int i=0; i<dataSrc.size(); ++i)
{
	for (int j=0; j<BlockT::size; ++j)
	{
		dataSrc.template get<0>(i)[j] = i;
	};
}
// Copy to device
starts.hostToDevice<0>();
dataSrc.hostToDevice<0>();
// Now perform the compaction
BlockMapGpuFunctors::BlockFunctor<128>::compact(starts, poolSize, dataSrc, dataDst);

// Now retrieve the dataDst vector
dataDst.deviceToHost<0>();

// Debug output
for (int i=0; i<dataDst.size(); ++i)
{
	std::cout << "dataDst[" << i << "][0] = " << dataDst.template get<0>(i)[0] << std::endl;
}
// Validation
int counter = 0;
for (int pool=0; pool<numPools; ++pool)
{
	int numElInPool = pool+1;
	for (int i=0; i<numElInPool; ++i)
	{
		int target = pool*10 + i;
		BOOST_REQUIRE_EQUAL(dataDst.template get<0>(counter)[0], target);
		++counter;
	}
}
}

BOOST_AUTO_TEST_CASE (testReorder)
{
typedef float ScalarT;
typedef DataBlock<ScalarT, 64> BlockT;
const unsigned int numEl = 10;
const unsigned int permutationOffset = 2;
openfpm::vector_gpu<aggregate<unsigned int>> srcIndices;
openfpm::vector_gpu<aggregate<BlockT>> dataSrc;
openfpm::vector_gpu<aggregate<BlockT>> dataDst;

srcIndices.resize(numEl);
dataSrc.resize(numEl);
dataDst.resize(numEl);

// Now fill the permutation
for (int i=0; i < numEl; ++i)
{
	srcIndices.template get<0>(i) = (i+permutationOffset)%numEl; // We periodically shift the array by permutationOffset
}

// Fill the data
for (int i=0; i<dataSrc.size(); ++i)
{
	for (int j=0; j<BlockT::size; ++j)
	{
		dataSrc.template get<0>(i)[j] = i;
	};
}
// Copy to device
srcIndices.hostToDevice<0>();
dataSrc.hostToDevice<0>();
// Now perform the compaction
BlockMapGpuFunctors::BlockFunctor<128>::reorder(srcIndices, dataSrc, dataDst);

// Now retrieve the dataDst vector
dataDst.deviceToHost<0>();

// Debug output
for (int i=0; i<dataDst.size(); ++i)
{
	std::cout << "dataDst[" << i << "][0] = " << dataDst.template get<0>(i)[0] << std::endl;
}
// Validation
for (int i=0; i<numEl; ++i)
{
	int target = (i+permutationOffset)%numEl;
	BOOST_REQUIRE_EQUAL(dataDst.template get<0>(i)[0], target);
}
}

BOOST_AUTO_TEST_CASE (testSeg_reduce)
{
typedef float ScalarT;
typedef DataBlock<ScalarT, 64> BlockT;
typedef DataBlock<unsigned char, 64> MaskBlockT;
const unsigned int numSegments = 10;
const unsigned int p=0, pMask=1, pSegment=0;
openfpm::vector_gpu<aggregate<unsigned int>> segments;
openfpm::vector_gpu<aggregate<BlockT, MaskBlockT>> dataSrc;
openfpm::vector_gpu<aggregate<BlockT, MaskBlockT>> dataDst; // Mask is included here for signature issues, but never filled

segments.resize(numSegments+1);
dataDst.resize(numSegments);

// Now fill the segments...
unsigned int curStart = 0;
for (int i=0; i <= numSegments; ++i)
{
	curStart += i; // Each segment contains i+1 elements
	segments.template get<pSegment>(i) = curStart;
//            curStart += 5; // Each segment contains RHS elements
}
dataSrc.resize(segments.template get<pSegment>(numSegments));

//        // Debug
//        std::cout << "dataSrc.size() = " << dataSrc.size() << std::endl;

// Fill the data
for (int i=0; i<dataSrc.size(); ++i)
{
	for (int j=0; j<BlockT::size; ++j)
	{
		dataSrc.template get<p>(i)[j] = 1;
		dataSrc.template get<pMask>(i)[j] = 1;
	};
}
// Copy to device
segments.hostToDevice<pSegment>();
dataSrc.hostToDevice<p, pMask>();
// Now perform the compaction
typedef boost::mpl::vector<sadd_<p>> vv_reduce;
BlockMapGpuFunctors::BlockFunctor<128>::seg_reduce<pSegment, vv_reduce, boost::mpl::int_<0>>(segments, dataSrc, dataDst);

// Now retrieve the dataDst vector
dataDst.deviceToHost<0>();

// Debug output
for (int i=0; i<dataDst.size(); ++i)
{
	std::cout << "dataDst[" << i << "][0] = " << dataDst.template get<0>(i)[0] << std::endl;
}
// Validation
for (int i=0; i<dataDst.size(); ++i)
{
	BOOST_REQUIRE_EQUAL(dataDst.template get<0>(i)[0], i+1);
}
}

BOOST_AUTO_TEST_CASE (testSolve_conflicts)
{
typedef float ScalarT;
typedef DataBlock<ScalarT, 64> BlockT;
typedef DataBlock<unsigned char, 64> MaskBlockT;
const unsigned int p=0, pMask=1, pInd=0;
openfpm::vector_gpu<aggregate<unsigned int>> keys, mergeIndices, tmpIndices, keysOut;
openfpm::vector_gpu<aggregate<BlockT, MaskBlockT>> dataOld, dataNew, tmpData, dataOut;
mgpu::ofp_context_t ctx;

// Keys
keys.resize(15);
keys.template get<pInd>(0) = 0;
keys.template get<pInd>(1) = 1;
keys.template get<pInd>(2) = 2;
keys.template get<pInd>(3) = 2;
keys.template get<pInd>(4) = 3;
keys.template get<pInd>(5) = 3;
keys.template get<pInd>(6) = 3;
keys.template get<pInd>(7) = 4;
keys.template get<pInd>(8) = 5;
keys.template get<pInd>(9) = 6;
keys.template get<pInd>(10) = 7;
keys.template get<pInd>(11) = 8;
keys.template get<pInd>(12) = 9;
keys.template get<pInd>(13) = 10;
keys.template get<pInd>(14) = 10;

// Merge Indices
mergeIndices.resize(15);
mergeIndices.template get<pInd>(0) = 5;
mergeIndices.template get<pInd>(1) = 0;
mergeIndices.template get<pInd>(2) = 6;
mergeIndices.template get<pInd>(3) = 7;
mergeIndices.template get<pInd>(4) = 1;
mergeIndices.template get<pInd>(5) = 8;
mergeIndices.template get<pInd>(6) = 9;
mergeIndices.template get<pInd>(7) = 10;
mergeIndices.template get<pInd>(8) = 2;
mergeIndices.template get<pInd>(9) = 11;
mergeIndices.template get<pInd>(10) = 3;
mergeIndices.template get<pInd>(11) = 12;
mergeIndices.template get<pInd>(12) = 4;
mergeIndices.template get<pInd>(13) = 13;
mergeIndices.template get<pInd>(14) = 14;

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

// Copy to device
keys.hostToDevice<pInd>();
mergeIndices.hostToDevice<pInd>();
dataOld.hostToDevice<p, pMask>();
dataNew.hostToDevice<p, pMask>();

// Now perform the compaction
BlockMapGpuFunctors::BlockFunctor<128>::solve_conflicts<
		decltype(keys),
		decltype(dataOld),
		sadd_<p>
		>(
				keys, mergeIndices,
				dataOld, dataNew,
				tmpIndices, tmpData,
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
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(2)[0], 2);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(3)[0], 3);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(4)[0], 1);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(5)[0], 1);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(6)[0], 1);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(7)[0], 1);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(8)[0], 1);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(9)[0], 1);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(10)[0], 2);

BOOST_REQUIRE_EQUAL(dataOut.template get<p>(0)[1], 5);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(1)[1], 0);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(2)[1], 13);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(3)[1], 18);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(4)[1], 10);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(5)[1], 2);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(6)[1], 11);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(7)[1], 3);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(8)[1], 12);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(9)[1], 4);
BOOST_REQUIRE_EQUAL(dataOut.template get<p>(10)[1], 27);
}

BOOST_AUTO_TEST_SUITE_END() //SparseGridGpu_functors_tests

