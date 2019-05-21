//
// Created by tommaso on 14/05/19.
//

#ifndef SPARSE_GRID_UNIT_TESTS_CU_
#define SPARSE_GRID_UNIT_TESTS_CU_

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "SparseGridGpu/SparseGridGpu.hpp"
#include "SparseGridGpu/SparseGridGpu_ker.cuh"
#include "SparseGridGpu/DataBlock.cuh"

#include <limits>

template<unsigned int p, typename SparseGridType, typename VectorOutType>
__global__ void copyBlocksToOutput(SparseGridType sparseGrid, VectorOutType output)
{
    int pos = blockIdx.x * blockDim.x + threadIdx.x;
    output.template get<p>(pos) = sparseGrid.template get<p>(pos);
}

template<unsigned int p, typename SparseGridType>
__global__ void insertValues(SparseGridType sparseGrid)
{
    sparseGrid.init();

//    __syncthreads();

    int pos = blockIdx.x * blockDim.x + threadIdx.x;

    if (pos == 1102)
    {printf("DTECTED POS \n");}

    sparseGrid.template insert<p>(pos) = pos;

//    __syncthreads();

    sparseGrid.flush_block_insert();
}

template<unsigned int p, unsigned int k, typename SparseGridType>
__global__ void insertValuesBlocked(SparseGridType sparseGrid)
{
    typedef BlockTypeOf<typename SparseGridType::AggregateType, p> BlockT;

    sparseGrid.init();

    __shared__ BlockT *blocks[k]; //todo: here assuming max blockDim.x = 4*BlockT::size

    int pos = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int dataBlockId = pos / BlockT::size;
    unsigned int offset = pos % BlockT::size;

    unsigned int dataBlockNum = dataBlockId % k;

    if (offset == 0) // Just one thread per data block
    {
        blocks[dataBlockNum] = sparseGrid.insertBlock<p>(dataBlockId);
        blocks[dataBlockNum]->existBitMask = ~0; // All elements are set
    }

    __syncthreads();

//    (*(blocks[dataBlockNum]))[offset] = pos;
    blocks[dataBlockNum]->block[offset] = pos;

    sparseGrid.flush_block_insert();
}

template<typename type_t, unsigned int blockLength>
struct maximum_block_t_loc  : public std::binary_function<type_t, type_t, type_t> {
    MGPU_HOST_DEVICE type_t operator()(type_t & a, type_t & b) const {
        type_t res;
        for (int i=0; i<blockLength; ++i)
        {
            auto ae = a.exist(i)?a[i]:-std::numeric_limits<typename type_t::scalarType>::max();
            auto be = b.exist(i)?b[i]:-std::numeric_limits<typename type_t::scalarType>::max();

//            res[i] = (!a.exist(i)) && (!b.exist(i)) ? a[i] : max(ae, be);
            res[i] = (!a.exist(i)) && (!b.exist(i)) ? 69 : max(ae, be); //debug
        }
        res.existBitMask = a.existBitMask | b.existBitMask;
        return res;
    }
};

template<unsigned int prp, unsigned int blockLength>
struct smax_block_loc_
{
    typedef boost::mpl::int_<prp> prop;

    template<typename red_type> using op_red = maximum_block_t_loc<red_type, blockLength>;

    template<typename red_type>
    __device__ __host__ static red_type red(red_type & r1, red_type & r2)
    {
        red_type res;
        for (int i=0; i<blockLength; ++i)
        {
            auto ae = r1.exist(i)?r1[i]:-std::numeric_limits<typename red_type::scalarType>::max();
            auto be = r2.exist(i)?r2[i]:-std::numeric_limits<typename red_type::scalarType>::max();
            auto maximum = (ae < be) ? be : ae;
            res[i] = (!r1.exist(i)) && (!r2.exist(i)) ? r1[i] : maximum;
        }
        res.existBitMask = r1.existBitMask | r2.existBitMask;
        return res;
    }

    static bool is_special()
    {
        return false;
    }

    //! is not special reduction so it does not need it
    template<typename seg_type, typename output_type>
    __device__ __host__ static void set(seg_type seg_next, seg_type seg_prev, output_type & output, int i)
    {}
};

BOOST_AUTO_TEST_SUITE(SparseGridGpu_tests)

    BOOST_AUTO_TEST_CASE(SparseGridGpu_testBackground)
    {
        typedef aggregate<DataBlock<float>> AggregateSGT;
        typedef aggregate<float> AggregateOutT;
        SparseGridGpu<AggregateSGT> sparseGrid;
        sparseGrid.template setBackground<0>(666);

        const unsigned int gridSize = 10;
        const unsigned int blockSize = 128;

        // Get output
        openfpm::vector_gpu<AggregateOutT> output;
        output.resize(gridSize * blockSize);
        copyBlocksToOutput<0> <<< gridSize, blockSize >>> (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
            std::cout << "output(" << i << ") = " << output.template get<0>(i) << std::endl;
            match &= output.template get<0>(i) == sparseGrid.template get<0>(i);
        }

        BOOST_REQUIRE_EQUAL(match, true);
    }

    BOOST_AUTO_TEST_CASE(SparseGridGpu_testInsert)
    {
        typedef aggregate<DataBlock<float>> AggregateT;
        typedef aggregate<float> AggregateOutT;
        SparseGridGpu<AggregateT> sparseGrid;
        sparseGrid.template setBackground<0>(666);

        const unsigned int gridSize = 10;
        const unsigned int blockSizeInit = 128; // Should be multiple of BlockT::size
        const unsigned int blockSizeInsert = 128;
        const unsigned int gridSizeRead = gridSize + 1;
        const unsigned int blockSizeRead = 128;

        // Prealloc insert buffer
        sparseGrid.setGPUInsertBuffer(gridSize, blockSizeInit);
//        sparseGrid.setGPUInsertBuffer(gridSize*blockSizeInsert, 1);

        // Initialize the insert buffer
//        sparseGrid.initializeGPUInsertBuffer<0>(gridSize, blockSizeInit);
        sparseGrid.initializeGPUInsertBuffer<0>();

        // Insert values
        insertValues<0> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel());
//        insertValuesBlocked<0, 2> <<< gridSize, blockSizeInsert >>> (sparseGrid.toKernel());

        // Flush inserts
        mgpu::ofp_context_t ctx;
        sparseGrid.flush<smax_block_loc_<0, BlockTypeOf<AggregateT, 0>::size>>(ctx, flush_type::FLUSH_ON_DEVICE);
//        ctx.synchronize();

        // Get output
        openfpm::vector_gpu<AggregateOutT> output;
        output.resize(gridSizeRead * blockSizeRead);

//        std::cout << demangle(typeid(decltype(sparseGrid.toKernel())).name()) << std::endl;

        copyBlocksToOutput<0> <<< gridSizeRead, blockSizeRead >>> (sparseGrid.toKernel(), output.toKernel());

        output.template deviceToHost<0>();
        sparseGrid.template deviceToHost<0>();

        // Compare
        bool match = true;
        for (size_t i = 0; i < output.size(); i++)
        {
            std::cout << "sparseGrid(" << i << ") = " << sparseGrid.template get<0>(i)
                    << " == "
                    << output.template get<0>(i) << " = output(" << i << ")"
                            << std::endl;
            match &= output.template get<0>(i) == sparseGrid.template get<0>(i);
        }

        BOOST_REQUIRE_EQUAL(match, true);
    }

BOOST_AUTO_TEST_SUITE_END()

#endif /* SPARSE_GRID_UNIT_TESTS_CU_ */
