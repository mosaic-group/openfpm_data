/*
 * map_vector_cuda_funcs_tests.cu
 *
 *  Created on: Aug 17, 2018
 *      Author: i-bird
 */


#define BOOST_GPU_ENABLED __host__ __device__

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "util/cuda_util.hpp"
#include "Vector/map_vector.hpp"

template<typename vector_vector_type, typename vector_out_type>
__global__ void vv_test_size(vector_vector_type vvt, vector_out_type out)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vvt.size()) return;

//	out.template get<0>(p) = vvt.template get<0>(p).size();
}

BOOST_AUTO_TEST_SUITE( vector_cuda_tests )

BOOST_AUTO_TEST_CASE ( test_vector_of_vector_gpu )
{
	typedef openfpm::vector<Box<3,float>,CudaMemory,typename memory_traits_inte<Box<3,float>>::type,memory_traits_inte> proc_boxes;

	openfpm::vector<aggregate<proc_boxes>,CudaMemory,typename memory_traits_inte<aggregate<proc_boxes>>::type,memory_traits_inte> vb_int_proc;

	vb_int_proc.resize_no_device(5);

	for (size_t i = 0 ; i< vb_int_proc.size() ; i++)
	{
		vb_int_proc.template get<0>(i).resize(5);

		for (size_t j = 0 ; j < vb_int_proc.template get<0>(i).size() ; j++)
		{
			for (size_t k = 0 ; k < 3 ; k++)
			{
				vb_int_proc.template get<0>(i).template get<0>(j)[k] = i+j;
				vb_int_proc.template get<0>(i).template get<1>(j)[k] = 100+i+j;
			}
		}

		vb_int_proc.template get<0>(i).template hostToDevice<0,1>();
	}

	vb_int_proc.template hostToDevice<0>();

	openfpm::vector_gpu<aggregate<unsigned int>> out;
	out.resize(vb_int_proc.size());

	auto ite = vb_int_proc.getGPUIterator();

	auto test = vb_int_proc.toKernel();

	std::cout << std::string(demangle(typeid(decltype(test)).name())) << std::endl;

//	vv_test_size<decltype(vb_int_proc.toKernel()),decltype(out.toKernel())><<<ite.wthr,ite.thr>>>(vb_int_proc.toKernel(),out.toKernel());
}

BOOST_AUTO_TEST_SUITE_END()
