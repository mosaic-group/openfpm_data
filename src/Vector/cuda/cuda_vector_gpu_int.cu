#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Space/Shape/Box.hpp"
#include "Vector/map_vector.hpp"
#include "NN/CellList/cuda/CellList_gpu.hpp"
#include "util/cuda/moderngpu/kernel_scan.hxx"
#define DISABLE_MPI_WRITTERS
#include "VTKWriter/VTKWriter.hpp"

template<typename vector_pos_type, typename vector_prop_type>
__global__ void test_launch(vector_pos_type set_points, vector_prop_type prop, Box<3,float> domain)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= set_points.size())
	{return;}

	float pos[3];
	pos[0] = set_points.template get<0>(p)[0];
	pos[1] = set_points.template get<0>(p)[1];
	pos[2] = set_points.template get<0>(p)[1];

	float scalar = prop.template get<0>(p);

	float v[3];

	v[0] = prop.template get<1>(p)[0];
	v[1] = prop.template get<1>(p)[1];
	v[2] = prop.template get<1>(p)[2];

	printf("Point p %f %f %f     scalar: %f  vector: %f %f %f \n",pos[0],pos[1],pos[2],scalar,v[0],v[1],v[2]);
}

template<typename grid_type>
__global__ void test_launch_grid(grid_type grid, Box<3,float> domain, ite_gpu<3> ite_gpu)
{
	GRID_ID_3(ite_gpu);

	float scalar = grid.template get<0>(key);

	float v[3];

	v[0] = grid.template get<1>(key)[0];
	v[1] = grid.template get<1>(key)[1];
	v[2] = grid.template get<1>(key)[2];

	printf("Grid point %d %d %d     scalar: %f  vector: %f %f %f \n",(int)key.get(0),(int)key.get(1),(int)key.get(2),scalar,v[0],v[1],v[2]);
}

__global__ void test_launch_cuda_native_vector(float * scalar, float * vector, int size, int capacity)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;

	if (id >= size) {return;}

	float s = scalar[id];

	float v[3];

	v[0] = vector[id + 0*capacity];
	v[1] = vector[id + 1*capacity];
	v[2] = vector[id + 2*capacity];

	printf("Vector point from CUDA %d     scalar: %f  vector: %f %f %f \n",id,s,v[0],v[1],v[2]);
}

constexpr int NN_num = 4;

template<typename celllist_type, typename vector_pos_type>
__global__ void test_launch_cell_list(celllist_type cell, vector_pos_type pos)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;

	int nn_part = 0;

	int idx = 0;
	int NN[NN_num];

	Point<3,float> xp = pos.template get<0>(id);

	auto NN_it = cell.getNNIterator(cell.getCell(xp));

	while (NN_it.isNext())
	{
		auto q = NN_it.get();

		if (idx < NN_num)
		{
			NN[idx] = q;
			idx++;
		}

		nn_part++;

		++NN_it;
	}

	printf("CELLLIST %d  nn_part: %d NN: %d %d %d %d \n",id,nn_part,NN[0],NN[1],NN[2],NN[3]);
}

BOOST_AUTO_TEST_SUITE( vector_gpu_func_verlet )

BOOST_AUTO_TEST_CASE (gpu_verlet)
{
	openfpm::vector_gpu<Point<3,float>> pos;
	openfpm::vector_gpu<aggregate<float,float[3]>> prop;


	pos.resize(100);
	prop.resize(100);

	for (size_t i = 0 ; i < 100 ; i++)
	{
		pos.template get<0>(i)[0] = (float)rand() / RAND_MAX;
		pos.template get<0>(i)[1] = (float)rand() / RAND_MAX;
		pos.template get<0>(i)[2] = (float)rand() / RAND_MAX;

		prop.template get<0>(i) = 5.0;
		prop.template get<1>(i)[0] = 3.0;
		prop.template get<1>(i)[1] = 4.0;
		prop.template get<1>(i)[2] = 5.0;
	}

	pos.template hostToDevice<0>();
	prop.template hostToDevice<0,1>();

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	auto ite = pos.getGPUIterator();

	test_launch<<<ite.wthr,ite.thr>>>(pos.toKernel(),prop.toKernel(),domain);

	//////////// Cuda interoperability

	test_launch_cuda_native_vector<<<ite.wthr,ite.thr>>>((float *)prop.template getDeviceBuffer<0>(),(float *)prop.template getDeviceBuffer<1>(),prop.size(),prop.capacity());

	//////////// Cell-list

	openfpm::vector_gpu<Point<3,float>> pos_sort;
	openfpm::vector_gpu<aggregate<float,float[3]>> prop_sort;

	pos_sort.resize(pos.size());
	prop_sort.resize(prop.size());

	size_t g_m = pos.size();

	mgpu::ofp_context_t context(false);

	const size_t (& sz)[3] = {10,10,10};

	CellList_gpu<3,float,CudaMemory,shift_only<3,float>> cl(domain,sz,2);

	cl.construct(pos,pos_sort,prop,prop_sort,context,g_m);

	test_launch_cell_list<<<ite.wthr,ite.thr>>>(cl.toKernel(),pos.toKernel());

	// scan example

	openfpm::vector_gpu<aggregate<int,int>> numbers;

	numbers.resize(100);

	for (size_t i = 0 ; i < 100 ; i++)
	{
		numbers.template get<0>(i) = (float)rand() /RAND_MAX * 10;
	}

	numbers.hostToDevice<0>();

	mgpu::scan((int *)numbers.template getDeviceBuffer<0>(), numbers.size(), (int *)numbers.template getDeviceBuffer<1>() , context);

	numbers.deviceToHost<1>();

	for (size_t i = 0 ; i < 100 ; i++)
	{
		std::cout << "Element: " << numbers.template get<0>(i) << "  scan: " << numbers.template get<1>(i) << std::endl;
	}

	//////////////// VTK

	// VTKWriter for a set of points
	VTKWriter<boost::mpl::pair<openfpm::vector_gpu<Point<3,float>>,
							   openfpm::vector_gpu<aggregate<float,float[3]>>>,
	                           VECTOR_POINTS> vtk_writer;
	vtk_writer.add(pos,prop,pos.size());

	openfpm::vector<std::string> prp_names;
	prp_names.add("scalar");
	prp_names.add("vector");

	file_type ft = file_type::ASCII;

	pos.template deviceToHost<0>();
	prop.template deviceToHost<0,1>();

	// Write the VTK file
	vtk_writer.write("particles.vtk",prp_names,"particles",ft);
}

BOOST_AUTO_TEST_CASE (gpu_m2p)
{

}

BOOST_AUTO_TEST_SUITE_END()
