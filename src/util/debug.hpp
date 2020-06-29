/*
 * debug.hpp
 *
 *  Created on: Mar 12, 2020
 *      Author: i-bird
 */

#ifndef DEBUG_HPP_
#define DEBUG_HPP_

#include "VTKWriter/VTKWriter.hpp"

template<typename T>
void load_sizes(T * sizes, int d, int sza)
{
	sizes[d] = sza;

	return;
}

template<typename T, typename szAType, typename ... szType>
void load_sizes(T * sizes, int d , szAType sza, szType ... sz)
{
	sizes[d] = sza;

	load_sizes(sizes,d+1,sz ...);
}

template<typename T, typename ... szType>
void print_array(std::string file,T * arr, szType ... sz)
{
	size_t sizes[sizeof...(sz)];

	int d = 0;

	load_sizes(sizes,d,sz ...);

	// ok now we create a grid

	grid_cpu<sizeof...(sz),aggregate<T>> grd(sizes);
	grd.setMemory();

	auto it = grd.getIterator();
	auto & lin = grd.getGrid();

	while (it.isNext())
	{
		auto p = it.get();

		grd.template get<0>(p) = arr[lin.LinId(p)];

		++it;
	}

	// Create a writer and write
	VTKWriter<boost::mpl::pair<grid_cpu<sizeof...(sz),aggregate<T>>,double>,VECTOR_GRIDS> vtk_g;

	Point<sizeof...(sz),double> offset;
	Box<sizeof...(sz),size_t> dom;
	Point<sizeof...(sz),double> spacing;

	for (int i = 0 ; i < sizeof...(sz) ; i++)
	{
		offset[i] = 0;
		spacing[i] = 1;
		dom.setLow(i,0);
		dom.setHigh(i,sizes[i]);
	}

	vtk_g.add(grd,offset,spacing,dom);

	openfpm::vector<std::string> prp_names;
	prp_names.add("scalar");

	vtk_g.write(file.c_str(), prp_names, "grids", file_type::ASCII);

}

template<typename T>
void print_grid(std::string file,T & grid)
{
	// Create a writer and write
	VTKWriter<boost::mpl::pair<T,double>,VECTOR_GRIDS> vtk_g;

	Point<T::dims,double> offset;
	Box<T::dims,size_t> dom;
	Point<T::dims,double> spacing;

	for (int i = 0 ; i < T::dims ; i++)
	{
		offset[i] = 0;
		spacing[i] = 1;
		dom.setLow(i,0);
		dom.setHigh(i,grid.getGrid().size(i));
	}

	vtk_g.add(grid,offset,spacing,dom);

	openfpm::vector<std::string> prp_names;

	vtk_g.write(file.c_str(), prp_names, "grids", file_type::ASCII);

}


#endif /* DEBUG_HPP_ */
