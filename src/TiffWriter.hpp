#ifndef TIFFWRITER_HPP
#define TIFFWRITER_HPP

#include <iostream>
#include "tiffio.h"
#include "map_grid.hpp"
#include <string>

/*! \brief This class is able to save grid into tiff files
 *
 * This class is able to save grid into tiff files
 *
 */

template<unsigned int dim, typename T>
class TiffWriter
{
	/*! \brief Save grid into tiff files
	 *
	 * Save grid into tiff files
	 *
	 */

	template<typename grid, typename Mem> int write(grid data, std::string file)
	{
		// Grid can have several properties we can save only scalar fields
		// for each properties save one scalar fields

		// Tiff files can be saved up to 5D

		if (dim > 5)
		{
			std::cerr << "Error Tiff writer support until 5D images" << "\n";
		}

		// Open the tiff image

		uint32 width;
		uint32 height;
		TIFF *tif = TIFFOpen(file.c_str(),"w");

		// if the file is open

		if(tif)
		{
			// set width and height for 2D

			width = data.getGrid().size(0);
			height = data.getGrid().size(1);

			TIFFSetField(tif,TIFFTAG_IMAGEWIDTH, &width);
			TIFFSetField(tif,TIFFTAG_IMAGELENGTH, &height);

			// Create the tiff line, in case the grid is CPU, we have only
			// one property and is a scalar, we can directly copy the line

			typename boost::fusion::result_of::at<T::type,0>::type first_element_type;

			if (typeid(grid).name() == "grid_cpu" && T::num_prop == 1 && boost::is_array<first_element_type>::type::value == true)
			{
				// Get the grid key iterator

				grid_key_dx_iterator<dim> key = data.getIterator();

				// write all lines

				for(int i = 0; i < height ; i++)
				{
					// select the correct lines

					key.set(1,i);

					// we have only one scalar properties, get the buffer pointer
					void * buf = &data.template get<0>(key);

					TIFFWriteScanline(tif,buf,i, 0);
				}
			}
			else
			{
				// we have to create the a scan line for each properties and index array
				// each property and index array became a channel

				// Get the grid key iterator

				grid_key_dx_iterator<dim> key = data.getIterator();

				// count how many properties and how many indexes we have

				const int n_prp = total_prop<T>;

				// write all lines

				for(int i = 0; i < height ; i++)
				{
					// select the correct lines

					key.set(1,i);

					// we have only one scalar properties, get the buffer pointer
					void * buf = &data.template get<0>(key);

					TIFFWriteScanline(tif,buf,i, 0);
				}
			}
		}
	}

};

#endif
